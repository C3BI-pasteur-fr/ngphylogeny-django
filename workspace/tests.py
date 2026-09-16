from datetime import timedelta
from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.core import mail
from django.core.cache import cache
from django.test import TestCase, override_settings
from django.utils import timezone

from blast.models import BlastRun
from galaxy.models import Server
from workflows.models import Workflow
from workspace.emails import build_job_completion_email, send_job_completion_email
from workspace.models import WorkspaceHistory
from workspace.reports import (WEEKLY_TO_MONTHLY_SPAN_DAYS,
                                build_report_context, build_report_web_context,
                                gather_all_time, gather_last_7_days,
                                gather_period_totals, render_report_html)
from workspace.tasks import deleteoldgalaxyhistory, send_daily_report


class DeleteOldGalaxyHistoryTest(TestCase):
    """
    Regression tests for deleteoldgalaxyhistory(): it used to mark a
    WorkspaceHistory (and its associated Workflow) as deleted=True
    unconditionally, even when the actual Galaxy-side
    delete_history()/delete_workflow() call silently failed (both
    deletegalaxyhistory() and deletegalaxyworkflow() swallow their own
    exceptions and just log a warning). Since every future run only
    looks at deleted=False rows, a transient Galaxy failure at cleanup
    time made the row look "cleaned up" forever - the data could still
    exist on Galaxy with no way for this task to ever notice or retry.
    Only caught by code review, not by any existing test.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    def _make_old_history(self, workflow=None):
        h = WorkspaceHistory.objects.create(
            history='hist1', name='test', email='', monitored=True,
            finished=True, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server, workflow=workflow,
            history_content_json='{"some": "content"}',
            history_info_json='{"more": "content"}')
        # created_date is auto_now_add - .save()/.create() always stamp it
        # with "now", so backdate it via a queryset update (bypasses
        # auto_now_add) to simulate an old, finished history.
        WorkspaceHistory.objects.filter(pk=h.pk).update(
            created_date=timezone.now() - timedelta(days=15))
        h.refresh_from_db()
        return h

    def test_marks_deleted_and_clears_json_on_success(self):
        h = self._make_old_history()
        with patch('workspace.tasks.deletegalaxyhistory', return_value=True):
            deleteoldgalaxyhistory()
        h.refresh_from_db()
        self.assertTrue(h.deleted)
        self.assertEqual(h.history_content_json, "")
        self.assertEqual(h.history_info_json, "")

    def test_leaves_deleted_false_when_galaxy_history_delete_fails(self):
        h = self._make_old_history()
        with patch('workspace.tasks.deletegalaxyhistory', return_value=False):
            deleteoldgalaxyhistory()
        h.refresh_from_db()
        self.assertFalse(h.deleted)
        # Not touched either - still there to retry against next run.
        self.assertEqual(h.history_content_json, '{"some": "content"}')

    def test_leaves_deleted_false_when_associated_workflow_delete_fails(self):
        wf = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='wfid1',
            name='PhyML OneClick', category='duplicated',
            description='PhyML OneClick',
            slug='wfid1_PhyML OneClick_copy')
        h = self._make_old_history(workflow=wf)
        with patch('workspace.tasks.deletegalaxyworkflow',
                   return_value=False), \
             patch('workspace.tasks.deletegalaxyhistory',
                   return_value=True):
            deleteoldgalaxyhistory()
        h.refresh_from_db()
        wf.refresh_from_db()
        self.assertFalse(h.deleted)
        self.assertFalse(wf.deleted)

    def test_one_failing_history_does_not_block_the_others(self):
        h_fail = self._make_old_history()
        h_fail.history = 'hist-fail'
        h_fail.save()
        h_ok = self._make_old_history()
        h_ok.history = 'hist-ok'
        h_ok.save()
        # Re-backdate: .save() above didn't touch created_date (not
        # auto_now), so both are still old - just confirm the setup.
        h_fail.refresh_from_db()
        h_ok.refresh_from_db()

        def fake_delete_history(historyid):
            return historyid != 'hist-fail'

        with patch('workspace.tasks.deletegalaxyhistory',
                   side_effect=fake_delete_history):
            deleteoldgalaxyhistory()

        h_fail.refresh_from_db()
        h_ok.refresh_from_db()
        self.assertFalse(h_fail.deleted)
        self.assertTrue(h_ok.deleted)


class DailyReportTest(TestCase):
    """
    Tests for workspace/reports.py (data gathering + chart/HTML rendering)
    and workspace.tasks.send_daily_report, backing the daily HTML email
    report of workflow usage (CELERY_BEAT_SCHEDULE's
    'workspace-send-daily-report', 8am UTC).
    """

    def setUp(self):
        # build_report_web_context() caches its result (see
        # workspace/reports.py) - clear it so tests don't see a stale
        # value left over from a previous test.
        cache.clear()
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    _history_counter = 0

    def _make_history(self, category, workflow=None, workflow_steps='',
                       days_ago=0):
        DailyReportTest._history_counter += 1
        h = WorkspaceHistory.objects.create(
            history='h-%s-%s-%d-%d' % (
                category, workflow_steps, days_ago,
                DailyReportTest._history_counter),
            name='test', email='', monitored=True, finished=True,
            source_ip='127.0.0.1', workflow_category=category,
            workflow_steps=workflow_steps, workflow=workflow,
            galaxy_server=self.server)
        if days_ago:
            WorkspaceHistory.objects.filter(pk=h.pk).update(
                created_date=timezone.now() - timedelta(days=days_ago))
        return h

    def test_gather_last_7_days_buckets_by_day_and_category(self):
        wf = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='wf1',
            name='FastME OneClick', category='duplicated',
            description='FastME OneClick', slug='wf1-copy')
        self._make_history('OneClick', workflow=wf, days_ago=0)
        self._make_history('OneClick', workflow=wf, days_ago=0)
        self._make_history('Tool', workflow_steps='MAFFT', days_ago=3)
        # Outside the 7-day window - must not be counted here.
        self._make_history('OneClick', workflow=wf, days_ago=10)

        days, by_day_category, by_day_oneclick_workflow = gather_last_7_days()

        today = timezone.localdate()
        self.assertEqual(by_day_category[today]['OneClick'], 2)
        self.assertEqual(by_day_oneclick_workflow[today]['FastME OneClick'], 2)
        three_days_ago = today - timedelta(days=3)
        self.assertEqual(by_day_category[three_days_ago]['Tool'], 1)
        self.assertEqual(sum(sum(c.values()) for c in by_day_category.values()), 3)

    def test_gather_last_7_days_includes_blast_runs_as_a_category(self):
        """
        BlastRun isn't a WorkspaceHistory row at all (blast is a
        parallel, self-contained app - see CLAUDE.md) - gather_last_7_days()
        merges its counts into the same by_day_category dict under a
        synthetic 'blast' key so it shows up as just one more category,
        same as OneClick/Advanced/A La Carte/Single Tool.
        """
        today_run = BlastRun.objects.create(query_id='', query_seq='')
        old_run = BlastRun.objects.create(query_id='', query_seq='')
        BlastRun.objects.filter(pk=old_run.pk).update(
            date=timezone.now() - timedelta(days=10))

        _, by_day_category, _ = gather_last_7_days()

        today = timezone.localdate()
        self.assertEqual(by_day_category[today]['blast'], 1)
        # Outside the 7-day window - must not be counted here.
        self.assertEqual(
            sum(c.get('blast', 0) for c in by_day_category.values()), 1)

    def test_gather_all_time_includes_everything_regardless_of_age(self):
        self._make_history('Tool', workflow_steps='BMGE', days_ago=0)
        self._make_history('Tool', workflow_steps='BMGE', days_ago=400)
        self._make_history('automaker', days_ago=0)

        by_category, by_workflow = gather_all_time()

        self.assertEqual(by_category['Tool'], 2)
        self.assertEqual(by_category['automaker'], 1)
        self.assertEqual(by_workflow['BMGE'], 2)
        self.assertEqual(by_workflow['A La Carte'], 1)

    def test_gather_all_time_includes_blast_runs_regardless_of_deleted(self):
        """
        Same "usage report, not a what's-still-retained report" reasoning
        already applied to WorkspaceHistory (see this module's docstring)
        - a soft-deleted BlastRun (cleaned up by deleteoldblastruns()'s
        14-day cutoff) still counts as a real, historical usage.
        """
        BlastRun.objects.create(query_id='', query_seq='')
        deleted_run = BlastRun.objects.create(query_id='', query_seq='')
        deleted_run.soft_delete()

        by_category, _ = gather_all_time()

        self.assertEqual(by_category['blast'], 2)

    def test_gather_period_totals_buckets_across_iso_weeks(self):
        # Two entries on the same day land in the same week's bucket; a
        # third, 3 weeks earlier, leaves at least one fully-empty week in
        # between that must still show up as a zero, not be skipped. Well
        # under WEEKLY_TO_MONTHLY_SPAN_DAYS, so this stays weekly.
        self._make_history('Tool', workflow_steps='MAFFT', days_ago=0)
        self._make_history('Tool', workflow_steps='BMGE', days_ago=0)
        self._make_history('Tool', workflow_steps='MAFFT', days_ago=21)

        granularity, weekly_totals = gather_period_totals()

        self.assertEqual(granularity, 'week')
        self.assertEqual(sum(n for _, n in weekly_totals), 3)
        self.assertIn(0, [n for _, n in weekly_totals])
        # oldest week first, and every consecutive pair is exactly one
        # week apart - no gaps silently dropped.
        self.assertLess(weekly_totals[0][0], weekly_totals[-1][0])
        for (d1, _), (d2, _) in zip(weekly_totals, weekly_totals[1:]):
            self.assertEqual((d2 - d1).days, 7)

    def test_gather_period_totals_switches_to_monthly_for_long_spans(self):
        # A history spanning years (see CLAUDE.md's "Restoring historical
        # workspace_workspacehistory data") would otherwise produce
        # hundreds of weekly bars, crushed illegible by .report-chart's
        # max-width: 100%% once squeezed into a normal page/email width.
        self._make_history('Tool', workflow_steps='MAFFT',
                            days_ago=WEEKLY_TO_MONTHLY_SPAN_DAYS + 30)
        self._make_history('Tool', workflow_steps='BMGE', days_ago=0)

        granularity, monthly_totals = gather_period_totals()

        self.assertEqual(granularity, 'month')
        self.assertEqual(sum(n for _, n in monthly_totals), 2)
        # oldest month first, one calendar month apart between entries -
        # not a fixed 30/31-day step.
        self.assertLess(monthly_totals[0][0], monthly_totals[-1][0])
        for (d1, _), (d2, _) in zip(monthly_totals, monthly_totals[1:]):
            expected_next = (d1.replace(year=d1.year + 1, month=1)
                              if d1.month == 12
                              else d1.replace(month=d1.month + 1))
            self.assertEqual(d2, expected_next)

    def test_single_tool_runs_use_workflow_steps_not_a_workflow_fk(self):
        """
        Regression guard: 'Tool' category WorkspaceHistory rows never get
        a Workflow FK (see tools/views.py's create_history() calls) - only
        workflow_steps carries the tool name. If _workflow_label ever
        started preferring workflow__name unconditionally, every
        single-tool run would collapse into a single 'Unknown' bucket.
        """
        self._make_history('Tool', workflow=None, workflow_steps='Gblocks')
        _, by_workflow = gather_all_time()
        self.assertEqual(by_workflow['Gblocks'], 1)
        self.assertNotIn('Unknown', by_workflow)

    def test_build_report_context_renders_charts_when_data_present(self):
        """
        Charts are inline (Content-ID) attachments, not base64 data: URIs
        - see workspace/reports.py's module docstring for why (many mail
          clients, Outlook included, don't render data: URI images in
          HTML email at all). build_report_context() returns the cid:
          name a chart was rendered under in the context, and the actual
          PNG bytes in a separate images dict for the caller to attach.
        """
        wf = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='wf2',
            name='PhyML OneClick', category='duplicated',
            description='PhyML OneClick', slug='wf2-copy')
        self._make_history('OneClick', workflow=wf)

        context, images = build_report_context()

        self.assertEqual(context['alltime_total'], 1)
        self.assertEqual(context['week_total'], 1)
        expected_charts = ['daily_category_chart', 'daily_oneclick_chart',
                            'alltime_category_chart', 'alltime_workflow_chart',
                            'weekly_chart']
        for key in expected_charts:
            cid = context[key]
            self.assertIsNotNone(cid, key)
            self.assertIn(cid, images)
            self.assertGreater(len(images[cid]), 0)
        # PNG magic bytes - these really are images, not placeholders.
        for png_bytes in images.values():
            self.assertTrue(png_bytes.startswith(b'\x89PNG\r\n\x1a\n'))
        self.assertEqual(len(images), len(expected_charts))

    def test_build_report_context_handles_no_data_at_all(self):
        context, images = build_report_context()
        self.assertEqual(context['alltime_total'], 0)
        self.assertEqual(context['week_total'], 0)
        self.assertIsNone(context['daily_category_chart'])
        self.assertIsNone(context['daily_oneclick_chart'])
        self.assertIsNone(context['alltime_category_chart'])
        self.assertIsNone(context['alltime_workflow_chart'])
        self.assertIsNone(context['weekly_chart'])
        self.assertEqual(images, {})
        # Must still render without error - the template has to handle
        # every chart being None gracefully.
        html = render_report_html(context)
        self.assertIn('NGPhylogeny.fr', html)

    @override_settings(NGPHYLO_REPORT_RECIPIENTS=[])
    def test_send_daily_report_noops_without_recipients(self):
        send_daily_report()
        self.assertEqual(len(mail.outbox), 0)

    @override_settings(NGPHYLO_REPORT_RECIPIENTS=['team@example.org'],
                        NGPHYLO_REPORT_FROM_EMAIL='ngphylo@example.org')
    def test_send_daily_report_sends_html_email_to_configured_recipients(self):
        self._make_history('Tool', workflow_steps='MAFFT')

        send_daily_report()

        self.assertEqual(len(mail.outbox), 1)
        sent = mail.outbox[0]
        self.assertEqual(sent.to, ['team@example.org'])
        self.assertEqual(sent.from_email, 'ngphylo@example.org')
        self.assertEqual(len(sent.alternatives), 1)
        html_body, mimetype = sent.alternatives[0]
        self.assertEqual(mimetype, 'text/html')
        self.assertIn('MAFFT', html_body)
        # multipart/related, not multipart/mixed - and the HTML actually
        # references the charts as cid:, not as data: URIs.
        self.assertEqual(sent.mixed_subtype, 'related')
        self.assertNotIn('data:image', html_body)
        self.assertIn('cid:chart_alltime_category', html_body)
        # The charts are attached as inline images with matching
        # Content-IDs, not as ordinary (non-inline) file attachments.
        self.assertGreater(len(sent.attachments), 0)
        for attachment in sent.attachments:
            self.assertEqual(attachment.get_content_type(), 'image/png')
            self.assertEqual(attachment['Content-Disposition'].split(';')[0],
                              'inline')
            content_id = attachment['Content-ID'].strip('<>')
            self.assertIn('cid:%s' % content_id, html_body)

    def test_build_report_web_context_uses_data_uris_not_cid(self):
        """
        Unlike the emailed report, the web page is rendered directly in a
        browser, which has no trouble with data: URI images (it's only
        mail clients like Outlook that don't render them - see
        reports.py's module docstring) - so it embeds charts that way
        instead of needing a separate cid:-matched attachment mechanism.
        """
        self._make_history('Tool', workflow_steps='MAFFT')
        context = build_report_web_context()
        self.assertTrue(context['alltime_category_chart'].startswith(
            'data:image/png;base64,'))
        self.assertIsNone(context['daily_oneclick_chart'])

    def test_build_report_web_context_is_cached(self):
        """
        Rendering 5 matplotlib charts on every single page view was the
        actual slow part a real user hit - build_report_web_context()
        caches its result rather than rebuilding it from scratch on every
        call within REPORT_WEB_CACHE_TTL.
        """
        self._make_history('Tool', workflow_steps='MAFFT')

        with patch('workspace.reports.build_report_context',
                   wraps=build_report_context) as mock_build:
            first = build_report_web_context()
            second = build_report_web_context()
            self.assertEqual(mock_build.call_count, 1)
        self.assertEqual(first, second)

        with patch('workspace.reports.build_report_context',
                   wraps=build_report_context) as mock_build:
            build_report_web_context(force_refresh=True)
            self.assertEqual(mock_build.call_count, 1)


class DailyReportViewTest(TestCase):
    """
    Access-control tests for workspace.views.daily_report_view (URL name
    'daily_report', /workspace/report) - must be admin/staff only.
    """

    def setUp(self):
        cache.clear()
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    def test_anonymous_user_is_redirected_to_login(self):
        response = self.client.get('/workspace/report')
        self.assertEqual(response.status_code, 302)
        self.assertIn('/admin/login/', response.url)

    def test_non_staff_user_is_redirected_to_login(self):
        User.objects.create_user('regularuser', password='pw')
        self.client.login(username='regularuser', password='pw')
        response = self.client.get('/workspace/report')
        self.assertEqual(response.status_code, 302)
        self.assertIn('/admin/login/', response.url)

    def test_staff_user_sees_the_report(self):
        User.objects.create_user('staffuser', password='pw', is_staff=True)
        self.client.login(username='staffuser', password='pw')
        WorkspaceHistory.objects.create(
            history='view-test', name='test', email='', monitored=True,
            finished=True, source_ip='127.0.0.1', workflow_category='Tool',
            workflow_steps='MAFFT', galaxy_server=self.server)

        response = self.client.get('/workspace/report')

        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        self.assertIn('Daily Workflow Report', content)
        self.assertIn('MAFFT', content)
        self.assertIn('data:image/png;base64,', content)

    def test_refresh_param_bypasses_the_cache(self):
        User.objects.create_user('staffuser', password='pw', is_staff=True)
        self.client.login(username='staffuser', password='pw')

        with patch('workspace.views.build_report_web_context',
                   wraps=build_report_web_context) as mock_build:
            self.client.get('/workspace/report')
            self.client.get('/workspace/report')
            self.assertEqual(
                [c.kwargs.get('force_refresh', False)
                 for c in mock_build.call_args_list],
                [False, False])

            self.client.get('/workspace/report?refresh=1')
            self.assertTrue(
                mock_build.call_args_list[-1].kwargs.get('force_refresh'))


class JobCompletionEmailTest(TestCase):
    """
    Tests for workspace/emails.py (the HTML job-completion email sent by
    workspace.tasks.updateworkspacestatus once a monitored history's jobs
    are all done) - NGPhylogeny.fr/Institut Pasteur branded, with inline
    (Content-ID) header/footer logos rather than base64 data: URIs (see
    the module docstring - same reasoning as the daily report).
    """

    def _built(self, error=False):
        return build_job_completion_email('hist123', 'user@example.org', error)

    @override_settings(NGPHYLO_REPORT_FROM_EMAIL='ngphylogeny@pasteur.fr')
    def test_success_email_content_and_structure(self):
        msg = self._built(error=False)

        self.assertEqual(msg.to, ['user@example.org'])
        self.assertEqual(msg.from_email, 'ngphylogeny@pasteur.fr')
        self.assertIn('finished', msg.subject)
        self.assertNotIn('error', msg.subject.lower())

        self.assertEqual(len(msg.alternatives), 1)
        html_body, mimetype = msg.alternatives[0]
        self.assertEqual(mimetype, 'text/html')
        self.assertIn('Institut Pasteur', html_body)
        self.assertIn('Finished successfully', html_body)
        self.assertNotIn('Finished with errors', html_body)
        self.assertIn('hist123', html_body)
        self.assertIn('doi.org/10.1093/nar/gkz303', html_body)
        # No C3BI/CNRS branding - just Institut Pasteur.
        self.assertNotIn('C3BI', html_body)
        self.assertNotIn('CNRS', html_body)

        # multipart/related with two inline logo attachments, matching
        # cid: references in the HTML - not data: URIs, and not ordinary
        # (non-inline) file attachments - see reports.py/tasks.py's
        # daily-report email for why this exact structure matters.
        self.assertEqual(msg.mixed_subtype, 'related')
        self.assertNotIn('data:image', html_body)
        self.assertEqual(len(msg.attachments), 2)
        for attachment in msg.attachments:
            self.assertEqual(attachment.get_content_type(), 'image/png')
            self.assertEqual(attachment['Content-Disposition'].split(';')[0],
                              'inline')
            content_id = attachment['Content-ID'].strip('<>')
            self.assertIn('cid:%s' % content_id, html_body)

    def test_error_email_shows_error_status(self):
        msg = self._built(error=True)
        self.assertIn('error', msg.subject.lower())
        html_body, _ = msg.alternatives[0]
        self.assertIn('Finished with errors', html_body)
        self.assertNotIn('Finished successfully', html_body)

    @override_settings(NGPHYLO_HTTPS_HOST='ngphylogeny.fr')
    def test_results_link_uses_https_when_configured(self):
        msg = self._built()
        html_body, _ = msg.alternatives[0]
        self.assertIn('https://ngphylogeny.fr/workspace/history/hist123',
                       html_body)

    @override_settings(NGPHYLO_HTTPS_HOST=None)
    def test_results_link_falls_back_to_http(self):
        msg = self._built()
        html_body, _ = msg.alternatives[0]
        # Only the site links should fall back to http - the citation's
        # DOI link is always https regardless of site config.
        self.assertIn('http://ngphylogeny.fr/workspace/history/hist123',
                       html_body)
        self.assertNotIn('https://ngphylogeny.fr', html_body)

    def test_send_job_completion_email_actually_sends(self):
        send_job_completion_email('hist456', 'someone@example.org', False)
        self.assertEqual(len(mail.outbox), 1)
        self.assertEqual(mail.outbox[0].to, ['someone@example.org'])

    @override_settings(NGPHYLO_REPORT_FROM_EMAIL='authorized@pasteur.fr')
    def test_sender_reuses_the_report_from_email_setting(self):
        """
        Regression test: this used to hardcode 'ngphylogeny@pasteur.fr' as
        the sender, which the real SMTP account isn't authorized to send
        as (institutional "Send As" restriction) - hit for real sending a
        live test email, exactly like the daily report hit the same class
        of issue earlier. Reuses NGPHYLO_REPORT_FROM_EMAIL instead of a
        second, independently-hardcoded address.
        """
        msg = self._built()
        self.assertEqual(msg.from_email, 'authorized@pasteur.fr')
