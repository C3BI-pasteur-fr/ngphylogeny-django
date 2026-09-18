import json
from datetime import timedelta
from unittest.mock import Mock, patch

import requests
from django.contrib.auth.models import User
from django.core import mail, signing
from django.core.cache import cache
from django.template.loader import render_to_string
from django.test import RequestFactory, TestCase, override_settings
from django.urls import reverse

from django.utils import timezone

from blast.models import BlastRun
from galaxy.models import GalaxyUser, Server
from tools.models import Tool
from workflows.models import Workflow
from workspace.emails import build_job_completion_email, send_job_completion_email
from workspace.models import WorkspaceHistory
from workspace.reports import (WEEKLY_TO_MONTHLY_SPAN_DAYS,
                                build_report_context, build_report_web_context,
                                gather_all_time, gather_last_7_days,
                                gather_period_totals, render_report_html)
from workspace.tasks import (deleteoldgalaxyhistory, send_daily_report,
                              updateworkspacestatus)
from workspace.views import PERMALINK_SALT, build_citations, resolve_dataset_tools


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


class UpdateWorkspaceStatusStaleRunTest(TestCase):
    """
    Regression tests for updateworkspacestatus()'s new staleness check:
    a WorkspaceHistory used to stay finished=False forever if its Galaxy
    jobs got stuck (a hung cluster node, a tool that never returns) -
    launchmonitorworkspaces just keeps polling it, and
    deleteoldgalaxyhistory only ever looks at finished=True rows, so
    nothing ever cleaned it up either. Same class of gap
    blast.tests.CheckBlastRunsTest already covers for BLAST
    (PASTEUR_RUN_STALE_AFTER) - this is that same fix for regular
    workflow/tool runs (WORKFLOW_RUN_STALE_AFTER, 24h).

    Goes through galaxy_connection() for real (not request.galaxy/
    connection_galaxy - this is a Celery task, not a view) - same
    Server + anonymous GalaxyUser DB fixture established elsewhere this
    session (e.g. tools.tests.GetToolNameViewTest), with the specific
    bioblend client methods actually called patched directly.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)

    def _make_history(self, age, history_content):
        h = WorkspaceHistory.objects.create(
            history='hist1', name='test', email='', monitored=True,
            finished=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server,
            history_content_json=json.dumps(history_content),
            history_info_json=json.dumps({'id': 'hist1', 'name': 'test'}))
        WorkspaceHistory.objects.filter(pk=h.pk).update(
            created_date=timezone.now() - age)
        h.refresh_from_db()
        return h

    def _fake_show_history(self, history_content):
        def fake(history_id, contents=False, **kwargs):
            if contents:
                return history_content
            return {'id': history_id, 'name': 'test'}
        return fake

    def test_cancels_stale_jobs_and_marks_running_datasets_as_error(self):
        history_content = [
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'still going', 'state': 'running',
             'visible': True, 'extension': 'fasta'},
        ]
        h = self._make_history(timedelta(hours=25), history_content)
        jobs = [
            {'id': 'job-finished', 'state': 'ok'},
            {'id': 'job-stuck', 'state': 'running'},
        ]
        with patch('bioblend.galaxy.histories.HistoryClient.show_history',
                   side_effect=self._fake_show_history(history_content)), \
             patch('bioblend.galaxy.jobs.JobsClient.get_jobs',
                   return_value=jobs) as get_jobs, \
             patch('bioblend.galaxy.jobs.JobsClient.cancel_job') as cancel_job:
            updateworkspacestatus(h.history)

        get_jobs.assert_called_once_with(history_id='hist1')
        cancel_job.assert_called_once_with('job-stuck')

        h.refresh_from_db()
        self.assertTrue(h.finished)
        saved_content = json.loads(h.history_content_json)
        by_id = {f['id']: f['state'] for f in saved_content}
        self.assertEqual(by_id['d1'], 'ok')  # untouched - already done
        self.assertEqual(by_id['d2'], 'error')  # was running, now error

    def test_leaves_a_fresh_still_running_history_alone(self):
        history_content = [
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'still going', 'state': 'running',
             'visible': True, 'extension': 'fasta'},
        ]
        h = self._make_history(timedelta(hours=1), history_content)
        with patch('bioblend.galaxy.histories.HistoryClient.show_history',
                   side_effect=self._fake_show_history(history_content)), \
             patch('bioblend.galaxy.jobs.JobsClient.get_jobs') as get_jobs, \
             patch('bioblend.galaxy.jobs.JobsClient.cancel_job') as cancel_job:
            updateworkspacestatus(h.history)

        get_jobs.assert_not_called()
        cancel_job.assert_not_called()
        h.refresh_from_db()
        self.assertFalse(h.finished)


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

    def test_blast_query_length_histogram_renders_when_data_present(self):
        """
        BlastRun.query_length (blast/models.py) is set at submission time
        and re-derived at cleanup time if it was ever missed (see
        blast/tasks.py) - build_report_context() renders it as an
        all-time histogram, separate from the day/category breakdowns
        above since it isn't bucketed by day or category at all. Rows
        with query_length still NULL (not yet set/re-derived) are
        excluded rather than counted as a zero-length search.
        """
        BlastRun.objects.create(query_id='', query_seq='', query_length=120)
        BlastRun.objects.create(query_id='', query_seq='', query_length=350)
        BlastRun.objects.create(query_id='', query_seq='', query_length=None)

        context, images = build_report_context()

        self.assertEqual(context['blast_length_count'], 2)
        cid = context['blast_length_chart']
        self.assertIsNotNone(cid)
        self.assertIn(cid, images)
        self.assertTrue(images[cid].startswith(b'\x89PNG\r\n\x1a\n'))

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
        self.assertIsNone(context['blast_length_chart'])
        self.assertEqual(context['blast_length_count'], 0)
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


class RunningJobsViewTest(TestCase):
    """
    Access-control + content tests for workspace.views.running_jobs_view
    (URL name 'running_jobs', /workspace/running) - admin/staff only,
    same access-control pattern as daily_report_view
    (DailyReportViewTest above). Lists still-running WorkspaceHistory and
    BlastRun rows together, oldest first - finished/deleted rows of
    either kind must not show up, since this is meant to answer "what's
    running right now", not a usage history.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    def test_anonymous_user_is_redirected_to_login(self):
        response = self.client.get('/workspace/running')
        self.assertEqual(response.status_code, 302)
        self.assertIn('/admin/login/', response.url)

    def test_non_staff_user_is_redirected_to_login(self):
        User.objects.create_user('regularuser', password='pw')
        self.client.login(username='regularuser', password='pw')
        response = self.client.get('/workspace/running')
        self.assertEqual(response.status_code, 302)
        self.assertIn('/admin/login/', response.url)

    def test_staff_user_sees_running_jobs_oldest_first_with_step_counts(self):
        User.objects.create_user('staffuser', password='pw', is_staff=True)
        self.client.login(username='staffuser', password='pw')

        older = WorkspaceHistory.objects.create(
            history='hist-older', name='Older run', email='', monitored=True,
            finished=False, deleted=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server,
            history_content_json=json.dumps([
                {'id': 'd1', 'hid': 1, 'name': 'a', 'state': 'ok'},
                {'id': 'd2', 'hid': 2, 'name': 'b', 'state': 'running'},
            ]))
        WorkspaceHistory.objects.filter(pk=older.pk).update(
            created_date=timezone.now() - timedelta(hours=2))

        WorkspaceHistory.objects.create(
            history='hist-newer', name='Newer run', email='', monitored=True,
            finished=False, deleted=False, source_ip='127.0.0.1',
            workflow_category='Tool', workflow_steps='MAFFT',
            galaxy_server=self.server,
            history_content_json=json.dumps([
                {'id': 'd1', 'hid': 1, 'name': 'a', 'state': 'queued'},
            ]))

        # Should NOT show up: finished, or deleted.
        WorkspaceHistory.objects.create(
            history='hist-finished', name='Finished run', email='',
            monitored=True, finished=True, deleted=False,
            source_ip='127.0.0.1', workflow_category='OneClick',
            workflow_steps='', galaxy_server=self.server)
        WorkspaceHistory.objects.create(
            history='hist-deleted', name='Deleted run', email='',
            monitored=True, finished=False, deleted=True,
            source_ip='127.0.0.1', workflow_category='OneClick',
            workflow_steps='', galaxy_server=self.server)

        blast_running = BlastRun.objects.create(
            query_id='query1', status=BlastRun.RUNNING, server=BlastRun.NCBI,
            blastprog='blastn', deleted=False)
        BlastRun.objects.filter(pk=blast_running.pk).update(
            date=timezone.now() - timedelta(hours=1))
        # Should NOT show up: a finished BLAST run.
        BlastRun.objects.create(
            query_id='query2', status=BlastRun.FINISHED, server=BlastRun.NCBI,
            blastprog='blastn', deleted=False)

        response = self.client.get('/workspace/running')

        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        self.assertIn('Older run', content)
        self.assertIn('Newer run', content)
        self.assertIn('query1', content)
        # Regression guard: BlastRun.BLASTSERVERS used to label NCBI
        # runs as "Pasteur" too (see blast/models.py) - this page's own
        # server_labels lookup would have silently inherited that.
        self.assertIn('NCBI BLAST', content)
        self.assertNotIn('Finished run', content)
        self.assertNotIn('Deleted run', content)
        self.assertNotIn('query2', content)
        # Oldest (2h ago) first, then the BLAST run (1h ago), then the
        # just-created workflow run.
        self.assertLess(content.index('Older run'), content.index('query1'))
        self.assertLess(content.index('query1'), content.index('Newer run'))
        # 1 of "Older run"'s 2 datasets is 'ok'.
        self.assertIn('1 / 2', content)


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


class HistoryStepChainTemplateTest(TestCase):
    """
    templates/workspace/include/history_contents_refreshable.html's
    graphical step-chain (above the detailed dataset table, included by
    history_contents_provenance_ajax.html into workspace/history.html) -
    one box per tool step, colored by status. Rendered directly against
    the template (not through HistoryDetailView) since Galaxy-connected
    views (galaxy.decorator.connection_galaxy) have no existing
    mocking pattern anywhere in this codebase to build a real
    request/response test on - this covers what's actually new here
    (the template logic/JS it emits), the same way
    DailyReportTest.test_build_report_context_handles_no_data_at_all
    above calls render_report_html() directly rather than through a
    view.
    """

    def _render(self, history_content):
        rf = RequestFactory()
        request = rf.get('/workspace/history/fakehist123')
        request.session = {}

        class FakeObj:
            pass
        obj = FakeObj()
        obj.history_content = history_content
        obj.history_info = {'id': 'fakehist123', 'name': 'Test run'}
        obj.finished = False
        obj.name = 'Test run'
        obj.email = ''
        obj.workflow = None

        return render_to_string('workspace/history.html',
                                 {'object': obj, 'request': request,
                                  'csrf_token': 'faketoken'},
                                 request=request)

    def test_step_chain_container_present_with_multiple_datasets(self):
        html = self._render([
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'MAFFT alignment', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
        ])
        self.assertIn('id="workflow-step-chain"', html)
        self.assertIn('var historySteps', html)

    def test_history_steps_emitted_oldest_first_with_id_name_state(self):
        """
        The table below sorts newest-first (dictsortreversed) - the step
        chain deliberately uses the opposite order (dictsort) since a
        left-to-right chain of boxes reads naturally as progress through
        the pipeline, from its first step to its most recent.
        """
        html = self._render([
            {'id': 'd2', 'hid': 2, 'name': 'MAFFT alignment', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
        ])
        start = html.index('var historySteps')
        end = html.index('];', start)
        snippet = html[start:end]
        self.assertLess(
            snippet.index('input.fasta'), snippet.index('MAFFT alignment'),
            'expected the oldest dataset (hid 1) before the newest (hid 2)')

    def test_history_steps_are_escaped(self):
        """
        Regression guard: dataset names come from Galaxy, not from a
        trusted source - a name containing a quote or HTML must not
        break out of the JS string literal.
        """
        html = self._render([
            {'id': 'd1', 'hid': 1, 'name': 'weird"</script><b>name',
             'state': 'ok', 'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'second', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
        ])
        self.assertNotIn('weird"</script>', html)
        self.assertIn('</script>', html)  # the real closing tags survive


class HistoryTableRedesignTest(TestCase):
    """
    Regression tests for the dataset table's redesign (card container,
    color-coded status pills reusing the step-chain's own palette, line-
    art SVG icon buttons with a fixed set of action slots so the same
    action lines up in the same column across rows) - the old plain
    Bootstrap-striped table with bare glyphicon buttons should be
    entirely gone, not just visually superseded.
    """

    def _render(self, history_content):
        rf = RequestFactory()
        request = rf.get('/workspace/history/fakehist123')
        request.session = {}

        class FakeObj:
            pass
        obj = FakeObj()
        obj.history_content = history_content
        obj.history_info = {'id': 'fakehist123', 'name': 'Test run'}
        obj.finished = False
        obj.name = 'Test run'
        obj.email = ''
        obj.workflow = None

        return render_to_string('workspace/history.html',
                                 {'object': obj, 'request': request,
                                  'csrf_token': 'faketoken'},
                                 request=request)

    def test_ok_row_shows_a_done_pill_and_no_old_glyphicons(self):
        html = self._render([
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'MAFFT alignment', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
        ])
        self.assertIn('history-table-card', html)
        self.assertIn('status-pill-ok', html)
        self.assertIn('<span class="dot"></span>Done', html)
        # The whole point: no leftover glyphicon-based buttons/status
        # from the previous design.
        self.assertNotIn('glyphicon-ok', html)
        self.assertNotIn('glyphicon-download-alt', html)
        self.assertNotIn('glyphicon-eye-open', html)

    def test_action_slots_stay_aligned_for_a_plain_dataset(self):
        """
        A dataset with no special extension gets the 4 always-present
        actions (session/params/stdout/download/display is 5, actually -
        see below) plus two empty spacer slots (shared tree/MSA viewer
        slot, iTOL slot) so it still lines up with rows that do have
        those actions.
        """
        html = self._render([
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'some.txt', 'state': 'ok',
             'visible': True, 'extension': 'txt'},
        ])
        # d2 is .txt - no MSAViewer/tree-viewer/iTOL - just the two
        # trailing empty spacer slots.
        self.assertIn(
            '<span class="action-slot"></span>\n              '
            '<span class="action-slot wide"></span>', html)

    def test_error_row_shows_an_error_pill_and_the_messages_button(self):
        html = self._render([
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'broken output', 'state': 'error',
             'visible': True, 'extension': 'txt'},
        ])
        self.assertIn('status-pill-error', html)
        self.assertIn('<span class="dot"></span>Error', html)
        self.assertIn('icon-btn danger" title="Show messages"', html)


class HistoryPartialRefreshTemplateTest(TestCase):
    """
    Regression test: the history detail page used to poll with a full
    location.reload() every 10s, then (still too slow-feeling) a plain
    $('#history-refreshable-region').load(...) every 10s - .load() blanks
    the region immediately, so the step chain/table visibly disappeared
    and slowly rebuilt (several AJAX round trips for tool-name
    resolution) on every tick. It now loads into a hidden
    #history-refreshable-staging container instead (?staging=1) and lets
    history_contents_refreshable.html's own script build the whole thing
    out of sight, only swapping the finished result into
    #history-refreshable-region once fully resolved. Covers: the outer
    page wiring up window.historyRefreshTimer/the staging container the
    JS targets; a non-staging render (the very first, direct {% include
    %}) targeting the live region directly and never containing the swap
    script (running it there would empty() the live region using the
    otherwise-untouched, empty staging container - see the template's own
    comments); a staging render targeting the staging container and
    containing the swap; and the refreshed fragment's own "the run is
    finished" check knowing to stop polling - re-evaluated fresh on every
    reload rather than only once at the initial page load, unlike before.
    """

    def _obj(self, finished):
        class FakeObj:
            pass
        obj = FakeObj()
        # 2+ datasets - the "please wait" state (see
        # HistoryUnifiedWaitStateTest below) is what history_contents_
        # refreshable.html shows for fewer than that.
        obj.history_content = [
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'MAFFT alignment', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
        ]
        obj.history_info = {'id': 'fakehist123', 'name': 'Test run'}
        obj.finished = finished
        obj.name = 'Test run'
        obj.email = ''
        obj.workflow = None
        return obj

    def _request(self):
        request = RequestFactory().get('/workspace/history/fakehist123')
        request.session = {}
        return request

    def test_outer_page_wires_up_the_refreshable_region_and_timer(self):
        request = self._request()
        html = render_to_string(
            'workspace/history.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertIn('id="history-refreshable-region"', html)
        self.assertIn('id="history-refreshable-staging"', html)
        self.assertIn('window.historyRefreshTimer', html)
        self.assertIn(
            "$('#history-refreshable-staging').load(", html)
        # function refresh() { location.reload(); } is what this
        # replaced - checked for the actual function body, not just the
        # bare phrase "location.reload()", since an explanatory comment
        # in the new code legitimately mentions it too.
        self.assertNotIn('function refresh() {\n\tlocation.reload();', html)

    def test_non_staging_fragment_targets_the_live_region_with_no_swap_script(self):
        request = self._request()
        html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertIn("var $root = $('#history-refreshable-region');", html)
        self.assertNotIn(
            '$.when(window.__historyStepChainReady, '
            'window.__historyTableReady).done(function () {', html)

    def test_staging_fragment_targets_the_staging_container_and_has_the_swap(self):
        request = self._request()
        html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken', 'staging': True},
            request=request)
        self.assertIn(
            "var $root = $('#history-refreshable-staging');", html)
        self.assertIn(
            '$.when(window.__historyStepChainReady, '
            'window.__historyTableReady).done(function () {', html)
        self.assertIn(
            "$('#history-refreshable-region').empty()."
            "append($root.contents());", html)

    def test_unfinished_fragment_does_not_stop_polling(self):
        request = self._request()
        html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertNotIn('clearInterval(window.historyRefreshTimer)', html)

    def test_finished_fragment_stops_polling(self):
        request = self._request()
        html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj(True), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertIn('clearInterval(window.historyRefreshTimer)', html)

    def test_info_refresh_banner_is_gone(self):
        """
        The "This page is will be refreshed in N sec." banner + countdown
        (#info-refresh) were removed outright - they described the old
        full-page-reload behavior and were never updated for the
        background-staged refresh above, which doesn't visibly reload
        anything for the banner to announce.
        """
        request = self._request()
        html = render_to_string(
            'workspace/history.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertNotIn('info-refresh', html)
        self.assertNotIn('countdown_span', html)

    def test_manual_refresh_button_is_gone(self):
        """
        The manual "Refresh" button (onclick="refresh()") is redundant
        now that the page already refreshes itself in the background
        every 10s - removed, though the refresh() function itself stays
        (still what drives that automatic polling).
        """
        request = self._request()
        html = render_to_string(
            'workspace/history.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertNotIn('onclick="refresh()"', html)
        self.assertIn('function refresh() {', html)

    def test_table_grouping_never_merges_unresolved_rows(self):
        """
        Regression test: the table's grouping script starts id_tool at ''
        and an unresolved dataset's own cell stays at its default data=''
        too (never overwritten - there's nothing to set) - a bare
        `id_tool == toolId` comparison matches '' against '' from the
        very first row onward whenever nothing resolves, collapsing every
        such row into one and $(el).remove()-ing the rest instead of
        leaving each its own visible row with a blank Tool label. Hit for
        real: a bug in the load-reduction refactor above left
        dataset_tool_ids/tool_names completely unresolved on
        HistoryDetailView's first render, and this compounded it further
        - not just blank Tool labels, but most rows vanishing from the
        table outright.
        """
        request = self._request()
        html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken', 'staging': True},
            request=request)
        self.assertIn("toolId !== '' && id_tool === toolId", html)

    def test_citations_are_embedded_server_side_not_fetched_by_ajax(self):
        """
        Regression test: citations used to be fetched via a client-side
        $.getJSON call from within this fragment (itself a fix for an
        even older bug - see git history - where it was only ever
        fetched once from the outer shell, missing anything from a page
        loaded during the "please wait" state). That AJAX call was one
        of three independent, redundant Galaxy-backed AJAX round trips
        this fragment made every 10s poll (see CLAUDE.md's step-chain
        section for the real production incident - pod restarts - this
        contributed to) - now HistoryContentRefreshView.get_context_data
        resolves citations server-side (reusing the same dataset_tool_ids
        it already computes for the step chain/table) and this fragment
        just writes the result into the DOM, no AJAX call left at all.
        """
        request = self._request()
        fragment_html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj(False), 'request': request,
             'csrf_token': 'faketoken', 'staging': True,
             'citations': ['<b>Citation A</b>', '<b>Citation B</b>']},
            request=request)
        self.assertNotIn('$.getJSON(', fragment_html)
        # escapejs turns <b> into <b> (defense in depth against
        # a citation string breaking out of the JS string literal) -
        # check for the actual escaped form, not the raw text.
        self.assertIn('Citation A\\u003C/b\\u003E', fragment_html)
        self.assertIn('Citation B\\u003C/b\\u003E', fragment_html)
        self.assertIn('$("#pub-container").html(citationItems', fragment_html)


class HistoryUnifiedWaitStateTest(TestCase):
    """
    Regression test: workspace/include/history_wait.html used to be a
    completely separate template (full-page location.reload() every 10s,
    its own differently-styled name/email inputs, a Url field, its own
    copy of the info-refresh banner) shown instead of
    history_contents_provenance_ajax.html whenever a history had fewer
    than 2 datasets yet (a run just starting). Deleted - that state is
    now just a branch inside history_contents_refreshable.html (see its
    own top comment), so it automatically gets the exact same shell
    (styled name/email panel, no Url field, no info-refresh banner) and
    the exact same background-staged 10s polling as the real step-chain/
    table view, with nothing bespoke left to drift out of sync again.
    """

    def _obj(self, history_content):
        class FakeObj:
            pass
        obj = FakeObj()
        obj.history_content = history_content
        obj.history_info = {'id': 'fakehist123', 'name': 'Test run'}
        obj.finished = False
        obj.name = 'Test run'
        obj.email = ''
        obj.workflow = None
        return obj

    def _request(self):
        request = RequestFactory().get('/workspace/history/fakehist123')
        request.session = {}
        return request

    def test_no_datasets_yet_shows_the_wait_message_via_the_same_shell(self):
        request = self._request()
        html = render_to_string(
            'workspace/history.html',
            {'object': self._obj([]), 'request': request,
             'csrf_token': 'faketoken'},
            request=request)
        self.assertIn(
            'Analysis is being initialized on the Galaxy server', html)
        # Same shell as the ready state: styled panel, no Url field, no
        # info-refresh banner, same background-staged polling wired up.
        self.assertIn('history-meta-panel', html)
        self.assertIn('id="history-refreshable-region"', html)
        self.assertIn('id="history-refreshable-staging"', html)
        self.assertIn('window.historyRefreshTimer', html)
        self.assertNotIn('info-refresh', html)
        # The old Url field's own clipboard-copy button - a more precise
        # marker than the bare word "Url", which also legitimately shows
        # up in this file's own explanatory CSS comments.
        self.assertNotIn('data-clipboard-text', html)

    def test_exactly_one_dataset_still_shows_the_wait_message(self):
        request = self._request()
        html = render_to_string(
            'workspace/history.html',
            {'object': self._obj([
                {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
                 'visible': True, 'extension': 'fasta'},
            ]), 'request': request, 'csrf_token': 'faketoken'},
            request=request)
        self.assertIn(
            'Analysis is being initialized on the Galaxy server', html)
        self.assertNotIn('id="workflow-step-chain"', html)

    def test_wait_state_resolves_its_readiness_deferreds_immediately(self):
        """
        Nothing async happens in the wait state - the staging-swap script
        (history_contents_refreshable.html) would otherwise wait forever
        on deferreds that only the step-chain/table-building scripts
        (skipped here) ever resolve.
        """
        request = self._request()
        html = render_to_string(
            'workspace/include/history_contents_refreshable.html',
            {'object': self._obj([]), 'request': request,
             'csrf_token': 'faketoken', 'staging': True},
            request=request)
        self.assertIn('window.__historyStepChainReady.resolve();', html)
        self.assertIn('window.__historyTableReady.resolve();', html)
        self.assertIn(
            '$.when(window.__historyStepChainReady, '
            'window.__historyTableReady).done(function () {', html)


class ResolveDatasetToolsTest(TestCase):
    """
    Unit tests for resolve_dataset_tools/build_citations - the helpers
    HistoryContentRefreshView now uses to resolve dataset->tool_id/
    tool_id->name/citations once per poll, shared across the step chain,
    the table, and the citations list (previously three independent,
    redundant rounds of Galaxy-backed AJAX calls every 10s - see
    CLAUDE.md's step-chain section for the real production incident,
    pod restarts, this contributed to). Tested directly against a mocked
    `gi` rather than through the decorated view/HTTP layer - see
    HistoryContentRefreshViewTest below for that.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        with patch('tools.models.requests.get',
                   return_value=Mock(status_code=200, json=lambda: {
                       'id': 'mafft', 'name': 'ignored', 'version': '1.0',
                       'inputs': [], 'outputs': [],
                   })):
            self.tool = Tool.objects.create(
                galaxy_server=self.server, id_galaxy='mafft',
                name='MAFFT', description='Alignment', version='1.0')

    def test_resolves_tool_ids_and_prefers_the_local_tool_name(self):
        gi = Mock()
        gi.histories.show_dataset_provenance.side_effect = (
            lambda history_id, dataset_id, follow: {'tool_id': 'mafft'})

        dataset_tool_ids, tool_names = resolve_dataset_tools(
            gi, self.server, 'hist1', ['d1', 'd2'])

        self.assertEqual(dataset_tool_ids, {'d1': 'mafft', 'd2': 'mafft'})
        self.assertEqual(tool_names, {'mafft': 'MAFFT'})
        # The whole point: a tool already known locally shouldn't need a
        # Galaxy API call at all to get its name.
        gi.tools.show_tool.assert_not_called()

    def test_falls_back_to_galaxy_for_an_unknown_tool(self):
        gi = Mock()
        gi.histories.show_dataset_provenance.return_value = {'tool_id': 'unknown_tool'}
        gi.tools.show_tool.return_value = {'name': 'Unknown Tool'}

        dataset_tool_ids, tool_names = resolve_dataset_tools(
            gi, self.server, 'hist1', ['d1'])

        self.assertEqual(tool_names, {'unknown_tool': 'Unknown Tool'})
        gi.tools.show_tool.assert_called_once_with(tool_id='unknown_tool')

    def test_a_dataset_that_fails_to_resolve_is_just_left_out(self):
        gi = Mock()
        gi.histories.show_dataset_provenance.side_effect = requests.exceptions.Timeout('boom')

        dataset_tool_ids, tool_names = resolve_dataset_tools(
            gi, self.server, 'hist1', ['d1'])

        self.assertEqual(dataset_tool_ids, {})
        self.assertEqual(tool_names, {})

    def test_build_citations_includes_ngphylo_plus_resolved_tools(self):
        # self.tool has no Citation rows of its own here - just checks
        # ngphylo's own citation is always present regardless.
        refs = build_citations({'d1': 'mafft'})
        self.assertEqual(len(refs), 1)
        self.assertIn('NGPhylogeny.fr', refs[0])

    def test_build_citations_extends_with_a_resolved_tools_own_citations(self):
        from tools.models import Citation
        Citation.objects.create(
            tool=self.tool, reference='@article{a,title={MAFFT paper}}')
        refs = build_citations({'d1': 'mafft'})
        self.assertEqual(len(refs), 2)
        self.assertTrue(any('MAFFT paper' in r for r in refs))


class HistoryContentRefreshViewTest(TestCase):
    """
    Integration test for HistoryContentRefreshView going through the
    actual decorated view (galaxy.decorator.connection_galaxy), not just
    rendering the template directly like this file's other history-
    detail-page tests. tools.tests.GetToolNameViewTest already
    established that this is workable for a Galaxy-connected view: a
    real Server + anonymous GalaxyUser DB fixture, then patching the
    specific bioblend client methods actually called - no need to mock
    HTTP directly, and no need for CELERY_TASK_ALWAYS_EAGER either,
    since updateworkspacestatus.delay() is patched out (this test cares
    about what the view renders from the history_content already stored
    on the row, not about that Celery task refreshing it).
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)
        with patch('tools.models.requests.get',
                   return_value=Mock(status_code=200, json=lambda: {
                       'id': 'mafft', 'name': 'ignored', 'version': '1.0',
                       'inputs': [], 'outputs': [],
                   })):
            Tool.objects.create(
                galaxy_server=self.server, id_galaxy='mafft',
                name='MAFFT', description='Alignment', version='1.0')

        history_content = [
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'MAFFT alignment', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
        ]
        self.history = WorkspaceHistory.objects.create(
            history='hist1', name='Test run', email='', monitored=True,
            finished=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server,
            history_content_json=json.dumps(history_content),
            history_info_json=json.dumps({'id': 'hist1', 'name': 'Test run'}))

    def test_resolves_and_embeds_tool_names_with_no_client_ajax_left(self):
        with patch('bioblend.galaxy.histories.HistoryClient.show_dataset_provenance',
                   return_value={'tool_id': 'mafft'}), \
             patch('workspace.views.updateworkspacestatus.delay'):
            response = self.client.get(
                reverse('history_content_refresh', kwargs={'history_id': 'hist1'}))

        self.assertEqual(response.status_code, 200)
        html = response.content.decode()
        self.assertIn('MAFFT', html)
        # The whole point of this change: the rendered fragment no
        # longer contains any of the client-side AJAX calls it used to
        # make to resolve this same data every poll.
        self.assertNotIn('get_dataset_tool', html)
        self.assertNotIn('get_tool_name', html)
        self.assertNotIn('get_dataset_citations', html)
        self.assertIn('NGPhylogeny.fr', html)

    def test_history_detail_view_also_resolves_tool_names_on_first_load(self):
        """
        Regression test: dataset_tool_ids/tool_names/citations used to be
        computed only by HistoryContentRefreshView.get_context_data() -
        HistoryDetailView (the very first render, before any poll) never
        set them at all. Harmless for a still-running history (the first
        poll, moments later, would backfill it) but permanently broken
        for an *already-finished* one: object.finished being true clears
        the client's own polling timer before it ever fires once (see
        history_contents_provenance_ajax.html), so /refresh is never
        called and the step chain/table/citations stay unresolved
        forever. Real bug, hit live on a finished production run - see
        CLAUDE.md's step-chain section. Fixed by moving the computation
        up into WorkspaceHistoryObjectMixin.get_context_data(), shared by
        both views.
        """
        WorkspaceHistory.objects.filter(pk=self.history.pk).update(finished=True)
        with patch('bioblend.galaxy.histories.HistoryClient.show_dataset_provenance',
                   return_value={'tool_id': 'mafft'}), \
             patch('workspace.views.updateworkspacestatus.delay'):
            response = self.client.get(
                reverse('history_detail', kwargs={'history_id': 'hist1'}))
        self.assertEqual(response.status_code, 200)
        html = response.content.decode()
        self.assertIn('MAFFT', html)
        self.assertIn('"d1": "mafft"', html)

    def test_malformed_history_content_does_not_crash_the_view(self):
        """
        Regression test: hit live on a real deploy-dev history - one
        specific history's history_content_json held plain strings, not
        the usual list of dataset dicts (AttributeError: 'str' object
        has no attribute 'get', 500ing the whole page from
        get_context_data's own `f.get('id')`). The template rendering
        this exact same data never crashed on it (Django template
        attribute lookups fail silently, not by raising) - this should
        degrade the same way, not crash.
        """
        WorkspaceHistory.objects.filter(pk=self.history.pk).update(
            history_content_json=json.dumps(['not', 'a', 'list', 'of', 'dicts']))
        with patch('workspace.views.updateworkspacestatus.delay'):
            response = self.client.get(
                reverse('history_detail', kwargs={'history_id': 'hist1'}))
        self.assertEqual(response.status_code, 200)

    def test_wait_state_does_not_call_galaxy_for_tool_resolution(self):
        """
        Fewer than 2 datasets - resolve_dataset_tools shouldn't be
        called at all (nothing to resolve yet, see the template's own
        "please wait" branch).
        """
        WorkspaceHistory.objects.filter(pk=self.history.pk).update(
            history_content_json=json.dumps([
                {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
                 'visible': True, 'extension': 'fasta'},
            ]))
        with patch('bioblend.galaxy.histories.HistoryClient.show_dataset_provenance') as prov, \
             patch('workspace.views.updateworkspacestatus.delay'):
            response = self.client.get(
                reverse('history_content_refresh', kwargs={'history_id': 'hist1'}))
        self.assertEqual(response.status_code, 200)
        prov.assert_not_called()
        self.assertIn('Analysis is being initialized', response.content.decode())


class TreePreviewTest(TestCase):
    """
    Regression tests for the inline phylotree.js tree preview
    (history_contents_refreshable.html, built by workspace.views.
    WorkspaceHistoryObjectMixin.get_context_data's tree_preview_dataset_id/
    tree_preview_url) - shown above the dataset table, below the
    workflow step chain, as soon as a finished nhx/nwk dataset exists in
    the history. Same setUp fixture as HistoryContentRefreshViewTest.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)
        # Resolves locally (Tool.objects...) rather than falling back to
        # a real Galaxy API call - same reasoning as
        # HistoryContentRefreshViewTest's own fixture: this test only
        # cares about the tree-preview wiring, not tool-name resolution,
        # so avoid a real outbound HTTP attempt to fake-galaxy.example.org.
        with patch('tools.models.requests.get',
                   return_value=Mock(status_code=200, json=lambda: {
                       'id': 'mafft', 'name': 'ignored', 'version': '1.0',
                       'inputs': [], 'outputs': [],
                   })):
            Tool.objects.create(
                galaxy_server=self.server, id_galaxy='mafft',
                name='MAFFT', description='Alignment', version='1.0')

    def _make_history(self, history_content):
        return WorkspaceHistory.objects.create(
            history='hist1', name='Test run', email='', monitored=True,
            finished=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server,
            history_content_json=json.dumps(history_content),
            history_info_json=json.dumps({'id': 'hist1', 'name': 'Test run'}))

    def test_tree_preview_wired_up_once_a_finished_tree_dataset_exists(self):
        self._make_history([
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'tree.nhx', 'state': 'ok',
             'visible': True, 'extension': 'nhx'},
        ])
        with patch('bioblend.galaxy.histories.HistoryClient.show_dataset_provenance',
                   return_value={'tool_id': 'mafft'}), \
             patch('workspace.views.updateworkspacestatus.delay'):
            response = self.client.get(
                reverse('history_content_refresh', kwargs={'history_id': 'hist1'}))
        self.assertEqual(response.status_code, 200)
        html = response.content.decode()
        self.assertIn('var treePreviewDatasetId = "d2";', html)
        self.assertIn(reverse('display_raw', kwargs={'file_id': 'd2'}), html)
        self.assertIn('id="tree-preview-anchor"', html)
        # Positioned below the step chain, above the dataset table.
        self.assertLess(html.index('id="workflow-step-chain"'),
                         html.index('id="tree-preview-anchor"'))
        self.assertLess(html.index('id="tree-preview-anchor"'),
                         html.index('id="myTable"'))

    def test_no_tree_preview_wiring_without_a_finished_tree_dataset(self):
        self._make_history([
            {'id': 'd1', 'hid': 1, 'name': 'input.fasta', 'state': 'ok',
             'visible': True, 'extension': 'fasta'},
            {'id': 'd2', 'hid': 2, 'name': 'tree.nhx', 'state': 'running',
             'visible': True, 'extension': 'nhx'},
        ])
        with patch('bioblend.galaxy.histories.HistoryClient.show_dataset_provenance',
                   return_value={'tool_id': 'mafft'}), \
             patch('workspace.views.updateworkspacestatus.delay'):
            response = self.client.get(
                reverse('history_content_refresh', kwargs={'history_id': 'hist1'}))
        self.assertEqual(response.status_code, 200)
        html = response.content.decode()
        self.assertIn('var treePreviewDatasetId = "";', html)
        self.assertIn('id="tree-preview-anchor"', html)


class WorkspacePermalinkTest(TestCase):
    """
    Tests for the workspace/histories permalink feature -
    PreviousHistoryListView.get_context_data()'s permalink_url and
    WorkspacePermalinkView (/workspace/permalink/<token>). The whole
    point: "Workspace" (previous_analyses) only ever shows what's in the
    *current* browser session's 'histories' list - this is the
    session-independent way back to that same list (a different
    browser, a cleared session, or a link handed to a collaborator).
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)
        self.history = WorkspaceHistory.objects.create(
            history='hist1', name='Test run', email='', monitored=True,
            finished=True, deleted=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server)

    def _set_session_histories(self, history_ids):
        session = self.client.session
        session['histories'] = history_ids
        session.save()

    def test_previous_analyses_has_no_permalink_with_an_empty_session(self):
        response = self.client.get(reverse('previous_analyses'))
        self.assertEqual(response.status_code, 200)
        self.assertIsNone(response.context['permalink_url'])

    def test_previous_analyses_shows_a_permalink_url(self):
        self._set_session_histories(['hist1'])
        response = self.client.get(reverse('previous_analyses'))
        permalink_url = response.context['permalink_url']
        self.assertIsNotNone(permalink_url)
        self.assertIn('/workspace/permalink/', permalink_url)
        self.assertContains(response, 'permalink-url')

    def test_following_the_permalink_restores_the_list_in_a_fresh_session(self):
        token = signing.dumps(['hist1'], salt=PERMALINK_SALT)

        # A brand new client - no cookies/session shared with whatever
        # session originally ran 'hist1' - is exactly the scenario this
        # feature exists for.
        from django.test import Client
        fresh_client = Client()
        response = fresh_client.get(
            reverse('workspace_permalink', kwargs={'token': token}),
            follow=True)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(fresh_client.session['histories'], ['hist1'])
        self.assertContains(response, 'Test run')

    def test_following_the_permalink_merges_with_existing_session_histories(self):
        WorkspaceHistory.objects.create(
            history='hist2', name='Second run', email='', monitored=True,
            finished=True, deleted=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server)
        self._set_session_histories(['hist2'])
        token = signing.dumps(['hist1'], salt=PERMALINK_SALT)

        self.client.get(
            reverse('workspace_permalink', kwargs={'token': token}))

        self.assertEqual(
            set(self.client.session['histories']), {'hist1', 'hist2'})

    def test_tampered_token_redirects_with_an_error_instead_of_crashing(self):
        response = self.client.get(
            reverse('workspace_permalink', kwargs={'token': 'not-a-real-token'}),
            follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'invalid')
