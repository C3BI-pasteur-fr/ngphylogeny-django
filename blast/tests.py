# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import copy
from datetime import timedelta
from unittest.mock import MagicMock, patch

from celery.exceptions import SoftTimeLimitExceeded
from django.conf import settings
from django.core import mail
from django.test import TestCase, override_settings
from django.utils import timezone

from .emails import build_blast_completion_email, send_blast_completion_email
from .models import BlastRun, BlastSubject
from .tasks import (PASTEUR_RUN_STALE_AFTER, deleteoldblastruns,
                     launch_ncbi_blast, launch_pasteur_blast, checkblastruns)


class BlastViewTest(TestCase):
    """
    Reactivated after a temporary disable (both submission paths had open
    issues - launch_ncbi_blast's missing timeout, now bounded to 10
    minutes, and the Pasteur BLAST server option, now env-var driven -
    see CLAUDE.md's "Kubernetes deployment" section for both). Just a
    smoke test that the real form is back, not the disabled notice - see
    NGPhylogeny_fr_django's git history for BlastViewDisabledTest, the
    regression test the disable itself had.
    """

    def test_get_shows_the_real_form(self):
        response = self.client.get('/blast/')
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'blast/blast.html')
        self.assertTemplateNotUsed(response, 'blast/blast_disabled.html')


class DeleteBlastRunViewTest(TestCase):
    """
    Regression test: clicking "delete" on a blast run used to crash with
    "TemplateDoesNotExist: blast/blastrun_confirm_delete.html". Django
    4.x's BaseDeleteView.post() was rewritten to go through FormMixin
    (get_form()/form_valid()/form_invalid()) and no longer calls
    self.delete() at all - the view's get() (forwarding to self.post()
    to skip DeleteView's confirmation page, since this view never had
    one) built an unbound form (get_form_kwargs() only binds
    request.POST/FILES when request.method is actually 'POST'), so
    form.is_valid() was always False and it fell through to rendering
    the (never-created, never-wanted) confirmation template instead of
    soft-deleting and redirecting.
    """

    def _run(self):
        return BlastRun.objects.create(query_id="", query_seq="")

    def test_get_soft_deletes_and_redirects(self):
        b = self._run()
        response = self.client.get('/blast/%s/delete' % b.id)
        self.assertRedirects(response, '/blast/')
        b.refresh_from_db()
        self.assertTrue(b.deleted)

    def test_post_soft_deletes_and_redirects(self):
        b = self._run()
        response = self.client.post('/blast/%s/delete' % b.id)
        self.assertRedirects(response, '/blast/')
        b.refresh_from_db()
        self.assertTrue(b.deleted)


def _pasteur_blastn_prog():
    """
    The settings.BLASTS['pasteur']['progs'] key for the blastn wrapper -
    derived by 'type', not hardcoded as a literal tool id/version string,
    since that id has already changed once (toolshed.pasteur.fr/repos/
    fmareuil/... -> toolshed.g2.bx.psu.edu/repos/devteam/...) independently
    of anything in this test file, silently breaking it until the id was
    updated here too.
    """
    for prog, cfg in settings.BLASTS['pasteur']['progs'].items():
        if cfg['type'] == 'blastn':
            return prog
    raise AssertionError(
        "No 'blastn' entry in settings.BLASTS['pasteur']['progs']")


def _blasts_with_pasteur_activated():
    blasts = copy.deepcopy(settings.BLASTS)
    blasts['pasteur']['activated'] = True
    return blasts


class BlastAjaxEndpointsTest(TestCase):
    """
    Regression test: available_blasts_dbs/blast_example's URL routes
    (blast/urls.py) take the Galaxy toolshed tool id as a path segment,
    inserted as-is by templates/blast/blast.html's JS (url.replace(...),
    not encodeURIComponent-escaped). Their prog regex ([\\w/\\.]+) didn't
    allow "+" - which Pasteur's current wrapper ids all contain (e.g.
    "...ncbi_blastn_wrapper/2.14.1+galaxy2") - so both routes 404ed for
    every Pasteur program: the database dropdown came back empty and the
    example-sequence button did nothing, with nothing but a 404 in the
    browser console to show why.
    """

    @override_settings(BLASTS=_blasts_with_pasteur_activated())
    def test_dbs_endpoint_resolves_for_a_prog_id_containing_a_plus(self):
        prog = _pasteur_blastn_prog()
        self.assertIn('+', prog)  # otherwise this isn't testing the bug

        response = self.client.get('/blast/dbs/pasteur/%s' % prog)

        self.assertEqual(response.status_code, 200)
        self.assertIn('nt', response.json())

    @override_settings(BLASTS=_blasts_with_pasteur_activated())
    def test_example_endpoint_resolves_for_a_prog_id_containing_a_plus(self):
        prog = _pasteur_blastn_prog()

        response = self.client.get('/blast/example/pasteur/%s' % prog)

        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.json())  # non-empty: the fasta content


class LaunchNcbiBlastTest(TestCase):
    """
    query_length used to only be derivable from query_seq - which
    deleteoldblastruns() later clears to free space, losing the length
    entirely. launch_ncbi_blast() now sets it right alongside query_seq,
    before anything else (including the alphabet check below) can
    short-circuit the run into ERROR - deliberately triggering that
    branch here (a protein sequence submitted to blastn) so this stays a
    fast, no-network test rather than needing to mock NCBIWWW.qblast().
    """

    def test_query_length_set_even_when_alphabet_check_fails(self):
        sequence = "MASGILVNVKEEVTCPICLE"
        b = BlastRun.objects.create(query_id="", query_seq="")

        launch_ncbi_blast(
            b.id, ">s1\n%s\n" % sequence, 'blastn', 'nr', 0.00001, 0.8, 10)

        b.refresh_from_db()
        self.assertEqual(b.status, BlastRun.ERROR)
        self.assertIn('wrong alphabet', b.message.lower())
        self.assertEqual(b.query_length, len(sequence))


class LaunchPasteurBlastTest(TestCase):
    """
    Regression test: launch_pasteur_blast() used to crash with "TypeError:
    a bytes-like object is required, not 'str'" - it wrote the (plain str)
    query sequence straight into a NamedTemporaryFile(), which defaults to
    binary mode. Same Python 2->3 bug class as
    workflows.tests.ProcessFileToUploadTest - only caught once Pasteur
    BLAST was actually activated and submitted for real (see CLAUDE.md's
    Pasteur BLAST activation section), not by any existing test.
    """

    @override_settings(BLASTS=_blasts_with_pasteur_activated())
    @patch('blast.tasks.galaxy_connection')
    def test_pasted_text_str_input_does_not_crash(self, mock_galaxy_connection):
        galaxycon = MagicMock()
        galaxycon.histories.create_history.return_value = {'id': 'fakehistoryid'}
        galaxycon.tools.upload_file.return_value = {
            'outputs': [{'id': 'fakefileid'}]
        }
        galaxycon.tools.run_tool.return_value = {
            'outputs': [{'id': 'fakeoutputid'}]
        }
        mock_galaxy_connection.return_value = galaxycon

        b = BlastRun.objects.create(query_id="", query_seq="")

        launch_pasteur_blast(
            b.id, ">s1\nACGTACGTAC\n", _pasteur_blastn_prog(), 'nt', 0.00001, 0.8, 10)

        b.refresh_from_db()
        self.assertEqual(b.status, BlastRun.PENDING)
        self.assertFalse(b.message)
        self.assertEqual(b.query_length, len("ACGTACGTAC"))
        galaxycon.tools.upload_file.assert_called_once()
        galaxycon.tools.run_tool.assert_called_once()

    @override_settings(BLASTS=_blasts_with_pasteur_activated())
    @patch('blast.tasks.deletegalaxyhistory')
    @patch('blast.tasks.galaxy_connection')
    def test_timeout_marks_error_and_cleans_up_galaxy_history(
            self, mock_galaxy_connection, mock_deletegalaxyhistory):
        """
        Regression test: launch_pasteur_blast() used to have no timeout
        at all (unlike launch_ncbi_blast, which got one after NCBI's
        qblast() turned out to hang indefinitely) - a hang in the
        (network-bound) create_history/upload_file/run_tool calls that
        submit the job would leave the run stuck forever, with the
        Galaxy history it already created orphaned on the Galaxy server
        since nothing else would ever clean it up before the 14-day
        deleteoldblastruns() cutoff.
        """
        galaxycon = MagicMock()
        galaxycon.histories.create_history.return_value = {'id': 'fakehistoryid'}
        galaxycon.tools.upload_file.side_effect = SoftTimeLimitExceeded()
        mock_galaxy_connection.return_value = galaxycon

        b = BlastRun.objects.create(query_id="", query_seq="")

        launch_pasteur_blast(
            b.id, ">s1\nACGTACGTAC\n", _pasteur_blastn_prog(), 'nt', 0.00001, 0.8, 10)

        b.refresh_from_db()
        self.assertEqual(b.status, BlastRun.ERROR)
        self.assertIn('too long', b.message)
        mock_deletegalaxyhistory.delay.assert_called_once_with('fakehistoryid')

    def test_history_fields_are_wide_enough_for_real_galaxy_ids(self):
        """
        Regression test: history/history_fileid used to be
        CharField(max_length=20) - too narrow for this Galaxy server's
        actual encoded dataset ids. launch_pasteur_blast() would run the
        blast job for real, then crash saving the result:
        "django.db.utils.DataError: value too long for type character
        varying(20)", leaving the run stuck showing PENDING in
        NGPhylogeny while it kept running/finished on Galaxy. Checked at
        the field level, not via an actual oversized save() - CI's test
        DB is sqlite (no NGPHYLO_DATABASE_HOST set - see settings/base.py
        and .gitlab-ci.yml's test job), which doesn't enforce CharField
        max_length at the DB layer the way the real deployment's Postgres
        does, so a save()-based test can't reproduce this crash here.
        """
        self.assertGreaterEqual(
            BlastRun._meta.get_field('history').max_length, 250)
        self.assertGreaterEqual(
            BlastRun._meta.get_field('history_fileid').max_length, 250)


class CheckBlastRunsTest(TestCase):
    """
    Regression test: checkblastruns() (the every-minute Celery-beat task
    that polls pending/running Pasteur runs) used to crash with
    "AttributeError: 'list' object has no attribute 'get'" - a real race
    with launch_pasteur_blast(), which saves a run as PENDING right after
    creating its Galaxy history but only sets history_fileid afterwards,
    once the (network-bound) file upload + tool run calls complete. If
    this task's schedule fires in that window, show_dataset(b.history,
    '') hits Galaxy's history *contents list* endpoint (trailing empty
    dataset id) instead of a single dataset, and gets a list back. The
    whole per-run loop also used to share one try/except, so this alone
    silently aborted checking every other pending/running run in the same
    pass too - both are covered here.
    """

    def _pasteur_run(self, history_fileid, status=BlastRun.PENDING):
        return BlastRun.objects.create(
            query_id="", query_seq="", server=BlastRun.PASTEUR,
            status=status, history='fakehistory',
            history_fileid=history_fileid)

    @patch('blast.tasks.galaxy_connection')
    def test_skips_runs_without_history_fileid_yet(self, mock_galaxy_connection):
        pending_submission = self._pasteur_run(history_fileid='')
        ready_run = self._pasteur_run(history_fileid='realfileid')

        galaxycon = MagicMock()
        galaxycon.histories.show_dataset.return_value = {
            'state': 'running', 'misc_info': ''}
        mock_galaxy_connection.return_value = galaxycon

        checkblastruns()

        galaxycon.histories.show_dataset.assert_called_once_with(
            'fakehistory', 'realfileid')
        pending_submission.refresh_from_db()
        self.assertEqual(pending_submission.status, BlastRun.PENDING)
        ready_run.refresh_from_db()
        self.assertEqual(ready_run.status, BlastRun.RUNNING)

    @patch('blast.tasks.galaxy_connection')
    def test_one_runs_failure_does_not_block_the_others(self, mock_galaxy_connection):
        broken_run = self._pasteur_run(history_fileid='brokenfileid')
        healthy_run = self._pasteur_run(history_fileid='healthyfileid')

        def fake_show_dataset(history, fileid):
            if fileid == 'brokenfileid':
                raise AttributeError("'list' object has no attribute 'get'")
            return {'state': 'running', 'misc_info': ''}

        galaxycon = MagicMock()
        galaxycon.histories.show_dataset.side_effect = fake_show_dataset
        mock_galaxy_connection.return_value = galaxycon

        checkblastruns()

        broken_run.refresh_from_db()
        self.assertEqual(broken_run.status, BlastRun.ERROR)
        healthy_run.refresh_from_db()
        self.assertEqual(healthy_run.status, BlastRun.RUNNING)

    @patch('blast.tasks.deletegalaxyhistory')
    @patch('blast.tasks.galaxy_connection')
    def test_gives_up_on_runs_stuck_past_the_staleness_cutoff(
            self, mock_galaxy_connection, mock_deletegalaxyhistory):
        """
        Regression test: unlike launch_ncbi_blast/launch_pasteur_blast's
        own soft_time_limit/time_limit (which only cover *submitting* the
        job), there used to be no timeout at all on the actual Galaxy-side
        blast computation that checkblastruns() polls - a run genuinely
        stuck "running" on Galaxy's/the cluster's side (not just a slow
        legitimate search) would poll forever, showing as permanently
        "Running" with no way to ever notice or recover.
        """
        stale_run = self._pasteur_run(
            history_fileid='stalefileid', status=BlastRun.RUNNING)
        BlastRun.objects.filter(pk=stale_run.pk).update(
            date=timezone.now() - PASTEUR_RUN_STALE_AFTER - timedelta(minutes=1))
        fresh_run = self._pasteur_run(
            history_fileid='freshfileid', status=BlastRun.RUNNING)

        galaxycon = MagicMock()
        galaxycon.histories.show_dataset.return_value = {
            'state': 'running', 'misc_info': ''}
        mock_galaxy_connection.return_value = galaxycon

        checkblastruns()

        stale_run.refresh_from_db()
        self.assertEqual(stale_run.status, BlastRun.ERROR)
        self.assertIn('longer than expected', stale_run.message)
        mock_deletegalaxyhistory.delay.assert_called_once_with('fakehistory')
        galaxycon.histories.show_dataset.assert_called_once_with(
            'fakehistory', 'freshfileid')
        fresh_run.refresh_from_db()
        self.assertEqual(fresh_run.status, BlastRun.RUNNING)


class DeleteOldBlastRunsTest(TestCase):
    """
    blast.tasks.deleteoldblastruns() (the daily 2am cleanup, 14-day
    cutoff on BlastRun.date). Covers: the deletegalaxyhistory call now
    being queued rather than blocking the rest of the batch, and
    query_seq/tree now being cleared to free space on old runs - safe
    for the daily report (workspace/reports.py), which only ever reads
    BlastRun's date/deleted/id, never query_seq/tree.
    """

    def _run(self, days_ago, server=BlastRun.NCBI, history=''):
        b = BlastRun.objects.create(
            query_id="", query_seq="ACGT", tree="(a,b);",
            server=server, history=history)
        BlastRun.objects.filter(pk=b.pk).update(
            date=timezone.now() - timedelta(days=days_ago))
        b.refresh_from_db()
        return b

    @patch('blast.tasks.deletegalaxyhistory')
    def test_clears_query_seq_and_tree_on_old_runs(self, mock_deletegalaxyhistory):
        old_run = self._run(days_ago=15)
        BlastSubject.objects.create(
            subject_id='s1', subject_seq='ACGT', subject_fullseq='ACGT',
            blastrun=old_run)

        deleteoldblastruns()

        old_run.refresh_from_db()
        self.assertTrue(old_run.deleted)
        self.assertEqual(old_run.query_seq, "")
        self.assertEqual(old_run.tree, "")
        self.assertEqual(BlastSubject.objects.filter(blastrun=old_run).count(), 0)

    @patch('blast.tasks.deletegalaxyhistory')
    def test_derives_query_length_before_clearing_query_seq(
            self, mock_deletegalaxyhistory):
        """
        Regression guard: query_length is normally set at submission
        time (blast/tasks.py's launch_ncbi_blast/launch_pasteur_blast),
        but a row predating that field (or from some other path that
        never set it) would otherwise lose the length entirely once
        query_seq is cleared here - re-derived from query_seq, before
        it's cleared, only when not already set.
        """
        never_set = self._run(days_ago=15)  # _run() leaves query_length unset
        already_set = self._run(days_ago=15)
        BlastRun.objects.filter(pk=already_set.pk).update(query_length=999)

        deleteoldblastruns()

        never_set.refresh_from_db()
        self.assertEqual(never_set.query_length, len("ACGT"))
        already_set.refresh_from_db()
        self.assertEqual(already_set.query_length, 999)

    @patch('blast.tasks.deletegalaxyhistory')
    def test_leaves_recent_runs_untouched(self, mock_deletegalaxyhistory):
        recent_run = self._run(days_ago=1)

        deleteoldblastruns()

        recent_run.refresh_from_db()
        self.assertFalse(recent_run.deleted)
        self.assertEqual(recent_run.query_seq, "ACGT")
        self.assertEqual(recent_run.tree, "(a,b);")

    @patch('blast.tasks.deletegalaxyhistory')
    def test_queues_galaxy_history_deletion_instead_of_blocking(
            self, mock_deletegalaxyhistory):
        pasteur_run = self._run(
            days_ago=15, server=BlastRun.PASTEUR, history='realhistoryid')
        ncbi_run = self._run(days_ago=15, server=BlastRun.NCBI, history='')

        deleteoldblastruns()

        mock_deletegalaxyhistory.assert_not_called()
        mock_deletegalaxyhistory.delay.assert_called_once_with('realhistoryid')
        pasteur_run.refresh_from_db()
        ncbi_run.refresh_from_db()
        self.assertTrue(pasteur_run.deleted)
        self.assertTrue(ncbi_run.deleted)


class BlastCompletionEmailTest(TestCase):
    """
    blast.emails.build_blast_completion_email() used to be a hand-built
    plain-text send_mail() call in blast/tasks.py - this now reuses the
    same branded HTML template/MIME wiring as the workflow job-completion
    email (workspace.emails), so BLAST notifications look the same as
    every other NGPhylogeny.fr email. Mirrors
    workspace.tests.JobCompletionEmailTest's coverage of that shared
    machinery.
    """

    def _run(self, status=BlastRun.FINISHED):
        return BlastRun.objects.create(
            query_id="", query_seq="", status=status)

    @override_settings(NGPHYLO_REPORT_FROM_EMAIL='ngphylogeny@pasteur.fr')
    def test_success_email_content_and_structure(self):
        b = self._run(status=BlastRun.FINISHED)
        msg = build_blast_completion_email(b, 'user@example.org')

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
        self.assertIn(str(b.id), html_body)
        self.assertIn('doi.org/10.1093/nar/gkz303', html_body)

        # Same multipart/related + inline logo structure as the workflow
        # job-completion email - see JobCompletionEmailTest's own version
        # of this check for why it matters.
        self.assertEqual(msg.mixed_subtype, 'related')
        self.assertNotIn('data:image', html_body)
        self.assertEqual(len(msg.attachments), 2)

    def test_error_email_shows_error_status(self):
        b = self._run(status=BlastRun.ERROR)
        msg = build_blast_completion_email(b, 'user@example.org')
        self.assertIn('error', msg.subject.lower())
        html_body, _ = msg.alternatives[0]
        self.assertIn('Finished with errors', html_body)
        self.assertNotIn('Finished successfully', html_body)

    @override_settings(NGPHYLO_HTTPS_HOST='ngphylogeny.fr')
    def test_results_link_uses_https_when_configured(self):
        b = self._run()
        msg = build_blast_completion_email(b, 'user@example.org')
        html_body, _ = msg.alternatives[0]
        self.assertIn('https://ngphylogeny.fr/blast/%s' % b.id, html_body)

    def test_send_blast_completion_email_actually_sends(self):
        b = self._run()
        send_blast_completion_email(b, 'someone@example.org')
        self.assertEqual(len(mail.outbox), 1)
        self.assertEqual(mail.outbox[0].to, ['someone@example.org'])
