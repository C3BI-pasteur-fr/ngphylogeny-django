# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import copy
from unittest.mock import MagicMock, patch

from django.conf import settings
from django.test import TestCase, override_settings

from .models import BlastRun
from .tasks import launch_pasteur_blast, checkblastruns


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


_PASTEUR_BLASTN = (
    'toolshed.pasteur.fr/repos/fmareuil/ncbi_blast_plus/'
    'ncbi_blastn_wrapper/2.6.0'
)


def _blasts_with_pasteur_activated():
    blasts = copy.deepcopy(settings.BLASTS)
    blasts['pasteur']['activated'] = True
    return blasts


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
            b.id, ">s1\nACGTACGTAC\n", _PASTEUR_BLASTN, 'nt', 0.00001, 0.8, 10)

        b.refresh_from_db()
        self.assertEqual(b.status, BlastRun.PENDING)
        self.assertFalse(b.message)
        galaxycon.tools.upload_file.assert_called_once()
        galaxycon.tools.run_tool.assert_called_once()

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
