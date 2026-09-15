# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import copy
from unittest.mock import MagicMock, patch

from django.conf import settings
from django.test import TestCase, override_settings

from .models import BlastRun
from .tasks import launch_pasteur_blast


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
