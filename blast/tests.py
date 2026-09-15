# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TestCase


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
