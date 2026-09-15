# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.test import TestCase


class BlastViewDisabledTest(TestCase):
    """
    BLAST analysis is temporarily disabled (both submission paths had
    open issues - see blast/views.py's BlastView.dispatch()). GET and
    POST should both get the disabled notice, not the real form or
    launch a run, and this must hold before any form processing (a POST
    with no data at all should not somehow reach form validation first).
    """

    def test_get_shows_disabled_notice(self):
        response = self.client.get('/blast/')
        self.assertEqual(response.status_code, 503)
        self.assertTemplateUsed(response, 'blast/blast_disabled.html')
        self.assertTemplateNotUsed(response, 'blast/blast.html')

    def test_post_does_not_launch_a_run(self):
        response = self.client.post('/blast/', {})
        self.assertEqual(response.status_code, 503)
        self.assertTemplateUsed(response, 'blast/blast_disabled.html')
