from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.test import TestCase

from .models import GalaxyUser, Server


class GalaxyUserGetGalaxyInstanceTest(TestCase):
    """
    Regression test: GalaxyUser.get_galaxy_instance used to construct its
    bioblend GalaxyInstance with bioblend's own default (timeout=None,
    wait forever) - a request handler blocked on one slow/unresponsive
    Galaxy call ties up one of a fixed, small number of uwsgi worker
    slots (see manifest.yaml's "--processes 4 --threads 2") until Galaxy
    answers or uwsgi's own 120s harakiri kills the worker outright. A
    real production incident (a galaxy.pasteur.fr 502 spell, now hit far
    more often since the history detail page's step chain/table polls
    Galaxy-backed AJAX endpoints every 10s - see CLAUDE.md) piled up
    enough stuck requests that the liveness probe itself (a plain GET
    /status) couldn't get served in time, and Kubernetes restarted the
    pod repeatedly (exit code 137/SIGKILL). Fixed by bounding every
    Galaxy call to GALAXY_REQUEST_TIMEOUT seconds.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        self.user = User.objects.create_user(username='testuser')

    def test_timeout_is_bounded_not_left_at_bioblends_own_default(self):
        gu = GalaxyUser.objects.create(
            user=self.user, galaxy_server=self.server, api_key='fakekey')
        gi = gu.get_galaxy_instance
        self.assertEqual(gi.timeout, GalaxyUser.GALAXY_REQUEST_TIMEOUT)
        # Bioblend's own default is None (wait forever) - pin the actual
        # bound to a real number too, not just "not None", so a future
        # change that raises it back toward "basically unbounded" (e.g.
        # accidentally passing 0 or a huge value) still fails loudly.
        self.assertEqual(GalaxyUser.GALAXY_REQUEST_TIMEOUT, 30)

    def test_raises_without_an_api_key(self):
        gu = GalaxyUser.objects.create(
            user=self.user, galaxy_server=self.server, api_key='')
        with self.assertRaises(ValueError):
            gu.get_galaxy_instance
