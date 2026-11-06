from unittest.mock import Mock, patch

import requests
from django.contrib.auth.models import AnonymousUser, User
from django.test import RequestFactory, TestCase
from django.urls import reverse

from account.models import UserProfile
from workspace.models import WorkspaceHistory
from .decorator import connection_galaxy
from .galaxylib import GalaxyInstance
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


class ConnectionGalaxySharedKeyTest(TestCase):
    """
    galaxy.decorator.connection_galaxy: every visitor - authenticated
    into an NGPhylogeny account or not - now authenticates to Galaxy
    through the same single, shared "anonymous" GalaxyUser (see that
    module's own docstring - NGPhylogeny holds one Galaxy identity,
    independent of any individual NGPhylogeny account). Used to
    get_or_create() a *personal* GalaxyUser for an authenticated
    request.user and redirect to a URL name ('galaxy_account') that
    never actually existed in this project's URLconf if it had no
    api_key - an uncaught NoReverseMatch 500 on every authenticated
    request with no personal key, never caught because nothing
    exercised that path until now.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        shared_owner = User.objects.create_user(username='ngphylo-shared')
        self.shared = GalaxyUser.objects.create(
            user=shared_owner, galaxy_server=self.server,
            api_key='shared-key', anonymous=True)
        self.user = User.objects.create_user(
            username='alice', password='secretpass')

    @staticmethod
    @connection_galaxy
    def _probe_view(request):
        from django.http import HttpResponse
        return HttpResponse(request.galaxy.key)

    def test_authenticated_user_with_no_personal_galaxyuser_still_works(self):
        request = RequestFactory().get('/')
        request.user = self.user
        request.session = {}
        response = self._probe_view(request)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content.decode(), 'shared-key')
        # The old behavior used to create one on the fly - confirm this
        # is genuinely gone, not just made to not crash.
        self.assertFalse(
            GalaxyUser.objects.filter(user=self.user).exists())

    def test_anonymous_user_gets_the_same_shared_key(self):
        request = RequestFactory().get('/')
        request.user = AnonymousUser()
        request.session = {}
        response = self._probe_view(request)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content.decode(), 'shared-key')


class AccountPageOwnHistoriesTest(TestCase):
    """
    /account (account.views.AccountDetailView) lists the logged-in
    account's own WorkspaceHistory rows - a real, stable list tied to
    the account itself (workspace.views.create_history() already sets
    WorkspaceHistory.user to request.user for every authenticated
    submission, nothing needed there), unlike the session-based
    "Workspace" list (workspace.views.PreviousHistoryListView), which
    only ever exists in the one browser session that ran each analysis.

    Moved here from account/tests.py once the page itself moved out of
    galaxy.views.UpdateApiKey (a *personal* Galaxy API key form,
    removed along with that whole concept - see
    ConnectionGalaxySharedKeyTest above) into account.views.
    AccountDetailView, a plain account/history page with no Galaxy
    form left on it at all.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        self.user = User.objects.create_user(
            username='alice', password='secretpass')
        self.other_user = User.objects.create_user(
            username='bob', password='secretpass')

    def _make_history(self, user, history_id, deleted=False):
        # WorkspaceHistory.save() always calls rename(), which looks up
        # a GalaxyUser (get_galaxy_user(), workspace/models.py) - now
        # falls back to the shared anonymous GalaxyUser rather than
        # raising when there's no personal one, but there's still none
        # of *any* kind in this fixture - patched out since rename()'s
        # own Galaxy-side behavior isn't what these tests are about.
        with patch('workspace.models.WorkspaceHistory.rename'):
            return WorkspaceHistory.objects.create(
                history=history_id, name='Analyse ' + history_id, email='',
                monitored=True, finished=True, source_ip='127.0.0.1',
                workflow_category='OneClick', workflow_steps='',
                galaxy_server=self.server, user=user, deleted=deleted)

    def test_anonymous_user_is_redirected_to_login(self):
        response = self.client.get(reverse('account'))
        self.assertEqual(response.status_code, 302)
        self.assertIn('/account/login', response.url)

    def test_shows_only_the_logged_in_users_own_histories(self):
        self._make_history(self.user, 'mine1')
        self._make_history(self.other_user, 'not-mine')

        self.client.login(username='alice', password='secretpass')
        response = self.client.get(reverse('account'))

        self.assertEqual(response.status_code, 200)
        histories = list(response.context['histories'])
        self.assertEqual([h.history for h in histories], ['mine1'])

    def test_deleted_histories_are_excluded(self):
        self._make_history(self.user, 'live1')
        self._make_history(self.user, 'gone1', deleted=True)

        self.client.login(username='alice', password='secretpass')
        response = self.client.get(reverse('account'))

        histories = list(response.context['histories'])
        self.assertEqual([h.history for h in histories], ['live1'])

    def test_no_histories_renders_the_empty_state_not_a_crash(self):
        self.client.login(username='alice', password='secretpass')
        response = self.client.get(reverse('account'))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(list(response.context['histories']), [])
        self.assertContains(response, 'no analyses on this account')

    def test_get_or_creates_a_user_profile(self):
        self.assertFalse(UserProfile.objects.filter(user=self.user).exists())
        self.client.login(username='alice', password='secretpass')
        response = self.client.get(reverse('account'))
        self.assertIsInstance(response.context['profile'], UserProfile)
        self.assertTrue(UserProfile.objects.filter(user=self.user).exists())

    def test_saving_an_authenticated_users_history_falls_back_to_the_shared_key(self):
        # Regression test for the actual crash this fallback exists to
        # prevent: WorkspaceHistory.save() -> rename() ->
        # get_galaxy_user() used to plain .get() a personal GalaxyUser
        # and raise DoesNotExist for any authenticated user with none -
        # i.e. every authenticated user now that no personal GalaxyUser
        # is ever created (see ConnectionGalaxySharedKeyTest). Deliberately
        # not patching rename() here, unlike _make_history() above.
        shared_owner = User.objects.create_user(username='ngphylo-shared')
        GalaxyUser.objects.create(
            user=shared_owner, galaxy_server=self.server,
            api_key='shared-key', anonymous=True)

        with patch('bioblend.galaxy.histories.HistoryClient.update_history',
                   return_value={}) as update_history:
            history = WorkspaceHistory.objects.create(
                history='mine1', name='Analyse mine1', email='',
                monitored=True, finished=True, source_ip='127.0.0.1',
                workflow_category='OneClick', workflow_steps='',
                galaxy_server=self.server, user=self.user)
            history.name = 'Renamed mine1'
            history.save()

        update_history.assert_called_once()

    def test_stale_personal_galaxyuser_with_no_api_key_is_ignored(self):
        # Caught live against the real local dev DB: two leftover
        # personal GalaxyUser rows (blank api_key) from the old, now-
        # removed connection_galaxy code path that used to
        # get_or_create() one for every authenticated visitor. A plain
        # "does a personal row exist" check picks that stale, keyless
        # row and crashes with ValueError('API key must be set')
        # instead of ever reaching the shared-key fallback.
        GalaxyUser.objects.create(
            user=self.user, galaxy_server=self.server, api_key='')
        shared_owner = User.objects.create_user(username='ngphylo-shared2')
        GalaxyUser.objects.create(
            user=shared_owner, galaxy_server=self.server,
            api_key='shared-key', anonymous=True)

        with patch('bioblend.galaxy.histories.HistoryClient.update_history',
                   return_value={}) as update_history:
            history = WorkspaceHistory.objects.create(
                history='mine2', name='Analyse mine2', email='',
                monitored=True, finished=True, source_ip='127.0.0.1',
                workflow_category='OneClick', workflow_steps='',
                galaxy_server=self.server, user=self.user)
            history.name = 'Renamed mine2'
            history.save()

        update_history.assert_called_once()


class GalaxyGetRequestCacheTest(TestCase):
    """
    galaxylib.get_request_cache - stored under response.url but looked up
    by a URL containing the API key, which stopped matching once bioblend
    moved the key to a header, so it never hit. Goes through the real
    Django cache, so each hit is also a real pickle round trip.
    """
    BASE = 'http://fake-galaxy.example.org/api'

    @staticmethod
    def _response(body=b'{"id": "w1"}', status_code=200):
        response = requests.Response()
        response.status_code = status_code
        response._content = body
        response.url = 'http://fake-galaxy.example.org/no-key-in-this-url'
        return response

    def _get_twice(self, url, response=None, gi=None, **kwargs):
        gi = gi or GalaxyInstance('http://fake-galaxy.example.org', key='k1')
        with patch('bioblend.galaxyclient.GalaxyClient.make_get_request',
                   return_value=(self._response() if response is None
                                 else response)) as get:
            gi.make_get_request(url, **kwargs)
            second = gi.make_get_request(url, **kwargs)
        return get.call_count, second

    def test_workflow_definition_is_served_from_cache(self):
        calls, second = self._get_twice(self.BASE + '/workflows/w1')
        self.assertEqual(calls, 1)
        self.assertEqual(second.json(), {'id': 'w1'})

    def test_tool_definition_is_served_from_cache(self):
        calls, _ = self._get_twice(self.BASE + '/tools/mafft')
        self.assertEqual(calls, 1)
        calls, _ = self._get_twice(
            self.BASE + '/tools/toolshed.example.org/repos/o/mafft/1.0')
        self.assertEqual(calls, 1)

    def test_live_state_is_never_cached(self):
        for path in ('/histories/h1', '/histories/h1/contents',
                     '/datasets/d1', '/jobs/j1', '/workflows',
                     '/tools/mafft/build'):
            calls, _ = self._get_twice(self.BASE + path)
            self.assertEqual(calls, 2, path)

    def test_maintenance_page_is_never_cached(self):
        calls, _ = self._get_twice(
            self.BASE + '/workflows/w1', self._response(b'<html>Maintenance'))
        self.assertEqual(calls, 2)

    def test_error_response_is_never_cached(self):
        calls, _ = self._get_twice(
            self.BASE + '/workflows/w1',
            self._response(b'{"err": 1}', status_code=403))
        self.assertEqual(calls, 2)

    def test_nocache_instances_always_reach_galaxy(self):
        gi = GalaxyInstance('http://fake-galaxy.example.org', key='k1')
        gi.nocache = True
        calls, _ = self._get_twice(self.BASE + '/workflows/w1', gi=gi)
        self.assertEqual(calls, 2)

    def test_entries_are_isolated_by_api_key_and_params(self):
        url = self.BASE + '/tools/mafft'
        gi1 = GalaxyInstance('http://fake-galaxy.example.org', key='k1')
        gi2 = GalaxyInstance('http://fake-galaxy.example.org', key='k2')
        with patch('bioblend.galaxyclient.GalaxyClient.make_get_request',
                   return_value=self._response()) as get:
            gi1.make_get_request(url)
            gi2.make_get_request(url)
            gi1.make_get_request(url, params={'io_details': True})
        self.assertEqual(get.call_count, 3)
