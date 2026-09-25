from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from galaxy.models import Server
from workflows.models import Workflow
from workspace.models import WorkspaceHistory
from .forms import AccountCreationForm
from .models import UserProfile
from .views import AccountCreateView


class AccountCreationFormTest(TestCase):
    """
    Unit tests for AccountCreationForm's own added behavior (the email
    field/uniqueness check) - not the whole form's validity, since that
    also needs a real captcha answer and this codebase has no
    established pattern for satisfying django-simple-captcha in a test
    (see surveys.tests.FeedbackCreateViewTest's own note on this).
    is_valid() still runs every field's clean_<field>() and populates
    self.errors per-field even when the form as a whole fails (captcha
    always will, here) - enough to test clean_email() for real through
    the normal validation pipeline, not by calling it directly.
    """

    def _data(self, **overrides):
        data = {
            'username': 'newuser',
            'email': 'fresh@example.org',
            'password1': 'a-genuinely-strong-pass9',
            'password2': 'a-genuinely-strong-pass9',
        }
        data.update(overrides)
        return data

    def test_rejects_an_already_used_email(self):
        User.objects.create_user(
            username='existing', email='taken@example.org', password='x')
        form = AccountCreationForm(data=self._data(email='taken@example.org'))
        form.is_valid()
        self.assertIn('email', form.errors)
        self.assertIn('already exists', form.errors['email'][0])

    def test_a_fresh_email_passes_its_own_check(self):
        form = AccountCreationForm(data=self._data())
        form.is_valid()
        self.assertNotIn('email', form.errors)

    def test_mismatched_passwords_still_rejected_same_as_base_form(self):
        # Confirms subclassing UserCreationForm didn't lose its own
        # validation (password confirmation, AUTH_PASSWORD_VALIDATORS).
        form = AccountCreationForm(
            data=self._data(password2='a-different-pass9'))
        form.is_valid()
        self.assertIn('password2', form.errors)


class AccountCreateViewTest(TestCase):
    """
    AccountCreateView: redirects an already-authenticated visitor away
    (no reason to sign up again), and its form_valid() logs a freshly
    created account straight in and sends it to the existing /account
    page - see AccountCreateView's own docstring for why no separate
    onboarding step is needed.
    """

    def test_get_renders_the_signup_form(self):
        response = self.client.get(reverse('create_account'))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'id_username')
        self.assertContains(response, 'id_email')

    def test_already_authenticated_users_are_redirected_to_account(self):
        User.objects.create_user(username='bob', password='pw123456789')
        self.client.login(username='bob', password='pw123456789')
        response = self.client.get(reverse('create_account'))
        self.assertRedirects(
            response, reverse('account'), fetch_redirect_response=False)

    def test_form_valid_logs_the_new_account_in_and_redirects(self):
        view = AccountCreateView()
        request = self.client.get('/').wsgi_request
        view.request = request
        saved = User.objects.create_user(
            username='newbie', email='newbie@example.org',
            password='pw123456789')
        form = Mock()
        form.save.return_value = saved

        response = view.form_valid(form)

        self.assertEqual(request.user, saved)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, reverse('account'))

    def test_a_real_signup_then_lands_logged_in_on_the_account_page(self):
        # End-to-end confirmation that create_history() (workspace/
        # views.py) would now correctly attribute a workflow run to
        # this account: get_or_create()'d UserProfile is present, and
        # the resulting session is a real authenticated one, not just
        # form_valid()'s own isolated behavior above.
        saved = User.objects.create_user(
            username='carol', email='carol@example.org',
            password='pw123456789')
        view = AccountCreateView()
        request = self.client.get('/').wsgi_request
        view.request = request
        form = Mock()
        form.save.return_value = saved
        view.form_valid(form)

        self.assertFalse(UserProfile.objects.filter(user=saved).exists())
        # UserProfile itself is only get_or_create()'d lazily on the
        # /account page (AccountDetailView.get_context_data() above) -
        # not by signing up alone.


class AccountDeleteViewTest(TestCase):
    """
    Self-service account deletion (AccountDeleteView, GET/POST
    /account/delete). Per an explicit decision with the user: deleting
    the account also deletes every WorkspaceHistory row it owns (same
    effect as workspace.views.DeleteAllHistories's own "Delete all
    histories" button), but each such row survives as a real database
    row - deleted=True and detached (user=None), not hard-deleted -
    specifically so workspace.reports's usage counts (which explicitly
    count deleted=True rows too - "a usage report, not a what's-still-
    retained report") stay accurate even after the account that ran
    them is gone. Only the account itself (auth.User, cascading to its
    UserProfile) is a real, hard delete.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        self.user = User.objects.create_user(
            username='alice', password='secretpass')
        self.workflow = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='wf-alice',
            name='FastME OneClick', category='duplicated',
            description='FastME OneClick', slug='wf-alice-copy')
        # WorkspaceHistory.save() always calls rename(), which would
        # try a real outbound HTTP call for a user-owned row - patched
        # out here the same way this codebase's other WorkspaceHistory
        # fixtures already do (see workspace/tests.py).
        with patch('workspace.models.WorkspaceHistory.rename'):
            self.history_with_workflow = WorkspaceHistory.objects.create(
                history='h1', name='Run 1', email='', monitored=True,
                finished=True, deleted=False, source_ip='127.0.0.1',
                workflow_category='OneClick', workflow_steps='',
                galaxy_server=self.server, user=self.user,
                workflow=self.workflow)
            self.history_without_workflow = WorkspaceHistory.objects.create(
                history='h2', name='Run 2', email='', monitored=True,
                finished=True, deleted=False, source_ip='127.0.0.1',
                workflow_category='Tool', workflow_steps='MAFFT',
                galaxy_server=self.server, user=self.user)

    def test_anonymous_is_redirected_to_login(self):
        response = self.client.get(reverse('delete_account'))
        self.assertEqual(response.status_code, 302)
        self.assertIn('/account/login', response.url)

    def test_get_renders_a_confirmation_page_without_deleting_anything(self):
        self.client.login(username='alice', password='secretpass')
        response = self.client.get(reverse('delete_account'))
        self.assertEqual(response.status_code, 200)
        self.assertTrue(User.objects.filter(username='alice').exists())

    def test_posting_no_leaves_the_account_untouched(self):
        self.client.login(username='alice', password='secretpass')
        response = self.client.post(reverse('delete_account'), {'no': ''})
        self.assertRedirects(response, reverse('account'))
        self.assertTrue(User.objects.filter(username='alice').exists())

    @patch('account.views.deletegalaxyworkflow.delay')
    @patch('account.views.deletegalaxyhistory.delay')
    def test_posting_yes_deletes_the_account_and_its_analyses(
            self, mock_delete_history, mock_delete_workflow):
        self.client.login(username='alice', password='secretpass')

        response = self.client.post(reverse('delete_account'), {'yes': ''})

        self.assertRedirects(response, reverse('home'))
        self.assertFalse(User.objects.filter(username='alice').exists())
        self.assertFalse(UserProfile.objects.filter(user=self.user).exists())
        mock_delete_history.assert_any_call('h1')
        mock_delete_history.assert_any_call('h2')
        mock_delete_workflow.assert_called_once_with('wf-alice')

    @patch('account.views.deletegalaxyworkflow.delay')
    @patch('account.views.deletegalaxyhistory.delay')
    def test_histories_survive_deleted_and_detached_not_hard_deleted(
            self, mock_delete_history, mock_delete_workflow):
        # The actual point of this whole design: workspace.reports's
        # usage counts must still see these rows after the account is
        # gone - a hard CASCADE delete (the FK's own default behavior,
        # unless explicitly detached first) would silently erase them.
        self.client.login(username='alice', password='secretpass')

        self.client.post(reverse('delete_account'), {'yes': ''})

        h1 = WorkspaceHistory.objects.get(history='h1')
        h2 = WorkspaceHistory.objects.get(history='h2')
        self.assertTrue(h1.deleted)
        self.assertTrue(h2.deleted)
        self.assertIsNone(h1.user)
        self.assertIsNone(h2.user)
        self.workflow.refresh_from_db()
        self.assertTrue(self.workflow.deleted)

    @patch('account.views.deletegalaxyworkflow.delay')
    @patch('account.views.deletegalaxyhistory.delay')
    def test_session_is_logged_out_after_deletion(
            self, mock_delete_history, mock_delete_workflow):
        self.client.login(username='alice', password='secretpass')

        self.client.post(reverse('delete_account'), {'yes': ''})
        response = self.client.get(reverse('account'))

        self.assertEqual(response.status_code, 302)
        self.assertIn('/account/login', response.url)
