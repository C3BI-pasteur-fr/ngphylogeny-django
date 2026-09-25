from django.conf import settings
from django.contrib import messages
from django.contrib.auth import login, logout
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.auth.views import LoginView
from django.shortcuts import redirect, render
from django.urls import reverse_lazy
from django.views.generic import CreateView, TemplateView, View

from galaxy.models import Server
from workflows.tasks import deletegalaxyworkflow
from workspace.models import WorkspaceHistory
from workspace.tasks import deletegalaxyhistory
from .forms import AccountCreationForm
from .models import UserProfile


class AccountLoginView(LoginView):
    """
    Plain django.contrib.auth.views.LoginView, except it also exposes
    settings.NGPHYLO_ACCOUNT_CREATION_ENABLED to the template
    (account/login.html's own "Create an account" link) - the stock
    LoginView has no context hook for this, and the link needs to
    reflect the real, current setting rather than a manually-maintained
    commented-out block that only updates when someone remembers to
    edit the template by hand.
    """
    template_name = 'account/login.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['account_creation_enabled'] = \
            settings.NGPHYLO_ACCOUNT_CREATION_ENABLED
        return context


class AccountCreateView(CreateView):
    """
    Public account sign-up (GET/POST /account/create) - see
    account/forms.py's AccountCreationForm for the actual validation.
    Logs the new account straight in and sends them to the existing
    /account page (name="account", AccountDetailView below).

    Gated off by default (settings.NGPHYLO_ACCOUNT_CREATION_ENABLED,
    driven by the NGPHYLO_ACCOUNT_CREATION_ENABLED env var /
    ACCOUNT_CREATION_ENABLED GitLab CI/CD variable - see settings/
    base.py) pending a real RGPD/privacy notice: this form collects
    personal data (email) with no consent checkbox or privacy-policy
    text anywhere in the app yet. Same code-level "quick disable" shape
    already established for BLAST (see CLAUDE.md's "BLAST analysis was
    briefly, temporarily disabled" section and templates/blast/
    blast_disabled.html), just settings-driven rather than a bare
    module constant, so it can be flipped per-deployment without a code
    change. Nothing else about the feature is removed -
    AccountCreationForm/the URL/form_valid() all still work exactly as
    before; this only gates dispatch() itself. The login page's own
    "Create an account" link (AccountLoginView below) reads the same
    setting, so it appears/disappears automatically alongside this.
    """
    form_class = AccountCreationForm
    template_name = 'account/create_account.html'
    success_url = reverse_lazy('account')

    def dispatch(self, request, *args, **kwargs):
        if not settings.NGPHYLO_ACCOUNT_CREATION_ENABLED:
            return render(
                request, 'account/create_account_disabled.html', status=503)
        if request.user.is_authenticated:
            return redirect('account')
        return super().dispatch(request, *args, **kwargs)

    def form_valid(self, form):
        response = super().form_valid(form)
        # login() needs to know which backend authenticated this user -
        # normally set by authenticate(), which isn't called here since
        # the account was just created directly (its password is
        # already known to be correct, nothing to authenticate against).
        # ModelBackend is this project's only configured backend
        # (AUTHENTICATION_BACKENDS isn't overridden anywhere in
        # settings/, so it's Django's own default).
        login(self.request, self.object,
              backend='django.contrib.auth.backends.ModelBackend')
        return response


class AccountDetailView(LoginRequiredMixin, TemplateView):
    """
    /account - a logged-in NGPhylogeny account's own page. Used to be
    galaxy.views.UpdateApiKey, a form for managing a *personal* Galaxy
    API key - removed along with that whole concept (every visitor now
    authenticates to Galaxy through one single, shared key regardless
    of their NGPhylogeny account - see galaxy.decorator.
    connection_galaxy) rather than leaving a form on the page that no
    longer affects anything. Now just account info plus this account's
    own analyses (WorkspaceHistory.user already gets set correctly for
    every authenticated submission by workspace.views.create_history() -
    nothing needed there).
    """
    template_name = 'account/user_info.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['profile'] = UserProfile.objects.get_or_create(
            user=self.request.user)[0]
        context['galaxy_server'] = Server.objects.filter(
            current=True).first()
        context['histories'] = (
            WorkspaceHistory.objects
            .filter(user=self.request.user, deleted=False)
            .order_by('-created_date'))
        return context


class AccountDeleteView(LoginRequiredMixin, View):
    """
    Self-service "delete my account" (GET/POST /account/delete). A plain
    View with its own get()/post() - not Django's DeleteView, which has
    a real Django 4.x quirk already hit and fixed once in this codebase
    for BLAST (see CLAUDE.md's "Deleting a BLAST run" section:
    BaseDeleteView.get() no longer shows a real confirmation page the
    way callers here expect, and its post()/form_valid() calls
    self.object.delete() directly, bypassing any custom delete() at
    all). Mirrors workspace.views.DeleteAllHistories's own simple
    get()-renders-a-confirm-page / post()-checks-"yes" shape instead,
    already established and working in this codebase.

    Per an explicit decision with the user: deleting the account also
    deletes every one of this account's own WorkspaceHistory rows - the
    same soft-delete-and-queue-Galaxy-cleanup effect as clicking
    "Delete all histories" on the Workspace page (workspace.views.
    DeleteAllHistories), not a separate/different cleanup path.

    Critically, each such row is also detached (user set to None)
    *before* the User row itself is deleted, not just marked deleted=True
    - WorkspaceHistory.user has on_delete=CASCADE, so without this,
    deleting the User row would hard-delete every one of these rows from
    the database outright, undoing the soft-delete and permanently
    losing them from workspace.reports's usage counts. That report
    explicitly counts deleted=True rows too (see its own module
    docstring: "a usage report, not a what's-still-retained report") -
    the row surviving, deleted and unlinked, is what keeps historical
    workflow/BLAST usage stats accurate even after the account that ran
    them is gone. The account itself (auth.User, and UserProfile via its
    OneToOneField's own CASCADE) is a real, hard delete - no PII is kept
    around for reporting; only total_users/the user-growth chart
    (workspace.reports.gather_user_growth(), which counts live User rows)
    will reflect one fewer account from this point forward - the normal,
    expected shape of a "current number of accounts" metric, not a loss
    of the workflow-usage history above.
    """
    template_name = 'account/delete_account_confirm.html'

    def get(self, request):
        return render(request, self.template_name)

    def post(self, request):
        if 'yes' not in request.POST:
            return redirect('account')

        user = request.user
        for h in WorkspaceHistory.objects.filter(user=user, deleted=False):
            if h.workflow is not None:
                deletegalaxyworkflow.delay(h.workflow.id_galaxy)
                h.workflow.deleted = True
                h.workflow.save()
            deletegalaxyhistory.delay(h.history)
            h.deleted = True
            h.user = None
            h.save()

        logout(request)
        user.delete()
        messages.add_message(
            request, messages.INFO,
            "Your account and its analyses have been deleted.")
        return redirect('home')
