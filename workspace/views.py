from __future__ import unicode_literals
import json

import requests
from django.contrib import messages
from django.contrib.admin.views.decorators import staff_member_required
from django.core import signing
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from django.utils import timezone
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import ensure_csrf_cookie
from django.views.generic import TemplateView, ListView, DeleteView, UpdateView, DetailView, View
from django.views.generic.edit import SingleObjectMixin
from bioblend.galaxy.client import ConnectionError
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect
from blast.models import BlastRun
from .tasks import deletegalaxyhistory
from workflows.tasks import deletegalaxyworkflow

from galaxy.decorator import connection_galaxy
from .emails import site_url
from .models import WorkspaceHistory
from .reports import CATEGORY_LABELS, build_report_web_context
from tools.models import Tool
from .tasks import updateworkspacestatus

# Salt for the workspace/histories permalink (PreviousHistoryListView/
# WorkspacePermalinkView below) - just namespaces the signed token so it
# can't be reinterpreted by some unrelated future use of
# django.core.signing in this app; doesn't add real secrecy on its own
# (nothing here does - see WorkspacePermalinkView's own docstring for
# why that's fine).
PERMALINK_SALT = 'workspace.permalink'

from utils import ip


@staff_member_required
def daily_report_view(request):
    """
    Web view of the same daily workflow-usage report emailed by
    workspace.tasks.send_daily_report (see workspace/reports.py) - admin/
    staff only (@staff_member_required redirects to the admin login page,
    same as the Django admin itself, for anyone not logged in as staff).

    Cached (see build_report_web_context()) - ?refresh=1 bypasses that
    for anyone who wants to force an immediate up-to-date render.
    """
    force_refresh = request.GET.get('refresh') == '1'
    return render(request, 'workspace/report_page.html',
                  build_report_web_context(force_refresh=force_refresh))


@staff_member_required
def running_jobs_view(request):
    """
    Live admin/staff-only view of everything currently running - both
    regular workflow/tool runs (WorkspaceHistory: monitored, not yet
    finished/deleted) and BLAST runs (BlastRun: still PENDING/RUNNING,
    not deleted) - oldest first, so a run approaching or past one of the
    two staleness cutoffs (workspace.tasks.WORKFLOW_RUN_STALE_AFTER,
    blast.tasks.PASTEUR_RUN_STALE_AFTER - both cancel/error a run out
    once it's been going too long, see CLAUDE.md) surfaces at the top
    instead of being buried under newer ones. Not cached, unlike
    daily_report_view - the whole point here is a live, up-to-the-moment
    picture, not a 15-minute-old snapshot.
    """
    now = timezone.now()
    rows = []

    for w in (WorkspaceHistory.objects
              .filter(monitored=True, finished=False, deleted=False)
              .values('history', 'name', 'created_date', 'email',
                      'workflow_category', 'workflow_steps', 'workflow__name',
                      'history_content_json')):
        # history_content_json is only ever dict-shaped in the normal
        # case - a genuinely malformed one shouldn't crash this page,
        # just show as 0 datasets (see WorkspaceHistoryObjectMixin.
        # get_context_data's own note on this - hit live once already).
        try:
            content = json.loads(w['history_content_json'] or '[]')
        except ValueError:
            content = []
        steps = [f for f in content if isinstance(f, dict)]
        done = sum(1 for f in steps if 'ok' in f.get('state', ''))
        running = sum(
            1 for f in steps
            if any(s in f.get('state', '') for s in ('new', 'queued', 'running')))
        rows.append({
            'kind': 'Workflow',
            'name': w['name'],
            'type': CATEGORY_LABELS.get(
                w['workflow_category'], w['workflow_category'] or 'Unknown'),
            'email': w['email'],
            'created_date': w['created_date'],
            'runtime': now - w['created_date'],
            'steps_total': len(steps),
            'steps_done': done,
            'steps_running': running,
            'url': reverse('history_detail', kwargs={'history_id': w['history']}),
        })

    server_labels = dict(BlastRun.BLASTSERVERS)
    status_labels = dict(BlastRun.RUNSTATUS)
    for b in (BlastRun.objects
              .filter(status__in=[BlastRun.PENDING, BlastRun.RUNNING], deleted=False)
              .values('id', 'query_id', 'date', 'email', 'server', 'blastprog', 'status')):
        is_running = b['status'] == BlastRun.RUNNING
        rows.append({
            'kind': 'BLAST',
            'name': b['query_id'] or 'BLAST run',
            'type': '%s BLAST (%s, %s)' % (
                server_labels.get(b['server'], b['server']), b['blastprog'],
                status_labels.get(b['status'], b['status'])),
            'email': b['email'],
            'created_date': b['date'],
            'runtime': now - b['date'],
            'steps_total': 1,
            'steps_done': 0,
            'steps_running': 1 if is_running else 0,
            'url': reverse('blast_view', kwargs={'pk': b['id']}),
        })

    rows.sort(key=lambda r: r['created_date'])

    return render(request, 'workspace/running_jobs.html', {'rows': rows, 'now': now})


@connection_galaxy
def create_history(request, name='', wf_category='', wf_steps=''):
    """
    Create a new galaxy history

    :param request:
    :param name: name of new history
    :return: galaxy id history
    """
    gi = request.galaxy
    server = request.galaxy_server
    if not name:
        name = 'NGPhylogeny analyse'
    history = gi.histories.create_history(name=name)

    if request.user.is_authenticated:
        current_user = request.user
    else:
        current_user = server.galaxyuser_set.get(anonymous=True).user

    wsph = WorkspaceHistory(history=history.get("id"),
                            name=history.get('name'),
                            user=current_user,
                            galaxy_server=server,
                            history_content_json = json.dumps(history),
                            history_info_json = json.dumps(history),
                            source_ip=ip.get_client_ip(request),
                            workflow_category=wf_category,
                            workflow_steps = wf_steps,
                            )
    wsph.save()

    # save the current history in session
    request.session.setdefault('histories', [])
    request.session['histories'].append(wsph.history)
    request.session["last_history"] = wsph.history
    request.session.modified = True
    return wsph


@connection_galaxy
def get_history(request):
    return request.session.get('last_history')


def get_or_create_history(request, name=''):
    """
    :param request:
    :param name: name of new history
    :return: history_id
    """
    history_id = get_history(request)
    if not history_id:
        # Create a new galaxy history - create_history() returns the
        # WorkspaceHistory model instance (other callers, e.g.
        # tools/views.py, need it for its FKs/wf_category/wf_steps), not
        # a plain id - unwrap .history here to actually satisfy this
        # function's own "return: history_id" contract.
        history_id = create_history(request, name).history

    return history_id


@connection_galaxy
def delete_history(request, history_id=None):
    """
    Delete history, update session
    :param request:
    :param history_id:
    :return:
    """

    last_history = get_history(request)

    if history_id == last_history:
        if history_id in request.session['histories']:
            request.session['histories'].remove(history_id)
            request.session.modified = True

        request.session['last_history'] = request.session['histories'][-1]

    WorkspaceHistory.objects.get(history=history_id).delete()

def resolve_dataset_tools(gi, galaxy_server, history_id, dataset_ids):
    """
    Resolve {dataset_id: tool_id} for every dataset in dataset_ids (via
    Galaxy's provenance API - there's no bulk equivalent bioblend
    exposes) and {tool_id: name} for every tool_id resolved that way.

    Tool names are looked up in this app's own Tool model (mirrors every
    tool NGPhylogeny actually runs - see CLAUDE.md's "App
    responsibilities") before ever falling back to a Galaxy API call:
    NGPhylogeny only ever submits its own preconfigured/imported
    workflows, so the tool that produced any given dataset is almost
    always already known locally, and a DB lookup is both cheaper and
    doesn't depend on Galaxy being reachable at all.

    Used to be done independently, every 10s poll, by three separate
    call sites (the step chain, the table, and the citations list - each
    re-deriving the exact same dataset->tool_id/tool_id->name mappings
    on their own) - see CLAUDE.md's step-chain section for the real
    production incident (pod restarts) this contributed to. Computed
    once per poll instead, in HistoryContentRefreshView.get_context_data,
    and shared across all three.

    A dataset/tool that can't be resolved (a transient Galaxy failure -
    see get_dataset_toolprovenance's own docstring for why that's common
    enough to matter here) is simply left out of the returned dicts
    rather than raising - callers already treat "not resolved yet" as a
    normal state (a step chain box/table row just keeps its fallback
    label).
    """
    dataset_tool_ids = {}
    for dataset_id in dataset_ids:
        try:
            provenance = gi.histories.show_dataset_provenance(
                history_id, dataset_id, follow=False)
        except (ConnectionError, requests.exceptions.RequestException):
            continue
        tool_id = provenance.get('tool_id')
        if tool_id:
            dataset_tool_ids[dataset_id] = tool_id

    tool_names = {}
    for tool_id in set(dataset_tool_ids.values()):
        local_tool = Tool.objects.filter(
            galaxy_server=galaxy_server, id_galaxy=tool_id).first()
        if local_tool:
            tool_names[tool_id] = local_tool.name
            continue
        try:
            tool = gi.tools.show_tool(tool_id=tool_id)
            tool_names[tool_id] = tool.get('name')
        except (ConnectionError, requests.exceptions.RequestException):
            pass

    return dataset_tool_ids, tool_names


def build_citations(dataset_tool_ids):
    """
    Citation list for every distinct tool in dataset_tool_ids.values() -
    shared logic between get_dataset_citations (the standalone AJAX
    endpoint, kept for any other caller) and HistoryContentRefreshView,
    which now computes this itself alongside the step chain/table so the
    history detail page's 10s poll doesn't also re-derive it separately
    via its own extra Galaxy calls (see resolve_dataset_tools above).
    """
    refs = [ngphylo_citation()]
    for tool_id in set(dataset_tool_ids.values()):
        try:
            t = Tool.objects.get(id_galaxy=tool_id)
            refs.extend(t.citations)
        except Tool.DoesNotExist:
            pass
    return refs


class WorkspaceHistoryObjectMixin(SingleObjectMixin):
    model = WorkspaceHistory
    pk_url_kwarg = 'history_id'

    def get_object(self, queryset=None):
        if queryset is None:
            queryset = self.get_queryset()
        server = self.request.galaxy_server
        hist_id = self.kwargs.get(self.pk_url_kwarg)
        if not hist_id:
            hist_id = self.request.session["last_history"]
        updateworkspacestatus.delay(hist_id)
        w = queryset.get(history=hist_id,
                         galaxy_server=server)

        w.history_content = json.loads(w.history_content_json)
        w.history_info = json.loads(w.history_info_json)
        return w

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Shared by every view built on this mixin - both HistoryDetailView
        # (the very first render) and HistoryContentRefreshView (every
        # later poll) include workspace/include/history_contents_
        # refreshable.html, which needs this same data either way. This
        # used to live only on HistoryContentRefreshView's own
        # get_context_data() - meaning an already-finished history (whose
        # client-side polling timer gets cleared before ever firing once
        # - see that template's own {% if object.finished %} handling)
        # loaded straight through HistoryDetailView and got a step chain/
        # table/citations list that stayed permanently unresolved, since
        # /refresh was never going to be called to backfill it. Real
        # production bug, not just a theoretical gap - a completed run's
        # tool column and citations list stayed blank forever.
        #
        # Matches the template's own "is there a real step chain/table
        # yet" condition (that same file's top-level {% if %}) - nothing
        # to resolve in the "please wait" state, so skip the Galaxy calls
        # entirely rather than doing pointless work before a run has even
        # really started.
        history_content = self.object.history_content or []
        if len(history_content) > 1:
            # Real production crash: history_content is only ever
            # dict-shaped in the *normal* case - a genuinely malformed/
            # unexpected history_content_json (seen live: a plain list
            # of strings) crashed this with AttributeError: 'str' object
            # has no attribute 'get', 500ing the whole page. Hit right on
            # a new workflow submission - Galaxy is at its busiest
            # scheduling many jobs at once right then, a plausible moment
            # for show_history(contents=True) to transiently return
            # something other than its usual dataset list (workspace.
            # tasks.initializeworkspacejob/updateworkspacestatus store
            # whatever it returns verbatim, with no validation). The
            # template rendering this exact same data
            # (dictsortreversed:"hid" etc.) never crashed on it, since
            # Django's template engine fails a bad attribute lookup
            # silently rather than raising - this filters out anything
            # that isn't actually dict-shaped (or has no 'id') to match
            # that same tolerance, rather than assuming the shape.
            dataset_ids = [f.get('id') for f in history_content
                           if isinstance(f, dict) and f.get('id')]
            dataset_tool_ids, tool_names = resolve_dataset_tools(
                self.request.galaxy, self.request.galaxy_server,
                self.object.history_info['id'], dataset_ids)
            context['dataset_tool_ids'] = dataset_tool_ids
            context['tool_names'] = tool_names
            context['citations'] = build_citations(dataset_tool_ids)

        # First finished newick/nhx dataset in the history, if any - the
        # same file the table's own {% if file.extension in "nhx,nwk" %}
        # branch shows the "Interactive Tree visualisation"/iTOL buttons
        # for (see history_contents_refreshable.html). Read straight off
        # history_content, already fetched above with no extra Galaxy
        # call - this is what history_contents_refreshable.html's inline
        # phylotree.js preview fetches (via display_raw) and renders,
        # once, the first time it appears.
        tree_dataset = next(
            (f for f in history_content
             if isinstance(f, dict) and f.get('extension') in ('nhx', 'nwk')
             and 'ok' in (f.get('state') or '')),
            None)
        if tree_dataset and tree_dataset.get('id'):
            context['tree_preview_dataset_id'] = tree_dataset['id']
            context['tree_preview_url'] = reverse(
                'display_raw', kwargs={'file_id': tree_dataset['id']})
        return context


@method_decorator(ensure_csrf_cookie, name="dispatch")
@method_decorator(connection_galaxy, name="dispatch")
class HistoryDetailView(WorkspaceHistoryObjectMixin, DetailView):
    """
        Display Galaxy like history information
    """
    template_name = 'workspace/history.html'


@method_decorator(connection_galaxy, name="dispatch")
class HistoryContentRefreshView(WorkspaceHistoryObjectMixin, DetailView):
    """
    Renders just the part of the history detail page that actually
    changes as a run progresses (workspace/include/
    history_contents_refreshable.html) - either the step-chain + dataset
    table once there are 2+ datasets, or a "please wait" message before
    that (that template's own top comment has the detail; this used to
    be a second, entirely separate template/timer - history_wait.html,
    now deleted - unified into this same one so both states share the
    exact same shell/polling). Polled client-side (jQuery's .load(),
    which - unlike a plain AJAX GET swapped in via .html() - executes
    the <script> tags in the response, so this reuses the exact same
    in-template JS the initial page load already runs, no separate
    client-side re-init logic needed) instead of the full page reload
    templates/workspace/history.html used to do every 10s.
    No @ensure_csrf_cookie here (unlike HistoryDetailView) - the initial
    full page load already guarantees the cookie exists by the time this
    is ever called from that page's own JS.

    ?staging=1 tells the template it's being loaded into the hidden
    #history-refreshable-staging container rather than rendered directly
    into the live #history-refreshable-region - see the template's own
    comments for why (builds the step chain/table out of sight, then
    swaps the finished result in, instead of visibly blanking and
    rebuilding the live page on every poll).
    """
    template_name = 'workspace/include/history_contents_refreshable.html'

    def get_context_data(self, **kwargs):
        # dataset_tool_ids/tool_names/citations are resolved by the
        # shared WorkspaceHistoryObjectMixin.get_context_data() above -
        # this just adds the one thing specific to being polled rather
        # than directly included on first load.
        context = super().get_context_data(**kwargs)
        context['staging'] = self.request.GET.get('staging') == '1'
        return context


@connection_galaxy
def get_dataset_toolprovenance(request, history_id, ):
    """
    Ajax: return tool id who produced the dataset
    """
    context = dict()
    if request.POST:
        gi = request.galaxy

        data_id = request.POST.get('dataset_id')
        if data_id:
            try:
                dataset_provenance = gi.histories.show_dataset_provenance(
                    history_id,
                    data_id,
                    follow=False)
                context.update({'tool_id': dataset_provenance.get("tool_id"),
                                'dataset_id': data_id})
            except (ConnectionError, requests.exceptions.RequestException):
                # Galaxy (or a proxy in front of it) can be transiently
                # unreachable/return a 502 - this endpoint is now polled
                # frequently (once per dataset, every 10s - see the
                # history detail page's step chain/table), so a single
                # hiccup shouldn't turn into an unhandled Django 500 (and
                # an admin error email under DEBUG=False) on every failed
                # poll. Both exception types matter here, not just
                # bioblend's own ConnectionError: bioblend's own retry
                # logic (bioblend.galaxy.client.Client._get) only catches
                # requests.exceptions.ConnectionError itself (converting
                # it to this one) - a plain read timeout
                # (requests.exceptions.ReadTimeout, the likely shape of a
                # slow/overloaded rather than fully unreachable Galaxy,
                # and specifically what GalaxyUser.get_galaxy_instance's
                # new timeout=30 is meant to turn a stuck request into)
                # isn't a ConnectionError subclass and would otherwise
                # propagate straight through uncaught. The client already
                # treats a failed request the same as any other rejected
                # promise - see history_contents_refreshable.html's own
                # .then(success, fail) handling - a non-2xx status is
                # enough for that.
                return HttpResponse(
                    json.dumps({'dataset_id': data_id}),
                    content_type='application/json', status=502)
    return HttpResponse(json.dumps(context), content_type='application/json')


@connection_galaxy
def get_dataset_citations(request, history_id):
    """
    Ajax: return citations of all tools used in the dataset
    """
    context = dict()
    refs = [ngphylo_citation()]
    tools = []
    gi = request.galaxy
    try:
        w = WorkspaceHistory.objects.get(history=history_id)
        w.history_content = json.loads(w.history_content_json)
        for file in w.history_content:
            try:
                dataset_provenance = gi.histories.show_dataset_provenance(
                    history_id,
                    file.get('id'),
                    follow=False)
            except (ConnectionError, requests.exceptions.RequestException):
                # Same transient-Galaxy-failure reasoning as
                # get_dataset_toolprovenance above (including why both
                # exception types are caught, not just bioblend's own
                # ConnectionError) - this makes one bioblend call per
                # dataset in the history, so it's even more likely than
                # that one to hit a flaky/overloaded Galaxy on a big
                # history. Skip just this dataset's provenance rather
                # than failing the whole citations fetch (and crashing
                # this now-every-10s-polled endpoint) over one bad
                # dataset.
                continue
            tools.append(dataset_provenance.get('tool_id'))
        tools = list(set(tools))
        for tid in tools:
            try:
                t = Tool.objects.get(id_galaxy=tid)
                refs.extend(t.citations)
            except Tool.DoesNotExist:
                pass
    except WorkspaceHistory.DoesNotExist:
        pass
    context.update({'citations': refs})
    return HttpResponse(json.dumps(context), content_type='application/json')


def ngphylo_citation():
    ref = """
    <a target="_blank" href="https://doi.org/10.1093/nar/gkz303">
    <div class="pub-date">2019</div>
    <div class="pub-content clear">
    <span class="pub-author-list">Lemoine, F. and Correia, D. and Lefort, V. and Doppelt-Azeroual, O. and Mareuil, F. and Cohen-Boulakia, S. and Gascuel, O.</span>
    <span class="pub-title">NGPhylogeny.fr: new generation phylogenetic services for non-specialists.</span>
    <span class="pub-journal-name">Nucleic acids research, 47:W260-W265</span>
    </div>
    </a>
    """
    return ref


@connection_galaxy
def get_dataset_citations_bibtex(request, history_id):
    """
    Ajax: return citations of all tools used in the dataset
    """
    refs = ""
    tools = []
    gi = request.galaxy
    try:
        w = WorkspaceHistory.objects.get(history=history_id)
        w.history_content = json.loads(w.history_content_json)
        for file in w.history_content:
            dataset_provenance = gi.histories.show_dataset_provenance(
                history_id,
                file.get('id'),
                follow=False)
            tools.append(dataset_provenance.get('tool_id'))
        tools = list(set(tools))
        for tid in tools:
            try:
                t = Tool.objects.get(id_galaxy=tid)
                for b in t.citation_set.all():
                    refs=refs+b.reference+"\n"
            except Tool.DoesNotExist:
                pass
    except WorkspaceHistory.DoesNotExist:
        pass
    return HttpResponse(refs, content_type='text/plain; charset=utf-8')

@connection_galaxy
def get_dataset_citations_txt(request, history_id):
    """
    Ajax: return citations of all tools used in the dataset
    """
    refs = ""
    tools = []
    gi = request.galaxy
    try:
        w = WorkspaceHistory.objects.get(history=history_id)
        w.history_content = json.loads(w.history_content_json)
        for file in w.history_content:
            dataset_provenance = gi.histories.show_dataset_provenance(
                history_id,
                file.get('id'),
                follow=False)
            tools.append(dataset_provenance.get('tool_id'))
        tools = list(set(tools))
        for tid in tools:
            try:
                t = Tool.objects.get(id_galaxy=tid)
                for b in t.citation_set.all():
                    refs=refs+b.txt()+"\n"
            except Tool.DoesNotExist:
                pass
    except WorkspaceHistory.DoesNotExist:
        pass
    return HttpResponse(refs, content_type='text/plain; charset=utf-8')


@method_decorator(connection_galaxy, name="dispatch")
class GalaxyErrorView(TemplateView):
    """
    Redirect to Galaxy server error page
    """
    template_name='display_galaxyerror.html'
    def get_context_data(self, *args, **kwargs):
        gi = self.request.galaxy
        context = super(GalaxyErrorView, self).get_context_data(*args, **kwargs)
        dsid = kwargs.get('id')
        ds = gi.datasets.show_dataset(dsid)
        state = ''
        name = ''
        errormessage = ''
        hid=''
        if ds is not None:
            state = ds.get('state')
            #errormessage = ds.get('misc_info')
            job = ds.get('creating_job')
            jinfo = gi.jobs.show_job(job,full_details=True)
            errormessage = jinfo.get('stderr')+jinfo.get('stdout')
            hid = ds.get('history_id')
            name = ds.get('name')
        context['state'] = state
        context['error'] = errormessage
        context['history_id'] = hid
        context['jobname'] = name
        return context

@method_decorator(connection_galaxy, name="dispatch")
class PreviousHistoryListView(ListView):
    """
    Display list of Previous analyses stored in the sessions cookies
    """
    queryset = WorkspaceHistory.objects.none()
    template_name = 'workspace/previous_analyses.html'
    context_object_name = 'histories'

    def get_queryset(self):
        self.queryset = WorkspaceHistory.objects.filter(history__in=self.request.session.get('histories', [])).filter(deleted=False).order_by("-created_date")

        # update session history
        self.request.session['histories'] = list(self.queryset.values_list('history', flat=True))

        return self.queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Permalink to *this* list - see WorkspacePermalinkView below.
        # Built from the just-cleaned session list (get_queryset() above
        # already dropped any deleted history), not the raw queryset, so
        # the token always matches what get_queryset() would itself
        # filter down to a second time when the permalink is followed.
        history_ids = self.request.session.get('histories', [])
        context['permalink_url'] = (
            site_url(reverse('workspace_permalink', kwargs={
                'token': signing.dumps(history_ids, salt=PERMALINK_SALT),
            }))
            if history_ids else None)
        return context


class WorkspacePermalinkView(View):
    """
    GET /workspace/permalink/<token> - restores the exact list of
    analyses a permalink (PreviousHistoryListView.get_context_data()
    above) was generated for into *this* browser's session, then
    redirects to the normal previous-analyses page. The whole point:
    "Workspace" is otherwise only ever readable from the one browser
    session that actually ran each analysis (see PreviousHistoryListView's
    own docstring) - clearing cookies, switching devices, or just coming
    back much later loses access to it entirely even though the
    underlying data is still there. This link is how to get it back (or
    hand the same list to a collaborator).

    The token itself is just a signed (not encrypted) list of history
    ids - django.core.signing guarantees it wasn't tampered with, not
    that it's secret. That's consistent with how an individual history
    is already reachable: WorkspaceHistoryObjectMixin.get_object()
    (history_detail, /workspace/history/<id>) has no ownership/session
    check at all - anyone who knows a 16-character Galaxy history id can
    already open that history directly. A permalink bundling several of
    those already-not-secret ids together doesn't introduce a new kind
    of exposure, just a convenient way to share/restore the same list.

    Merges with (rather than replacing) whatever's already in the
    current session, so following a permalink in a browser that already
    has its own, different set of analyses adds to that list instead of
    losing it.
    """

    def get(self, request, token):
        try:
            history_ids = signing.loads(token, salt=PERMALINK_SALT)
        except signing.BadSignature:
            messages.add_message(
                request, messages.ERROR,
                "This permalink is invalid or has been corrupted.")
            return redirect('previous_analyses')

        existing = set(request.session.get('histories', []))
        request.session['histories'] = list(existing.union(history_ids))
        return redirect('previous_analyses')


@method_decorator(connection_galaxy, name="dispatch")
class WorkspaceDeleteView(WorkspaceHistoryObjectMixin, DeleteView):
    """
    Delete Workspace and history
    """
    success_url = reverse_lazy('previous_analyses')

    # Overrides the DeletionMixin delete method to prevent actual
    # deletion, but instead mark it as deleted and delete the
   # galaxy history asyynchronously
    def delete(self, request, *args, **kwargs):
        """
        Calls the delete() method on the fetched object and then
        redirects to the success URL.
        """
        self.object = self.get_object()
        self.object.deleted = True
        self.object.save()
        if self.object.workflow is not None:
            deletegalaxyworkflow.delay(self.object.workflow.id_galaxy)
            self.object.workflow.deleted = True
            self.object.workflow.save()
        deletegalaxyhistory.delay(self.object.history)
        return HttpResponseRedirect(self.get_success_url())

class DeleteAllHistories(View):
    template_name = 'workspace/delete_all_histories_confirm.html'
    
    def get(self,request):
        return render(request, self.template_name, {})

    def post(self,request):
        if "yes" in request.POST:
            for e in WorkspaceHistory.objects.filter(history__in=request.session.get('histories', [])).filter(deleted=False):
                if e.workflow is not None:
                    deletegalaxyworkflow.delay(e.workflow.id_galaxy)
                    e.workflow.deleted = True
                    e.workflow.save()
                e.deleted = True
                e.save()
                deletegalaxyhistory.delay(e.history)
        return redirect('previous_analyses')

@method_decorator(connection_galaxy, name="dispatch")
class WorkspaceRenameView(HistoryDetailView, UpdateView):
    """
    Rename Workspace
    """
    fields = ['name']
    template_name = 'workspace/history.html'

    def get_success_url(self):
        return reverse_lazy('history_detail', args=(self.get_object().history,))

    
@method_decorator(connection_galaxy, name="dispatch")
class WorkspaceChangeEmailView(HistoryDetailView, UpdateView):
    """
    Change workspace contact Email
    """
    model = WorkspaceHistory
    fields = ['email']
    template_name = 'workspace/history.html'

    def get_success_url(self):
        return reverse_lazy('history_detail', args=(self.get_object().history,))

