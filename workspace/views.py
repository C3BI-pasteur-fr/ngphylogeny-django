from __future__ import unicode_literals
import json

from django.http import HttpResponse
from django.urls import reverse_lazy
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import ensure_csrf_cookie
from django.views.generic import TemplateView, ListView, DeleteView, UpdateView, DetailView, View
from django.views.generic.edit import SingleObjectMixin
from bioblend.galaxy.client import ConnectionError
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect
from tasks import deletegalaxyhistory
from workflows.tasks import deletegalaxyworkflow

from galaxy.decorator import connection_galaxy
from .models import WorkspaceHistory
from tools.models import Tool
from .tasks import updateworkspacestatus

from utils import ip

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

    if request.user.is_authenticated():
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
        # Create a new galaxy history
        history_id = create_history(request, name)

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

@method_decorator(ensure_csrf_cookie, name="dispatch")
@method_decorator(connection_galaxy, name="dispatch")
class HistoryDetailView(WorkspaceHistoryObjectMixin, DetailView):
    """
        Display Galaxy like history information
    """
    template_name = 'workspace/history.html'

    
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
            dataset_provenance = gi.histories.show_dataset_provenance(
                history_id,
                data_id,
                follow=False)
            context.update({'tool_id': dataset_provenance.get("tool_id"),
                            'dataset_id': data_id})
    return HttpResponse(json.dumps(context), content_type='application/json')


@connection_galaxy
def get_dataset_citations(request, history_id):
    """
    Ajax: return citations of all tools used in the dataset
    """
    context = dict()
    refs = []
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
                refs.extend(t.citations)
            except Tool.DoesNotExist:
                pass
    except WorkspaceHistory.DoesNotExist:
        pass
    context.update({'citations': refs})
    return HttpResponse(json.dumps(context), content_type='application/json')

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
                    refs=refs+unicode(b.reference)+unicode("\n")
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
                    refs=refs+unicode(b.txt())+unicode("\n")
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
            errormessage = ds.get('misc_info')
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

