import json

from django.http import HttpResponse
from django.urls import reverse_lazy
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import ensure_csrf_cookie
from django.views.generic import RedirectView, ListView, DeleteView, UpdateView, DetailView
from django.views.generic.edit import SingleObjectMixin

from galaxy.decorator import connection_galaxy
from models import WorkspaceHistory


@connection_galaxy
def create_history(request, name=''):
    """
    Create a new galaxy history

    :param request:
    :param name: name of new history
    :return: galaxy id history
    """
    gi = request.galaxy
    server = request.galaxy_server
    if name:
        name = 'NGPhylogeny analyse'
    history = gi.histories.create_history(name=name)

    if request.user.is_authenticated():
        current_user = request.user
    else:
        current_user = server.galaxyuser_set.get(anonymous=True).user

    wsph = WorkspaceHistory(history=history.get("id"),
                            name=history.get('name'),
                            user=current_user,
                            galaxy_server=server
                            )
    wsph.save()

    # save the current history in session
    request.session.setdefault('histories', [])
    request.session['histories'].append(wsph.history)
    request.session["last_history"] = wsph.history
    request.session.modified = True
    return wsph.history


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
def delete_history(history_id):
    """
    Delete history
    :param history_id:
    :return:
    """

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

        return queryset.get(history=hist_id,
                            galaxy_server=server)


@method_decorator(ensure_csrf_cookie, name="dispatch")
@method_decorator(connection_galaxy, name="dispatch")
class HistoryDetailView(WorkspaceHistoryObjectMixin, DetailView):
    """
        Display Galaxy like history information
    """
    template_name = 'workspace/history.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(HistoryDetailView, self).get_context_data(**kwargs)

        # first display history with id history contained in the url
        history_id = context.get('history_id', None) or self.kwargs.get('history_id', None)

        # if no history id try to retrieve
        if not history_id:
            history_id = get_history(self.request)

        if not history_id:
            return context

        gi = self.request.galaxy
        gi.nocache = True
        context['history_info'] = gi.histories.show_history(history_id)
        history_content = gi.histories.show_history(history_id, contents=True)

        context['history_content'] = history_content
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
            dataset_provenance = gi.histories.show_dataset_provenance(history_id,
                                                                      data_id,
                                                                      follow=False)

            context.update({'tool_id': dataset_provenance.get("tool_id"),
                            'dataset_id': data_id})

    return HttpResponse(json.dumps(context), content_type='application/json')


@method_decorator(connection_galaxy, name="dispatch")
class GalaxyErrorView(RedirectView):
    """
    Redirect to Galaxy server error page
    """

    def get_redirect_url(self, *args, **kwargs):
        return "%s/dataset/errors?id=%s" % (self.request.galaxy_server.url, kwargs.get('id'))


@method_decorator(connection_galaxy, name="dispatch")
class PreviousHistoryListView(ListView):
    """
    Display list of Previous analyses stored in the sessions cookies
    """
    queryset = WorkspaceHistory.objects.none()
    template_name = 'workspace/previous_analyses.html'
    context_object_name = 'histories'

    def get_queryset(self):
        self.queryset = WorkspaceHistory.objects.filter(history__in=self.request.session.get('histories', [])
                                                        ).order_by("-created_date")

        # update session history
        self.request.session['histories'] = list(self.queryset.values_list('history', flat=True))

        return self.queryset


@method_decorator(connection_galaxy, name="dispatch")
class WorkspaceDeleteView(WorkspaceHistoryObjectMixin, DeleteView):
    """
    Delete Workspace and history
    """
    success_url = reverse_lazy('previous_analyses')


@method_decorator(connection_galaxy, name="dispatch")
class WorkspaceRenameView(HistoryDetailView, UpdateView):
    """
    Rename Workspace
    """
    fields = ['name']
    template_name = 'workspace/history.html'

    def get_success_url(self):
        return reverse_lazy('history_detail', args=(self.get_object().history,))
