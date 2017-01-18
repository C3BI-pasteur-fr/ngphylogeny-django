from django.views.generic import TemplateView, RedirectView
from django.utils.decorators import method_decorator
from account.decorator import connection_galaxy
from account.models import GalaxyConf
from models import WorkspaceHistory
from django.shortcuts import get_object_or_404

@connection_galaxy
def create_history(request):
    """
    :param request:
    :return: history_id
    """

    gi = request.galaxy
    # Create a new galaxy history and delete older if the user is not authenticated
    history = gi.histories.create_history()

    galaxy_conf = GalaxyConf.objects.get(active=True)
    if request.user.is_anonymous:

        current_user = galaxy_conf.anonymous_user
    else:
        current_user = request.user

    wsph = WorkspaceHistory(history=history.get("id"),
                            name=history.get('name'),
                            user=current_user,
                            galaxy_server=galaxy_conf.galaxy_server)

    wsph.save()
    request.session["last_history"] = wsph.history

    return wsph.history

@connection_galaxy
def get_history(request):

    gi = request.galaxy
    return request.session.get('last_history')


def get_or_create_history(request):
    """
    :param request:
    :return: history_id
    """
    history_id = get_history(request)
    if not history_id:
        # Create a new galaxy history
        history_id = create_history(request)

    # save the current history
    request.session["last_history"] = history_id

    return history_id


@method_decorator(connection_galaxy, name="dispatch")
class HistoryDetailView(TemplateView):
    """
        Display Galaxy like history information
    """
    template_name = 'workspace/history.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(HistoryDetailView, self).get_context_data(**kwargs)

        # first display history with id history contained in the url
        history_id = context.get('history_id', None) or self.kwargs.get('history_id', None)

        # if no history id try to retrieve or create history
        if not history_id:
            history_id = get_history(self.request)

        if not history_id:
            return context

        gi = self.request.galaxy
        context['history_info'] = gi.histories.show_history(history_id)
        context['history_content'] = gi.histories.show_history(history_id, contents=True)
        return context


class GalaxyErrorView(RedirectView):

    def get_redirect_url(self, *args, **kwargs):
        galaxy_conf = get_object_or_404(GalaxyConf,active=True)

        return "%s/dataset/errors?id=%s" %(galaxy_conf.galaxy_server.url, kwargs.get('id'))