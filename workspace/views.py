from django.views.generic import TemplateView
from django.utils.decorators import method_decorator
from account.decorator import connection_galaxy
from account.models import GalaxyConf
from models import WorkspaceHistory


@connection_galaxy
def create_history(request):
    gi = request.galaxy
    # Create a new galaxy history and delete older if the user is not authenticated
    history = gi.histories.create_history()

    galaxy_conf = GalaxyConf.objects.get(active=True)
    if request.user.is_anonymous:
        current_user = None
    else:
        current_user = request.user

    wsph = WorkspaceHistory(history=history.get("id"),
                            name=history.get('name'),
                            user=current_user,
                            galaxy_server=galaxy_conf.galaxy_server)

    wsph.save()
    return wsph.history


@method_decorator(connection_galaxy, name="dispatch")
class HistoryDetailView(TemplateView):
    """
        Display Galaxy like history information
    """
    template_name = 'workspace/history.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(HistoryDetailView, self).get_context_data(**kwargs)

        gi = self.request.galaxy
        history_id = context.get('history_id', None)

        #if no history id try to retrieve last history
        if not history_id:
            history_id = self.request.session.get('last_history')

        if not history_id:
            # Create a new galaxy history
            history_id = create_history(self.request)

        #save the current history
        self.request.session["last_history"] = history_id

        context['history_info'] = gi.histories.show_history(history_id)
        context['history_content'] = gi.histories.show_history(history_id, contents=True)
        print context['history_content']
        return context

