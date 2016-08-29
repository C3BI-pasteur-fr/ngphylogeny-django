from django.views.generic import TemplateView
from django.utils.decorators import method_decorator
from account.decorator import connection_galaxy


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

        if not history_id:
            try:
                history_id = gi.histories.get_most_recently_used_history().get('id')
            except:
                # Create a new galaxy history and delete older if the user is not authenticated
                history_id = gi.histories.create_history().get('id')

        context['history_info'] = gi.histories.show_history(history_id)
        context['history_content'] = gi.histories.show_history(history_id, contents=True)
        print context['history_content']
        return context
