from django.views.generic import TemplateView


class HistoryDetailView(TemplateView):
    """
        Display Galaxy like histoiry information
    """
    template_name = 'history.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(HistoryDetailView, self).get_context_data(**kwargs)

        gi = self.request.galaxy

        history_content = gi.histories.show_history(context['history_id'], contents=True)

        context['history_info'] = gi.histories.show_history(context['history_id'])
        context['history_content'] = history_content

        return context
