from __future__ import absolute_import

from django.utils.decorators import method_decorator

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.models import Tool
from workflows.models import Workflow
from workflows.views.generic import WorkflowFormView, WorkflowListView

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickListView(WorkflowListView, UploadView):
    """
    OneClick workflows class-based view
    """
    template_name = 'workflows/workflows_oneclick_list.html'
    restricted_toolset = Tool.objects.filter(toolflag__name='oclik')

    # def get_context_data(self, **kwargs):
    #     context = super(WorkflowListView, self).get_context_data(**kwargs)
    #     if context.get('form') is None:
    #         context['form'] = UploadView.get_form()
    #     return context


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickFormView(WorkflowFormView):
    """
    Workflow oneclick form with the list of tools and
    Launch oneclick workflow
    """
    model = Workflow
    template_name = 'workflows/workflows_oneclick_form.html'
    object = None
    def form_valid(self, form):
        try:
            render = super(WorkflowOneClickFormView, self).form_valid(form)
        except Exception as galaxy_exception:
            raise galaxy_exception
        return render
