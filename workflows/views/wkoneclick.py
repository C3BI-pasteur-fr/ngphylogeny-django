from __future__ import absolute_import

from django.utils.decorators import method_decorator

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.models import Tool
from workflows.models import Workflow
from workflows.views.generic import WorkflowFormView, WorkflowListView

from workspace.tasks import monitorworkspace

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickListView(WorkflowListView):
    """
    OneClick workflows class-based view
    """
    template_name = 'workflows/workflows_oneclick_list.html'
    restricted_toolset = Tool.objects.filter(toolflag__name='oclik')

    def get_queryset(self):
        return super(WorkflowOneClickListView, self).get_queryset()

    def get_context_data(self, **kwargs):
        context = super(WorkflowListView, self).get_context_data(**kwargs)
        if context.get('form') is None:
            context['form'] = UploadView.form_class()

        return context


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickFormView(WorkflowFormView):
    """
    Workflow oneclick form with the list of tools and
    Launch oneclick workflow
    """
    model = Workflow
    template_name = 'workflows/workflows_oneclick_form.html'

    def form_valid(self, form):
        try:
            render = super(WorkflowOneClickFormView, self).form_valid(form)

        except Exception as galaxy_exception:
            raise galaxy_exception

        finally:
            self.clean_copy()
        return render

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowStartedView(WorkflowOneClickFormView):
    """
        Return form View to run the first One Click Workflow
    """

    def get_object(self, queryset=None, detail=None):
        wk_obj = self.model.objects.order_by('rank').first()
        wk_obj.json = self.request.galaxy.workflows.show_workflow(workflow_id=wk_obj.id_galaxy)
        return wk_obj
