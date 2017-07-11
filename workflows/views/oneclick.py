from __future__ import absolute_import

from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from django.views.generic import ListView

from galaxy.decorator import connection_galaxy
from workflows.models import Workflow, WorkflowStepInformation
from workflows.views.generic import WorkflowFormView


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickListView(ListView):
    """
    Generic class-based view
    """
    model = Workflow
    template_name = 'workflows/workflows_oneclick_choices.html'

    def get_queryset(self):

        self.queryset = Workflow.objects.filter(galaxy_server__current=True).select_related()
        gi = self.request.galaxy
        for workflow in self.queryset:
            workflow_json = gi.workflows.show_workflow(workflow_id=workflow.id_galaxy)

            # parse galaxy workflows json information
            workflow.detail = WorkflowStepInformation(workflow_json).sorted_tool_list

        return self.queryset


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickView(WorkflowFormView):
    """
    Workflow oneclick form with the list of tools and
    Launch oneclick workflow
    """

    def get_workflow(self):
        """create a work copy of workflow in user space, import workflow oneclick shared"""
        gi = self.request.galaxy

        wk = self.get_object()
        # import published shared workflow
        wf_import = gi.workflows.import_shared_workflow(wk.id_galaxy)

        # replace workflow with working copy of original workflow shared
        wk.pk = None
        wk.id_galaxy = wf_import.get('id')

        return wk

    def get_object(self,queryset=None):

        wk_obj = get_object_or_404(Workflow, slug=self.kwargs['slug'])
        wk_obj.json = self.request.galaxy.workflows.show_workflow(workflow_id=wk_obj.id_galaxy)

        return wk_obj


