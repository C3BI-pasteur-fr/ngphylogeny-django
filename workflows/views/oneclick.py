from __future__ import absolute_import

from django.utils.decorators import method_decorator

from data.views import UploadView, ImportPastedContentView
from galaxy.decorator import connection_galaxy
from tools.models import Tool
from workflows.models import Workflow
from workflows.views.generic import WorkflowFormView, WorkflowListView


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickListView(WorkflowListView ):
    """
    OneClick workflows class-based view
    """
    template_name = 'workflows/workflows_oneclick_choices.html'
    restrict_toolset = Tool.objects.filter(toolflag__name='oclik')

    def get_queryset(self):
        return super(WorkflowOneClickListView, self).get_queryset()

    def get_context_data(self, **kwargs):
        context = super(WorkflowListView, self).get_context_data(**kwargs)
        context['form'] = UploadView.form_class()
        context['textarea_form'] = ImportPastedContentView.form_class()

        return context

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickFormView(WorkflowFormView):
    """
    Workflow oneclick form with the list of tools and
    Launch oneclick workflow
    """

    def get_workflow(self):
        """create a work copy of workflow in user space, import workflow oneclick shared"""
        wk = self.get_object()
        gi = self.request.galaxy
        try:
            # import published shared workflow
            wf_import = gi.workflows.import_shared_workflow(wk.id_galaxy)
        except:
            # make a working copy workflow
            wk_cp = gi.workflows.export_workflow_dict(wk.id_galaxy)
            wf_import = gi.workflows.import_workflow_dict(wk_cp)

        if wf_import:
            # replace workflow with working copy of original workflow shared
            wk.pk = None
            wk.id_galaxy = wf_import.get('id')
            return wk
        return None


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowStartedView(WorkflowOneClickFormView):
    "redirect to workflows form class"

    def get_object(self, queryset=None, detail=None):
        wk_obj = Workflow.objects.order_by('rank').first()
        wk_obj.json = self.request.galaxy.workflows.show_workflow(workflow_id=wk_obj.id_galaxy)
        return wk_obj


