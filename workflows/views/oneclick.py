from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from django.views.generic import ListView

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history


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

            # parse galaxy workflows json informations
            workflow.detail = WorkflowStepInformation(workflow_json).sorted_tool_list

        return self.queryset


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickView(UploadView):
    """
    :return workflow oneclick form with tools used

    """
    object = Workflow
    template_name = 'workflows/workflows_oneclick_choices.html'

    def get_object(self, queryset=None):
        return get_object_or_404(Workflow, slug=self.kwargs['slug'])

    def get_context_data(self, **kwargs):

        context = super(WorkflowOneClickView, self).get_context_data(**kwargs)
        gi = self.request.galaxy
        workflow_obj = self.get_object()
        workflow_json = gi.workflows.show_workflow(workflow_id=workflow_obj.id_galaxy)

        # parse galaxy workflows json informations
        wk_galaxy = WorkflowStepInformation(workflow_json)

        # TODO use DOT library to generate Graph (SVG)
        # TODO >>>import graphiz

        context["steps"] = wk_galaxy.sorted_tool_list
        context["workflow"] = workflow_obj
        context["inputs"] = workflow_json['inputs'].keys()

        return context

    def form_valid(self, form):

        gi = self.request.galaxy

        # import published shared workflow
        shared_workflow_id = self.get_object().id_galaxy

        #wf_import = gi.workflows.import_shared_workflow(shared_workflow_id) pull-request #210
        payload = {'shared_workflow_id':shared_workflow_id}
        wf_import=gi.workflows._post(payload)
        workflow_id = wf_import.get('id')

        if wf_import:

            # upload user file
            upfile = self.upload_file(form)
            file_id = upfile.get('outputs')[0].get('id')
            history_id = create_history(self.request)

            # load workflow
            workflow_json = gi.workflows.show_workflow(wf_import.get('id'))
            wk_galaxy = WorkflowStepInformation(workflow_json)
            i_input = workflow_json['inputs'].keys()[0]

            # mappe l'input du workflow oneclick avec l'id du fichier soumi par l'utilisateur
            dataset_map = dict()
            dataset_map[i_input] = {'id': file_id, 'src': 'hda'}
            print dataset_map
            try:
                # lancement du workflow avec les parametres
                gi.workflows.run_workflow(workflow_id=workflow_id, history_id=history_id, dataset_map=dataset_map)#,params=wk_galaxy.params)
            except Exception, galaxy_exception:
                raise galaxy_exception
            finally:
                # supprime le workflow imported
                print gi.workflows.delete_workflow(workflow_id)

            self.success_url = reverse_lazy("history_detail", kwargs={'history_id': history_id}, )
        else:
            print "Fail import workflow"
        return HttpResponseRedirect(self.get_success_url())