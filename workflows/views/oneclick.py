from django.core.urlresolvers import reverse_lazy
from django.shortcuts import get_object_or_404,render
from django.utils.decorators import method_decorator
from django.views.generic import ListView

from account.decorator import connection_galaxy
from data.views import UploadView
from tools.models import Tool
from workflows.models import Workflow
from workspace.views import get_or_create_history


class WorkflowStepInformation(object):
    """
    Parse Galaxy Workflow json:
        :graph: {step_input: [step_output, step_output], ...}
        :params: {step_id: param_dict}
        :tool_list: list of tuple, [(step_id, queryset.tool)..]
        :sorted_tool_list: list of tuple, [(step_id, queryset.tool)..]
    """

    @staticmethod
    def __graph_search(graph, node):
        """Breadth-first search (BFS) is an algorithm
           :return: list of node
        """
        step_visited = []
        next_steps = [node]
        while next_steps:
            current_step = next_steps.pop(0)
            if not (current_step in step_visited):
                step_visited.append(current_step)
                next_steps.extend(graph.get(current_step, []))
        return step_visited


    @staticmethod
    def __sort_tools( sorted_step, unsorted_tool_list):
        """
        :param sorted_step: list of steps
        :param unsorted_tool_list: list of tuple (step, tool)
        :return: list of tuple sorted (step, tools)
        """

        sorted_tool_list = []
        for e, step in enumerate(sorted_step):
            for y, x in unsorted_tool_list:
                if step == y:
                    if x:
                        sorted_tool_list.append((e, x[0]))
        return sorted_tool_list

    def __init__(self, workflow_json):
        """
        :param workflow_json: Galaxy workflow json
        """
        self.workflow_json = workflow_json
        self.graph = {}
        self.params = {}
        self.tool_list = []
        self.sorted_tool_list =[]

        # get known tools
        query = Tool.objects.filter()

        # parse galaxy workflow information
        for step_id, step in self.workflow_json.get('steps').items():
            self.params[step_id] = step.get('tool_inputs')
            self.tool_list.append([step_id, query.filter(id_galaxy=step.get('tool_id'))])

            for input, step_output in step.get("input_steps", {}).items():
                self.graph.setdefault(str(step_output.get('source_step')), []).append(str(step_id))


        first_step = self.workflow_json['inputs'].keys()[0]
        sorted_step = self.__graph_search( self.graph,first_step )
        self.sorted_tool_list = self.__sort_tools(sorted_step, self.tool_list)


@connection_galaxy
def workflow_oneclick_form_view(request, slug):

    gi = request.galaxy
    workflow = get_object_or_404(Workflow, slug=slug)
    wf = gi.workflows.show_workflow(workflow.id_galaxy)

    workflow.json = wf
    i_inputs = wf['inputs'].keys()

    context = {"workflow": workflow, "inputs": i_inputs }
    return render(request, 'workflows/workflows_oneclik_form.html', context)


class WorkflowOneClickListView(ListView):
    """
    Generic class-based view
    """
    queryset = Workflow.objects.filter(galaxy_server__galaxyconf__active=True)


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowOneClickView(UploadView):
    """
    :return workflow oneclick form with tools used

    """
    object = Workflow
    template_name = 'workflows/workflows_oneclik_form.html'

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
        workflow_id = self.get_object().id_galaxy
        # wf_import = gi.workflows.import_shared_workflow(workflow_id)  #return nothing
        wf_export_dict = gi.workflows.export_workflow_json(workflow_id)
        wf_import = gi.workflows.import_workflow_json(wf_export_dict)

        if wf_import:

            # upload user file
            upfile = self.upload_file(form)
            file_id = upfile.get('outputs')[0].get('id')
            history_id = get_or_create_history(self.request)

            # load workflow
            workflow_json = gi.workflows.show_workflow(wf_import.get('id'))
            wk_galaxy = WorkflowStepInformation(workflow_json)
            i_input = workflow_json['inputs'].keys()[0]

            # mappe l'input du workflow oneclick avec l'id du fichier soumi par l'utilisateur
            dataset_map = dict()
            dataset_map[i_input] = {'id': file_id, 'src': 'hda'}

            # lancement du workflow avec les parametres
            gi.workflows.run_workflow(workflow_id=workflow_id, history_id=history_id,params=wk_galaxy.params, dataset_map=dataset_map)

            # supprime le workflow importe
            gi.workflows.delete_workflow(wf_import.get('id'))

            self.success_url = reverse_lazy("history_detail", kwargs={'history_id': history_id}, )
        else:
            print "Fail import workflow"
        return super(UploadView, self).form_valid(form)