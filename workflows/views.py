import tempfile

from django.shortcuts import render, redirect, get_object_or_404
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView
from formtools.wizard.views import SessionWizardView
from django.core.urlresolvers import reverse_lazy

from account.decorator import connection_galaxy
from workflows.models import Workflow
from tools.forms import ToolForm
from tools.models import Tool,ToolFlag
from tools.views import tool_exec
from workspace.views import get_or_create_history
from data.views import UploadView
from django.views.generic import ListView

from django.views.generic.detail import SingleObjectMixin

@connection_galaxy
def workflows_oneclick_mode_build(request):

    return render(request, 'workflows/workflows_oneclik.html')


class WorkflowOneClickListView(ListView):
    """
    Generic class-based view
    """
    queryset = Workflow.objects.filter(galaxy_server__galaxyconf__active=True)

@connection_galaxy
def workflow_oneclick_form_view(request, slug):

    gi = request.galaxy
    workflow = get_object_or_404(Workflow, slug=slug)
    wf = gi.workflows.show_workflow(workflow.id_galaxy)

    workflow.json = wf
    i_inputs = wf['inputs'].keys()

    context = {"workflow": workflow, "inputs": i_inputs }
    return render(request, 'workflows/workflows_oneclik_form.html', context)


class WorkflowOneClickView( UploadView ):

    queryset = Workflow.objects.filter(galaxy_server__galaxyconf__active=True)
    object = Workflow
    #context_object_name = 'workflow'
    template_name = 'workflows/workflows_oneclik_form.html'

    def get_context_data(self, **kwargs):

        context = super(WorkflowOneClickView, self).get_context_data(**kwargs)

        gi = self.request.galaxy
        workflow = self.get_object()
        workflow.json = gi.workflows.show_workflow(workflow.id_galaxy)
        print workflow.json
        query=Tool.objects
        start = workflow.json['inputs'].keys()[0]

        graph = {}
        tool_list = []
        for step_id, step in workflow.json.get('steps').items():

            tool_list.append((step_id, query.filter(id_galaxy=step.get('tool_id'))))

            for input, step_output in step.get("input_steps",{}).items():
                graph.setdefault(str(step_output.get('source_step')), []).append(str(step_id))


        def search_order(graph,step):
            """sort steps with naive BFS"""

            step_visited = []
            next_steps = [step]
            while next_steps:
                current_step = next_steps.pop(0)
                if not (current_step in step_visited):
                    step_visited.append(current_step)
                    next_steps.extend(graph.get(current_step,[]))
            return step_visited

        sorted_step = search_order(graph, start)
        sorted_tool_list =[]

        for e, step in enumerate(sorted_step):
            for y, x in tool_list:
                if step == y:
                    if x:
                        sorted_tool_list.append((e,x[0]))


        #TODO use DOT library to generate Graph (SVG)
        #TODO >>>import graphiz

        context["steps"] = sorted_tool_list
        context["workflow"] = workflow
        context["inputs"] = workflow.json['inputs'].keys()

        return context

    def get_object(self, queryset=None):
        return get_object_or_404(Workflow, slug=self.kwargs['slug'])

    def form_valid(self, form):

        upfile = self.upload_file(form)
        file_id = upfile.get('outputs')[0].get('id')

        history_id = get_or_create_history(self.request)
        workflow_id = self.get_object().id_galaxy

        # recupere l'input du worflow oneclick
        gi = self.request.galaxy
        wf = gi.workflows.show_workflow( workflow_id )
        i_input = wf['inputs'].keys()[0]

        # mappe l'input du workflow oneclick avec l'id du fichier soumi par l'utilisateur
        dataset_map = dict()
        dataset_map[i_input] = {'id':file_id, 'src': 'hda'}


        # lancement du workflow avec les parametres
        gi.workflows.run_workflow(workflow_id=workflow_id, history_id=history_id, dataset_map=dataset_map)

        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': history_id}, )

        return super(UploadView, self).form_valid(form)



@connection_galaxy
def launch_galaxy_workflow(request, slug):

    gi = request.galaxy
    history_id = get_or_create_history(request)
    workflow_id = get_object_or_404(Workflow, slug=slug).id_galaxy

    if request.method == 'POST':

            if request.FILES:

                input_submit = []

                for input_file_id in request.FILES:

                    uploaded_file = request.FILES.get(input_file_id)
                    tmp_file = tempfile.NamedTemporaryFile()
                    for chunk in uploaded_file.chunks():
                        tmp_file.write(chunk)
                    tmp_file.flush()

                    # send file to galaxy
                    outputs = gi.tools.upload_file(tmp_file.name, history_id, file_name=uploaded_file.name)
                    input_submit.append(outputs.get('outputs')[0].get('id'))

                # recupere les inputs du worflows
                wf = gi.workflows.show_workflow( workflow_id )
                i_inputs = wf['inputs'].keys()

                # mappe les inputs du workflow avec les id des fichiers soumit par l'utilisateur
                dataset_map = dict()
                for i, r in zip(i_inputs, input_submit):
                    dataset_map[i] = {'id': r, 'src': 'hda'}

                # lancement du worflow avec les parametres
                gi.workflows.run_workflow(workflow_id=workflow_id, history_id=history_id, dataset_map=dataset_map)
                return redirect('history_detail', history_id=history_id)

    return redirect("galaxy_workflow")




def form_class_list(gi, ids_tools):
    """
    Create ToolForm classes on the fly to be used by WizardForm
    :param gi:
    :param ids_tools:
    :return: list af Class form
    """
    tools = Tool.objects.filter(pk__in=ids_tools)
    tools_inputs_details = []
    for tool in tools:
        tool_inputs_details = gi.tools.show_tool(tool_id=tool.id_galaxy, io_details='true')
        tools_inputs_details.append(
            type(str(tool.name) + 'Form',
                 (ToolForm,),
                 {'tool_params': tool_inputs_details.get('inputs'),'tool_id': tool.id_galaxy }
                 )
        )
    return tools_inputs_details


@connection_galaxy
def workflow_form(request):

    tools = Tool.objects.all()
    form_list = form_class_list(request.galaxy, tools)
    selected_tools = request.session.get('selected_tools')

    if selected_tools:
         form_list = form_class_list(request.galaxy, selected_tools)

    return WorkflowWizard.as_view(form_list=form_list)(request)


class WorkflowWizard(SessionWizardView):

    template_name = 'workflows/workflows_form.html'

    def done(self, form_list, **kwargs):
        for form in form_list:
            #launch tools one by one
            output = tool_exec(self.request, form)
            print output

        return redirect(reverse_lazy("history_detail", kwargs={'history_id': output }))



@connection_galaxy
def workflows_advanced_mode_build(request):

    WORKFLOW_ADVANCED_MODE = [{"step": 0, "category": 'algn', "group": ""},
                              {"step": 1, "category": 'clean', "group": ""},
                              {"step": 2, "category": 'tree', "group": ""},
                              {"step": 3, "category": 'visu', "group": ""},
                              ]

    for step in WORKFLOW_ADVANCED_MODE:
        step['group'] = ToolFlag.objects.get(name=step.get('category'))

    if request.method == 'POST':

        tools = []
        for step in WORKFLOW_ADVANCED_MODE:
            tools.append(request.POST.get(step.get('category')))

        request.session['selected_tools'] = tools
        return redirect('workflow_form')

    context = {"workflow": WORKFLOW_ADVANCED_MODE}
    return render(request, 'workflows/workflows_advanced.html', context)


def workflows_alacarte_mode_build(request):

    WORKFLOW_ADVANCED_MODE = [{"step": 0, "category": 'algn', "group": ""},
                              {"step": 1, "category": 'clean', "group": ""},
                              {"step": 2, "category": 'tree', "group": ""},
                              {"step": 3, "category": 'visu', "group": ""},
                              ]

    for step in WORKFLOW_ADVANCED_MODE:
        step['group'] = ToolFlag.objects.get(name=step.get('category'))
    context = {"workflow": WORKFLOW_ADVANCED_MODE}

    return render(request, 'workflows/workflows_alacarte.html', context)
