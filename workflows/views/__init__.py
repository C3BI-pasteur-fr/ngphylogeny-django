import json
import tempfile

from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from formtools.wizard.views import SessionWizardView

from galaxy.decorator import connection_galaxy
from tools.models import ToolFlag
from tools.views import tool_exec
from workflows.models import Workflow, WorkflowStepInformation, WorkflowGalaxyFactory
from workspace.views import get_or_create_history


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



class WorkflowWizard(SessionWizardView):

    template_name = 'workflows/workflows_form.html'

    def done(self, form_list, **kwargs):
        for form in form_list:
            #launch tools one by one
            output = tool_exec(self.request, form)
            print output

        return redirect(reverse_lazy("history_detail", kwargs={'history_id': output }))


def form_class_list(galaxy_server, tools):
    """
    Create ToolForm classes on the fly to be used by WizardForm
    :param gi:
    :param tools:
    :return: list af Class form
    """
    tools_inputs_details = []
    for tool in tools:
        tools_inputs_details.append(tool.form_class(galaxy_server=galaxy_server))
    return tools_inputs_details


@connection_galaxy
def workflow_form(request, slug_workflow):

    gi = request.galaxy
    workflow= Workflow.objects.get(slug=slug_workflow)
    workflow_json = gi.workflows.show_workflow(workflow_id=workflow.id_galaxy)
    tools = [ t[1] for t in WorkflowStepInformation(workflow_json).sorted_tool_list]
    # parse galaxy workflows json informations
    form_list = form_class_list(request.galaxy_server, tools )
    selected_tools = request.session.get('selected_tools')

    if selected_tools:
         form_list = form_class_list(request.galaxy_server, selected_tools)

    ClassWizardView = type(str(slug_workflow)+"Wizard", (WorkflowWizard,), {'form_list': form_list} )

    return ClassWizardView.as_view()(request)




WORKFLOW_ADVANCED_MODE = [{"step": 0, "category": 'algn', "group": ""},
                          {"step": 1, "category": 'clean', "group": ""},
                          {"step": 2, "category": 'tree', "group": ""},
                          {"step": 3, "category": 'visu', "group": ""},
                         ]

@connection_galaxy
def workflows_advanced_mode_build(request):

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


@connection_galaxy
def workflows_alacarte_mode_build(request):


    for step in WORKFLOW_ADVANCED_MODE:
        step['group'] = ToolFlag.objects.get(name=step.get('category'))

    if request.method == 'POST':

        return workflow_build_view(request)


    context = {"workflow": WORKFLOW_ADVANCED_MODE}

    return render(request, 'workflows/workflows_alacarte.html', context)



def workflow_build_view(request):

    tools = []
    for step in WORKFLOW_ADVANCED_MODE:
        tools.append(request.POST.get(step.get('category')))

    list_t = tools

    from tools.models import Tool
    dict_tools = Tool.objects.in_bulk(list_t)
    list_tool = [dict_tools.get(int(t)) for t in list_t if t]


    if list_tool:
        gi = request.galaxy
        history_id = get_or_create_history(request)

        wkg = WorkflowGalaxyFactory(list_tool, gi, history_id)
        wkg.name = request.POST.get('name', 'generated workflow')

        return HttpResponse(json.dumps(eval(str(wkg))), content_type='application/json')







