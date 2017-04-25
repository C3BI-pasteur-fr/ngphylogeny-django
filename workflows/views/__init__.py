import tempfile

from django.core.urlresolvers import reverse_lazy
from django.shortcuts import render, redirect, get_object_or_404
from formtools.wizard.views import SessionWizardView

from galaxy.decorator import connection_galaxy
from tools.forms import ToolForm
from tools.models import Tool,ToolFlag
from tools.views import tool_exec
from workflows.models import Workflow
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
