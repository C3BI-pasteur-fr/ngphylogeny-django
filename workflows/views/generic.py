from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.shortcuts import redirect
from formtools.wizard.views import SessionWizardView

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.views import tool_exec
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history


class WorkflowFormView(UploadView):
    """
    Generic Workflow form view, upload one file and run workflow
    """

    object = Workflow
    template_name = 'workflows/workflows_form.html'

    def get_context_data(self, **kwargs):

        context = super(WorkflowFormView, self).get_context_data(**kwargs)

        # get workflows
        wk = self.get_object()

        if hasattr(wk, 'json'):
            # parse galaxy workflow json information
            wk_galaxy = WorkflowStepInformation(wk.json)
            context["inputs"] = wk.json['inputs'].keys()
            context["steps"] = wk_galaxy.sorted_tool_list

        context["workflow"] = wk
        return context

    def get_workflow(self):

        return self.get_object()

    def form_valid(self, form):

        gi = self.request.galaxy

        # create new history
        history_id = create_history(self.request)

        # upload user file
        upfile = self.upload_file(form, history_id)
        file_id = upfile.get('outputs')[0].get('id')

        workflow = self.get_workflow()

        i_input = workflow.json['inputs'].keys()[0]

        # input file
        dataset_map = dict()
        dataset_map[i_input] = {'id': file_id, 'src': 'hda'}
        print dataset_map

        try:
            # run workflow
            gi.workflows.invoke_workflow(workflow_id=workflow.id_galaxy,
                                         history_id=history_id,
                                         inputs=dataset_map)  # ,params=wk_galaxy.params)

        except Exception, galaxy_exception:
            raise galaxy_exception

        finally:
            # remove the working copy
            print gi.workflows.delete_workflow(workflow.id_galaxy)

        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': history_id}, )

        return HttpResponseRedirect(self.get_success_url())


class WorkflowWizard(SessionWizardView):
    """"""
    template_name = 'workflows/workflows_form.html'

    def done(self, form_list, **kwargs):
        for form in form_list:
            # launch tools one by one
            output = tool_exec(self.request, form)
            print output

        return redirect(reverse_lazy("history_detail", kwargs={'history_id': output}))


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
    """
    :param request:
    :param slug_workflow:
    :return:
    """
    gi = request.galaxy
    workflow = Workflow.objects.get(slug=slug_workflow)
    workflow_json = gi.workflows.show_workflow(workflow_id=workflow.id_galaxy)
    tools = [t[1] for t in WorkflowStepInformation(workflow_json).sorted_tool_list]

    # parse galaxy workflows json information
    # todo use get_form
    form_list = form_class_list(request.galaxy_server, tools)
    selected_tools = request.session.get('selected_tools')

    if selected_tools:
        form_list = form_class_list(request.galaxy_server, selected_tools)

    ClassWizardView = type(str(slug_workflow) + "Wizard", (WorkflowWizard,), {'form_list': form_list})

    return ClassWizardView.as_view()(request)
