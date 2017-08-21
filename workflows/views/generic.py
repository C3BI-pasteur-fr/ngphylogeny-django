from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.shortcuts import redirect
from formtools.wizard.views import SessionWizardView

from data.views import UploadView, ImportPastedContentView
from galaxy.decorator import connection_galaxy
from tools.forms import tool_form_factory
from tools.views import tool_exec
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history


class WorkflowFormView(UploadView, ImportPastedContentView):
    """
    Generic Workflow form view, upload one file and run workflow
    """
    form_class = UploadView.form_class
    form2_class = ImportPastedContentView.form_class

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
            context["steps"] = wk_galaxy.steps_tooldict

        context["workflow"] = wk
        return context

    def get_workflow(self):

        return self.get_object()

    def post(self, request, *args, **kwargs):

        # determine which form is being submitted
        # uses the name of the form's submit button
        if request.FILES:

            # get the primary form
            form_class = self.get_form_class()
            form_name = 'form'

        else:

            # get the secondary form
            form_class = self.form2_class
            form_name = 'form2'

        # get the form
        form = self.get_form(form_class)

        # validate
        if form.is_valid():

            return self.form_valid(form)
        else:
            return self.form_invalid(**{form_name: form})


    def form_valid(self, form):

        gi = self.request.galaxy

        # create new history
        history_id = create_history(self.request)

        # upload user file
        submitted_file = form.cleaned_data.get('file')
        if submitted_file:
            u_file = self.upload_file(submitted_file, history_id)

        # or upload user pasted content
        else:
            u_file = self.upload_content(form.cleaned_data['pasted_text'])

        file_id = u_file.get('outputs')[0].get('id')

        workflow = self.get_workflow()

        i_input = workflow.json['inputs'].keys()[0]
        # input file
        dataset_map = dict()
        dataset_map[i_input] = {'id': file_id, 'src': 'hda'}
        print dataset_map

        try:
            # run workflow
            gi.workflows.run_workflow(workflow_id=workflow.id_galaxy,
                                         history_id=history_id,
                                         dataset_map=dataset_map,
                                         # inputs=dataset_map
                                      )  # ,params=wk_galaxy.params)

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
        tools_inputs_details.append(tool_form_factory(tool, galaxy_server=galaxy_server))
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
