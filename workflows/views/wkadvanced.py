from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.utils.decorators import method_decorator
from django.utils.text import slugify
from formtools.wizard.views import SessionWizardView

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.forms import tool_form_factory
from tools.models import Tool
from workflows.models import Workflow, WorkflowStepInformation
from workflows.views.generic import WorkflowDetailView, WorkflowMixin
from workflows.views.oneclick import WorkflowOneClickListView
from workspace.views import create_history


def form_class_list(galaxy_server, tools):
    """
    Create ToolForm classes on the fly to be used by WizardForm
    :param gi:
    :param tools:
    :return: list af Class form
    """
    tool_form_class = []
    for tool in tools:
        tool_form_class.append(tool_form_factory(tool, galaxy_server=galaxy_server))
    return tool_form_class


from django.core.files.storage import FileSystemStorage


class WorkflowWizard(SessionWizardView, WorkflowMixin):
    """

    """
    template_name = 'workflows/workflows_advanced_form.html'
    file_storage = FileSystemStorage('/tmp')

    def done(self, form_list, **kwargs):

        history_id = create_history(self.request)
        params = {}

        workflow = self.get_object(detail=True)
        i_input = workflow.json['inputs'].keys()[0]
        # input file
        dataset_map = dict()

        for idx, tool_form in enumerate(form_list):
            params[str(idx)] = inputs()

            for key, form in tool_form.cleaned_data.items():
                if "file" in key:

                    output = self.request.galaxy.tools.upload_file(path=str(form.file), file_name=str(form.name),
                                                                   history_id=history_id)
                    galaxy_file = output.get('outputs')[0]
                    self.file_storage.delete(form)
                    dataset_map[i_input] = {'id': galaxy_file.get('id'), 'src': 'hda'}

                else:
                    # set the Galaxy parameter ( name, value)
                    params[str(idx)].set_param(tool_form.fields_ids_mapping.get(key), form)
            params[str(idx)] = params[str(idx)].to_dict()

        print dataset_map

        # run workflow
        self.request.galaxy.workflows.run_workflow(workflow_id=workflow.id_galaxy,
                                                   history_id=history_id,
                                                   dataset_map=dataset_map,
                                                   # inputs=dataset_map
                                                   params=params)

        return HttpResponseRedirect(reverse_lazy("history_detail", kwargs={'history_id': history_id}))


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowsAdvancedView(WorkflowDetailView):
    """
    Generic class-based view
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_advanced_form.html'
    restrict_toolset = Tool.objects.filter(toolflag__name='wadv')

    def get_context_data(self, **kwargs):
        context = super(WorkflowsAdvancedView, self).get_context_data(**kwargs)

        wk = self.get_object(detail=True)
        gi = self.request.galaxy
        workflow_json = gi.workflows.show_workflow(workflow_id=wk.id_galaxy)

        # parse galaxy workflows json information
        wk.detail = WorkflowStepInformation(workflow_json, tools=self.restrict_toolset).sorted_tool_list

        tools = []
        context['tool_list'] = []
        for t in wk.detail:
            tools.append(t[1])
            context['tool_list'].append(slugify(t[1].name))

        context['form_list'] = form_class_list(self.request.galaxy_server, tools)
        context['workflow_list'] = [wk, ]

        return context


class WorkflowAdvWizardClass(WorkflowWizard, WorkflowsAdvancedView):
    pass


class WorkflowsAdvancedRedirectView(WorkflowsAdvancedView):
    """
    """

    def get(self, request, *args, **kwargs):
        context = self.get_context_data()
        context['form_list'] = [UploadView.form_class] + context['form_list']

        return WorkflowAdvWizardClass.as_view(form_list=context['form_list'])(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        return self.get(request, *args, **kwargs)


class WorkflowsAdvancedList(WorkflowOneClickListView):
    """
        Workflow Advanced ListView
    """
    template_name = 'workflows/workflows_advanced_list.html'
    restrict_toolset = Tool.objects.filter(toolflag__name='wadv')
