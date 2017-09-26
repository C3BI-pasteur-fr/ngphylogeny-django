from django.utils.decorators import method_decorator
from django.utils.text import slugify
from django.views.generic import DetailView

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.forms import tool_form_factory
from tools.models import Tool
from workflows.models import Workflow, WorkflowStepInformation
from workflows.views.generic import WorkflowMixin, WorkflowWizard
from workflows.views.oneclick import WorkflowOneClickListView


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


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowsAdvancedFormView(DetailView, WorkflowMixin):
    """
    Generic class-based view
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_advanced_form.html'
    restrict_toolset = Tool.objects.filter(toolflag__name='wadv')

    def get_context_data(self, **kwargs):
        context = super(WorkflowsAdvancedFormView, self).get_context_data(**kwargs)

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


class WorkflowAdvWizardClass(WorkflowWizard, WorkflowsAdvancedFormView):
    template_name = 'workflows/workflows_advanced_form.html'


class WorkflowsAdvancedRedirectView(WorkflowsAdvancedFormView):
    """
    """

    def get(self, request, *args, **kwargs):
        context = self.get_context_data()
        context['form_list'] = [UploadView.form_class] + context['form_list']

        return WorkflowAdvWizardClass.as_view(form_list=context['form_list'])(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        return self.get(request, *args, **kwargs)


class WorkflowsAdvancedListView(WorkflowOneClickListView):
    """
        Workflow Advanced ListView
    """
    template_name = 'workflows/workflows_advanced_list.html'
    restrict_toolset = Tool.objects.filter(toolflag__name='wadv')
