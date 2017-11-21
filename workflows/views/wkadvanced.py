from django.utils.decorators import method_decorator
from django.utils.text import slugify
from django.views.generic import View
from django.views.generic.detail import SingleObjectMixin

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.forms import tool_form_factory
from tools.models import Tool
from workflows.views.generic import WorkflowWizard, WorkflowListView
from workflows.views.viewmixing import WorkflowDuplicateMixin, WorkflowDetailMixin


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
class WorkflowAdvancedListView(WorkflowListView):
    """
        Workflow Advanced ListView
    """
    template_name = 'workflows/workflows_advanced_list.html'
    restrict_toolset = Tool.objects.filter(toolflag__name='wadv')


class WorkflowAdvancedFormView(WorkflowDetailMixin, SingleObjectMixin):
    """
        Generic Workflow Advanced class-based view
    """

    template_name = 'workflows/workflows_advanced_form.html'
    context_object_name = "workflow_list"
    restrict_toolset = Tool.objects.filter(toolflag__name='wadv')

    def get_context_data(self, **kwargs):
        context = super(WorkflowAdvancedFormView, self).get_context_data(**kwargs)
        context['tool_list'] = []
        tools = []

        wk = self.get_object()
        self.fetch_workflow_detail(wk)

        for t in wk.detail:
            tools.append(t[1])
            context['tool_list'].append(slugify(t[1].name))

        context['form_list'] = form_class_list(self.request.galaxy_server, tools)
        context['workflow_list'] = [wk, ]

        return context


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowAdvWizardClass(WorkflowDuplicateMixin, WorkflowAdvancedFormView, WorkflowWizard):
    """
        Wizard Form
    """
    template_name = 'workflows/workflows_advanced_form.html'

    def done(self, form_list, **kwargs):
        render_wizard = super(WorkflowAdvWizardClass, self).done(form_list, **kwargs)

        if self.succes_url:
            # delete the workflow when the workflow has been run
            print (self.clean_copy())

        return render_wizard


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowAdvancedRedirectView(WorkflowAdvancedFormView, View):
    """
        Generic Workflow Advance class-based view
    """

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(object=self.object)
        context['form_list'] = [UploadView.form_class] + context['form_list']

        return WorkflowAdvWizardClass.as_view(form_list=context['form_list'])(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        return self.get(request, *args, **kwargs)
