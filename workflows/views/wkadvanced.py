from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.utils.decorators import method_decorator
from django.utils.text import slugify
from formtools.wizard.views import SessionWizardView

from galaxy.decorator import connection_galaxy
from tools.forms import tool_form_factory
from workflows.models import Workflow, WorkflowStepInformation
from workflows.views.generic import WorkflowDetailView
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

class WorkflowsAdvancedList(WorkflowOneClickListView):

    template_name = 'workflows/workflows_advanced_list.html'

class WorkflowWizard(SessionWizardView):
    """"""
    template_name = 'workflows/workflows_advanced_form.html'


    def done(self, form_list, **kwargs):

        history_id = create_history(self.request)

        params = {}
        for idx, tool_form in enumerate(form_list):
            params[str(idx)] =  inputs()

            for key, value in tool_form.cleaned_data.items():
                # set the Galaxy parameter ( name, value)
                params[str(idx)].set_param(tool_form.fields_ids_mapping.get(key), value)


        workflow = self.get_object(detail=True)

        i_input = workflow.json['inputs'].keys()[0]
        # input file
        dataset_map = dict()
        dataset_map[i_input] = {'id': "file_id", 'src': 'hda'}
        print dataset_map


        # run workflow
        self.request.galaxy.workflows.run_workflow(workflow_id=workflow.id_galaxy,
                                         history_id=history_id,
                                         dataset_map=dataset_map,
                                         # inputs=dataset_map
                                         params=params)


        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': history_id}, )

        return HttpResponseRedirect(self.get_success_url())



@method_decorator(connection_galaxy, name="dispatch")
class WorkflowsAdvancedView(WorkflowDetailView):
    """
    Generic class-based view
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_advanced_form.html'


    def get_context_data(self, **kwargs):
        context = super(WorkflowsAdvancedView, self).get_context_data(**kwargs)

        wk = self.get_object(detail=True)
        gi = self.request.galaxy
        workflow_json = gi.workflows.show_workflow(workflow_id=wk.id_galaxy)

        # parse galaxy workflows json information
        wk.detail = WorkflowStepInformation(workflow_json).sorted_tool_list

        tools = []
        context['tool_list'] = []
        for t in wk.detail:
            tools.append(t[1])
            context['tool_list'].append(slugify(t[1].name))

        context['form_list'] = form_class_list(self.request.galaxy_server, tools)
        context['workflow_list'] = [wk, ]

        return context

class WorkflowAdvWizardClass(WorkflowWizard, WorkflowsAdvancedView ):
    pass


class WorkflowsAdvancedRedirectView(WorkflowsAdvancedView):
    """
    """
    def get(self, request, *args, **kwargs):
        context = self.get_context_data()
        #named_forms = zip( context['tool_list'], context['form_list'])

        #WorkflowAdvWizardClass.form_list= context['form_list']
        #return ClassWizardView.as_view(named_forms, url_name='workflows_advanced_step', done_step_name='home')(request, *args, **kwargs)
        return WorkflowAdvWizardClass.as_view(form_list=context['form_list'])(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        return self.get(request, *args, **kwargs)