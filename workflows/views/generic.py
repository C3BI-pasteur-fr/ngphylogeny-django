from bioblend.galaxy.tools.inputs import inputs
from django.core.files.storage import FileSystemStorage
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.utils.decorators import method_decorator
from django.utils.functional import cached_property
from django.views.generic import ListView, DetailView
from formtools.wizard.views import SessionWizardView

from data.views import UploadView, ImportPastedContentView
from galaxy.decorator import connection_galaxy
from viewmixing import WorkflowDetailMixin
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowListView(WorkflowDetailMixin, ListView):
    """
    Generic class-based Listview
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_list.html'
    restrict_toolset = None

    @cached_property
    def workflow_list(self):
        workflow_queryset = Workflow.objects.filter(galaxy_server__current=True).select_related()

        for workflow in workflow_queryset:
            self.fetch_workflow_detail(workflow)
        return workflow_queryset

    def get_queryset(self):
        return self.workflow_list


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowFormView(WorkflowDetailMixin, UploadView, ImportPastedContentView, DetailView):
    """
    Generic Workflow form view, upload one file and run workflow
    """

    form_class = UploadView.form_class
    form2_class = ImportPastedContentView.form_class
    template_name = 'workflows/workflows_form.html'
    restrict_toolset = None

    def get_context_data(self, **kwargs):

        context = super(WorkflowFormView, self).get_context_data(**kwargs)

        # get workflows
        wk = self.get_object()
        self.fetch_workflow_detail(wk)

        context['form'] = UploadView.form_class()
        context['textarea_form'] = ImportPastedContentView.form_class()

        if hasattr(wk, 'json'):
            # parse galaxy workflow json information
            context["inputs"] = wk.json['inputs'].keys()

            # add workfow galaxy information
            context["steps"] = WorkflowStepInformation(wk.json, tools=self.restrict_toolset).steps_tooldict

        context["workflow"] = wk
        return context

    def post(self, request, *args, **kwargs):
        # determine which form is being submitted
        # uses the name of the form's submit button
        if request.FILES:

            # get the primary form (upload file)
            form_class = self.get_form_class()
            form_name = 'form'

        else:

            # get the secondary form (textaera)
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
        workflow = self.get_workflow()

        # create new history
        history_id = create_history(self.request, name="NGPhylogeny Analyse - " + workflow.name)

        # upload user file
        submitted_file = form.cleaned_data.get('input_file')
        if submitted_file:
            u_file = self.upload_file(submitted_file, history_id)

        # or upload user pasted content
        else:
            u_file = self.upload_content(form.cleaned_data['pasted_text'])

        file_id = u_file.get('outputs')[0].get('id')



        i_input = workflow.json['inputs'].keys()[0]
        # input file
        dataset_map = dict()
        dataset_map[i_input] = {'id': file_id, 'src': 'hda'}
        print dataset_map

        try:
            # run workflow
            self.outputs = gi.workflows.run_workflow(workflow_id=workflow.id_galaxy,
                                                     history_id=history_id,
                                                     dataset_map=dataset_map,
                                                     # inputs=dataset_map
                                                     )  # ,params=wk_galaxy.params)

        except Exception, galaxy_exception:
            raise galaxy_exception

        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': history_id}, )

        return HttpResponseRedirect(self.get_success_url())


class WorkflowWizard(SessionWizardView):
    """
    Generic Form Wizard for Advanced and Alacarte mode
    """
    template_name = 'workflows/workflows_wizard_form.html'
    file_storage = FileSystemStorage('/tmp')
    succes_url = ""

    def done(self, form_list, **kwargs):

        workflow = self.get_workflow()
        history_id = create_history(self.request, name="NGPhylogeny Analyse - " + self.workflow.name)

        i_input = workflow.json['inputs'].keys()[0]

        # input file
        dataset_map = {}
        # tool params
        params = {}

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

            # convert dict to json
            params[str(idx)] = params[str(idx)].to_dict()

        try:
            # run workflow
            self.request.galaxy.workflows.run_workflow(workflow_id=workflow.id_galaxy,
                                                       history_id=history_id,
                                                       dataset_map=dataset_map,
                                                       # inputs=dataset_map
                                                       params=params)

            self.succes_url = reverse_lazy("history_detail", kwargs={'history_id': history_id})

        except Exception, galaxy_exception:

            raise galaxy_exception

        return HttpResponseRedirect(self.succes_url)
