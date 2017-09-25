from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.utils.decorators import method_decorator
from django.views.generic import DetailView, ListView
from django.views.generic.edit import SingleObjectMixin

from data.views import UploadView, ImportPastedContentView
from galaxy.decorator import connection_galaxy
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history


class WorkflowMixin(SingleObjectMixin):
    """
    """
    model = Workflow
    object = Workflow
    template_name = 'workflows/workflows_form.html'

    def get_object(self, queryset=None, detail=True):
        """Get Workflow object"""
        wk_obj = super(WorkflowMixin, self).get_object()

        # add more information load from Galaxy
        if detail and wk_obj:
            wk_obj.json = self.request.galaxy.workflows.show_workflow(workflow_id=wk_obj.id_galaxy)

        return wk_obj

    def get_workflow(self):
        return self.get_object(detail=True)


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowDetailView(DetailView, WorkflowMixin):
    """
        Workflow Detail View
    """


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowFormView(UploadView, ImportPastedContentView, WorkflowMixin):
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
        wk = self.get_workflow()

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


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowListView(ListView):
    """
    Generic class-based view
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_list.html'
    restrict_toolset = None

    def get_queryset(self):
        self.queryset = Workflow.objects.filter(galaxy_server__current=True).select_related()

        gi = self.request.galaxy
        for workflow in self.queryset:
            workflow_json = gi.workflows.show_workflow(workflow_id=workflow.id_galaxy)

            # add galaxy workflows json information
            workflow.detail = WorkflowStepInformation(workflow_json, tools=self.restrict_toolset).sorted_tool_list

        return self.queryset
