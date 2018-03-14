from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.utils.decorators import method_decorator
from django.utils.functional import cached_property
from django.views.generic import ListView, DetailView

from data.views import UploadView, ImportPastedContentView
from galaxy.decorator import connection_galaxy
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history
from .viewmixing import WorkflowDetailMixin
from Bio import SeqIO

import StringIO


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowListView(WorkflowDetailMixin, ListView):
    """
    Generic class-based ListView
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_list.html'

    @cached_property
    def workflow_list(self):
        workflow_queryset = Workflow.objects.filter(
            galaxy_server__current=True).select_related()

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
    restricted_toolset = None

    def get_context_data(self, **kwargs):

        context = super(WorkflowFormView, self).get_context_data(**kwargs)

        # get workflows
        wk = self.get_object()
        self.fetch_workflow_detail(wk)

        if ('form') not in kwargs:
            context['form'] = UploadView.form_class()
            context['textarea_form'] = ImportPastedContentView.form_class()
        else:
            invalid_form = kwargs.get('form')
            if isinstance(invalid_form, ImportPastedContentView.form_class):
                context['textarea_form'] = invalid_form
                context['form'] = UploadView.form_class()
            else:
                context['form'] = invalid_form
                context['textarea_form'] = ImportPastedContentView.form_class()

        if hasattr(wk, 'json'):
            # parse galaxy workflow json information
            context["inputs"] = wk.json['inputs'].keys()

            # add workfow galaxy information
            context["steps"] = WorkflowStepInformation(
                wk.json, tools=self.restricted_toolset).steps_tooldict

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
            form_name = 'textarea_form'
        # get the form
        form = self.get_form(form_class)

        # validate
        if form.is_valid() and self.validate_form_inputs(form):
            return self.form_valid(form)
        else:

            return self.form_invalid(form)

    def validate_form_inputs(self, form):
        valid = True
        # Check uploaded file or pasted content
        submitted_file = form.cleaned_data.get('input_file')
        submitted_text = form.cleaned_data.get('pasted_text')

        if submitted_file:
            nbseq = 0
            print submitted_file
            for r in SeqIO.parse(submitted_file, "fasta"):
                nbseq += 1
            if nbseq == 0:
                form.add_error(
                    'input_file', ValueError("Input file format is not FASTA or file is empty"))
                valid = False
        elif submitted_text:
            nbseq = 0
            for r in SeqIO.parse(StringIO.StringIO(submitted_text), "fasta"):
                nbseq += 1
            if nbseq == 0:
                form.add_error(
                    'pasted_text', ValueError("Input file format is not FASTA or file is empty"))
                valid = False
        return valid

    def form_valid(self, form):
        gi = self.request.galaxy
        workflow = self.get_workflow()
        # create new history
        history_id = create_history(
            self.request, name="NGPhylogeny Analyse - " + workflow.name)
        # upload user file or pasted content
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
        print (dataset_map)
        try:
            # run workflow
            self.outputs = gi.workflows.run_workflow(workflow_id=workflow.id_galaxy,
                                                     history_id=history_id,
                                                     dataset_map=dataset_map,
                                                     # inputs=dataset_map
                                                     )  # ,params=wk_galaxy.params)
        except Exception as galaxy_exception:
            raise galaxy_exception
        self.success_url = reverse_lazy("history_detail", kwargs={
                                        'history_id': history_id}, )

        return HttpResponseRedirect(self.get_success_url())

