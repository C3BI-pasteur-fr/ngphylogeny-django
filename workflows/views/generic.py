from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect
from django.utils.decorators import method_decorator
from django.utils.functional import cached_property
from django.views.generic import ListView, DetailView
from django.shortcuts import get_object_or_404
from django.views.generic.base import RedirectView
from django.db.models import Q

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from workflows.models import Workflow, WorkflowStepInformation
from workspace.views import create_history
from tools.models import Tool
from blast.models import BlastRun
from workspace.tasks import initializeworkspacejob

import tempfile
import StringIO

from utils import biofile

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowListView(ListView):
    """
    Generic class-based ListView
    """
    model = Workflow
    context_object_name = "workflow_list"
    template_name = 'workflows/workflows_list.html'

    @cached_property
    def workflow_list(self):
        gi = self.request.galaxy
        workflow_queryset = Workflow.objects.filter(
            galaxy_server__current=True).filter(category='base').select_related()

        for workflow in workflow_queryset:
            workflow.fetch_details(gi, self.restricted_toolset)
        return workflow_queryset

    def get_queryset(self):
        return self.workflow_list    

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowFormView(UploadView, DetailView):
    """
    Generic Workflow form view, upload one file and run workflow
    """
    #model = Workflow
    template_name = 'workflows/workflows_form.html'
    restricted_toolset = None

  
    def get_context_data(self, **kwargs):
        gi = self.request.galaxy
        context = super(WorkflowFormView, self).get_context_data(**kwargs)
        # get workflows
        wk = self.get_object()
        wk.fetch_details(gi, self.restricted_toolset)
        if context.get('form') is None:
            context['form'] = UploadView.form()
        if hasattr(wk, 'json'):
            # parse galaxy workflow json information
            context["inputs"] = wk.json['inputs'].keys()
            # add workfow galaxy information
            context["steps"] = WorkflowStepInformation(
                wk.json, tools=self.restricted_toolset,
            ).steps_tooldict
        context["workflow"] = wk
        return context

    def post(self, request, *args, **kwargs):
        # determines which form is being submitted
        # uses the name of the form's submit button
        form = self.get_form()
        # validate
        if form.is_valid():
            return self.form_valid(form)
        else:
            return self.form_invalid(form)

    def form_valid(self, form):
        gi = self.request.galaxy
        wk = self.get_object()
        workflow = wk.duplicate(gi)
        workflow.fetch_details(gi, self.restricted_toolset)
        workflow.save()
        # create new history
    
        # upload user file or pasted content
        submitted_file = form.cleaned_data.get('input_file')
        pasted_text = form.cleaned_data.get('pasted_text')
        blast_run   = form.cleaned_data.get('blast_run')

        # We check that all the tools of the workflow oneclick are allowed to
        # run on this size of data
        if submitted_file:
            submitted_file.seek(0)
            nseq, length, seqaa = biofile.valid_fasta(submitted_file)
            submitted_file.seek(0)
        elif pasted_text:
            nseq, length, seqaa = biofile.valid_fasta(StringIO.StringIO(str(pasted_text)))
        elif blast_run != '--':
            b = BlastRun.objects.get(pk=blast_run)
            nseq, length, seqaa = biofile.valid_fasta(StringIO.StringIO(str(b.to_fasta())))
        else:
            form.add_error(
                'input_file',"No input file given")
            return self.form_invalid(form)

        if nseq == 0:
            form.add_error(
                'input_file', "Input file format is not FASTA or file is empty")
            return self.form_invalid(form)
        elif nseq <= 3:
            form.add_error(
                'input_file',"Input file should contain more than 3 sequences")
            return self.form_invalid(form)
            
        for k,v in workflow.json.get('steps',dict()).items():
             tid = v.get('tool_id',None)
             if tid :
                 t = Tool.objects.get(id_galaxy=tid)
                 if not t.can_run_on_data( nseq, length, -1, seqaa):
                     form.add_error(
                         'input_file',"Input data is too large for the workflow")
                     workflow.delete_from_galaxy(gi)
                     return self.form_invalid(form)

        wksph = create_history(
            self.request,
            name="NGPhylogeny Analyse - " + workflow.name,
            wf_category="OneClick",
            wf_steps=workflow.tooldesc,
        )

        if submitted_file:
            u_file = self.upload_file(submitted_file, wksph.history)
        # or upload user pasted content
        elif pasted_text:
            u_file = self.upload_content(pasted_text)
        elif blast_run :
            u_file = self.upload_content(b.to_fasta(),name="Blast_%s_%s" % (b.query_id,str(blast_run)))
            
        file_id = u_file.get('outputs')[0].get('id')
        i_input = workflow.json['inputs'].keys()[0]
        
        # input file
        dataset_map = dict()
        dataset_map[i_input] = {'id': file_id, 'src': 'hda'}
        try:
            # run workflow
            self.outputs = gi.workflows.run_workflow(
                workflow_id=workflow.id_galaxy,
                history_id=wksph.history,
                dataset_map=dataset_map,
            )
        except Exception as galaxy_exception:
            workflow.delete_from_galaxy(gi)
            raise galaxy_exception
        #finally:
        #    workflow.delete(gi)
        self.success_url = reverse_lazy("history_detail", kwargs={
                                        'history_id': wksph.history}, )
        # Start monitoring
        wksph.monitored = True
        wksph.workflow = workflow
        initializeworkspacejob.delay(wksph.history)
        wksph.save()
        return HttpResponseRedirect(self.get_success_url())


@method_decorator(connection_galaxy, name="dispatch")
class RerunWorkflow(RedirectView):

    permanent = False
    query_string = True
    pattern_name = 'workflow_maker_form'

    # We duplicate the workflow and go to the workflow maker form
    def get_redirect_url(self, *args, **kwargs):
        # Duplicate workflow using ID and go to wmaker url...
        gi = self.request.galaxy
        wk = get_object_or_404(Workflow, id_galaxy=kwargs['id'])
        workflow = wk.duplicate(gi)
        #workflow.fetch_details(gi, self.restricted_toolset)
        workflow.save()
        kwargs['id'] = workflow.id_galaxy
        return super(RerunWorkflow, self).get_redirect_url(*args, **kwargs)
