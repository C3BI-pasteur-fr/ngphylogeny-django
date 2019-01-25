from django.utils.decorators import method_decorator
from django.utils.text import slugify
from django.shortcuts import render
from django.views.generic import View
from django.views.generic.detail import SingleObjectMixin
from django.core.files.uploadedfile import InMemoryUploadedFile, TemporaryUploadedFile
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect

import tempfile

from galaxy.decorator import connection_galaxy
from tools.models import Tool
from tools.models import ToolFieldWhiteList
from tools.forms import ToolForm
from workspace.views import create_history, delete_history
from workflows.views.generic import WorkflowListView
from workflows.exceptions import WorkflowInvalidFormError
from workflows.models import Workflow
from workflows.exceptions import WorkflowInputFileFormatError
from blast.models import BlastRun
from workspace.tasks import initializeworkspacejob

from bioblend.galaxy.tools.inputs import inputs
from Bio import SeqIO

WORKFLOW_ADV_FLAG = "wadv"


def make_form(tool, request=None):
    """
    Instantiate one form, based on the request and the tool
    :param tool:
    :param request: the current post request
    :return: an instanciated and initialized form
    """

    tool_inputs_details = tool.fetch_tool_json()
    tool_field_white_list, created = ToolFieldWhiteList.objects.get_or_create(
        tool=tool, context="w")
    formname = str(slugify(tool.name).title().replace('-', '')) + 'Form'
    prefix = formname.lower()

    # If not post, it means that we want to create a form with
    # default values
    formdata = None
    if request is not None and request.POST:
        formdata = request.POST

    return ToolForm(data=formdata,
                    prefix=prefix,
                    tool_params=tool_inputs_details.get('inputs'),
                    tool_id=tool.id_galaxy,
                    tool_name=tool.name,
                    whitelist=tool_field_white_list.saved_params,
                    fields_ids_mapping={},
                    n=0)


def form_list(tools, request=None):
    """
    Instantiate a list of forms
    givent the tools and the request
    """
    tool_forms = []
    for tool in tools:
        tool_forms.append(make_form(tool, request))
    return tool_forms


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowAdvancedListView(WorkflowListView):
    """
        Workflow Advanced ListView
    """
    template_name = 'workflows/workflows_advanced_list.html'
    restricted_toolset = Tool.objects.filter(toolflag__name=WORKFLOW_ADV_FLAG)
     

@method_decorator(connection_galaxy, name="dispatch")
class WorkflowAdvancedFormView(SingleObjectMixin,
                               View):
    template_name = 'workflows/workflows_adv_singlepage_form.html'
    model = Workflow
    object = None
    context_object_name = "workflow"
    restricted_toolset = Tool.objects.filter(toolflag__name=WORKFLOW_ADV_FLAG)

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(object=self.object)
        return render(request, self.template_name, context)

    def get_context_data(self, **kwargs):
        gi = self.request.galaxy
        if not self.object:
            self.object = self.get_object()
        # Workflow
        self.object.fetch_details(gi, self.restricted_toolset)
        context = super(
            WorkflowAdvancedFormView,
            self).get_context_data(**kwargs)
        context['workflow_list'] = [self.object, ]
        context['tool_list'] = []
        tools = []

        if self.request.session.get('blastruns'):
            blastruns = []
            for b in BlastRun.objects.filter(pk__in=self.request.session['blastruns'], deleted=False, status='F').order_by('-date').all():
                x = (str(b.id),str(b.query_id))
                blastruns.append(x)
            context['blastruns'] = blastruns
        
        for t in self.object.detail:
            tools.append(t[1])
            context['tool_list'].append(slugify(t[1].name))

        context['form_list'] = form_list(tools, self.request)

        return context

    # Copies and checks format of input
    # file to upload on galaxy server
    def process_file_to_upload(self, file_to_upload):
        uploadfile_name = ""
        if isinstance(file_to_upload, InMemoryUploadedFile) or isinstance(file_to_upload, TemporaryUploadedFile):
            tmp_file = tempfile.NamedTemporaryFile()
            for chunk in file_to_upload.chunks():
                tmp_file.write(chunk)
            tmp_file.flush()
            uploadfile_name = str(file_to_upload.name)
        else:
            uploadfile_name = "Upload File"
            tmp_file = tempfile.NamedTemporaryFile()
            tmp_file.write(file_to_upload)
            tmp_file.flush()

        # Check that input file is Fasta and is not empty
        nseq, length = valid_fasta(tmp_file.name)
        if nseq < 4 :
            raise WorkflowInputFileFormatError(
                "Input data is malformed or contain less than 4 sequences"
            )
        return tmp_file, uploadfile_name, nseq, length

    def check_form_validity(self, request, context):
        for form in context['form_list']:
            if not form.is_valid():
                return False
        return True

    def analyze_forms(self, request, context, workflow, params, gi, wksph, nseq, length):
        steps = workflow.json['steps']
        step_id = u'0'
        for tool_form in context['form_list']:
            if not tool_form.is_valid():
                raise WorkflowInvalidFormError(
                    "One form is invalid %s " % (tool_form.prefix))
            tid = getattr(tool_form, 'tool_id', 'null')
            t=Tool.objects.get(id_galaxy=tid)
            tool_inputs = inputs()
            # mapping between form id to galaxy params names
            fields = tool_form.fields_ids_mapping
            inputs_data = set(tool_form.input_file_ids)
            # set the Galaxy parameter (name, value)
            nboot = 0
            boot = False
            for key, value in tool_form.cleaned_data.items():
                if key not in inputs_data:
                    if fields.get(key,"") == 'bootstrap|replicates':
                        nboot = value
                    if fields.get(key,"") == 'bootstrap|do_bootstrap':
                        boot = True
                    tool_inputs.set_param(fields.get(key), value)
            if not boot:
                nboot = 0
            if not t.can_run_on_data(nseq, length, nboot):
                raise WorkflowInputFileFormatError(
                    "Input data is too large for the workflow"
                )
            for inputfile in inputs_data:
                uploaded_file = ""
                if request.FILES:
                    uploaded_file = request.FILES.get(inputfile, '')
                if uploaded_file:
                    tmp_file = tempfile.NamedTemporaryFile()
                    for chunk in uploaded_file.chunks():
                        tmp_file.write(chunk)
                    tmp_file.flush()
                    # send file to galaxy
                    outputs = gi.tools.upload_file(
                        path=tmp_file.name,
                        file_name=uploaded_file.name,
                        history_id=wksph.history)
                    file_id = outputs.get('outputs')[0].get('id')
                    tool_inputs.set_dataset_param(
                        fields.get(inputfile), file_id)
                else:
                    # else paste content
                    content = tool_form.cleaned_data.get(inputfile)
                    if content:
                        tmp_file = tempfile.NamedTemporaryFile()
                        tmp_file.write(content)
                        tmp_file.flush()
                        # send file to galaxy
                        input_fieldname = tool_form.fields_ids_mapping.get(
                            inputfile)
                        outputs = gi.tools.upload_file(
                            path=tmp_file.name,
                            file_name=input_fieldname + " pasted_sequence",
                            history_id=wksph.history)
                        file_id = outputs.get('outputs')[0].get('id')
                        tool_inputs.set_dataset_param(fields.get(inputfile),
                                                      file_id)
            # workflow step
            # get from which step the tools are used
            for i, step in steps.items():
                if (getattr(tool_form, 'tool_id', 'null') ==
                        step.get('tool_id')):
                    step_id = i
                    break
            # convert inputs to dict
            params[step_id] = tool_inputs.to_dict()
            if not params[step_id]:
                del params[step_id]

    def post(self, request, *args, **kwargs):
        gi = request.galaxy

        # Get a copy of the workflow with full details
        wf = self.get_object()
        # If the workflow is not built by the user:
        # (ngphylogeny basic workflows) then we make a copy
        # otherwise we take the wf as is
        # (it will be deleted at the end)
        if wf.category == 'automaker':
            workflow = wf
        else:
            workflow = wf.duplicate(gi)

        workflow.fetch_details(gi, self.restricted_toolset)

        context = self.get_context_data(object=self.object)

        # input file
        dataset_map = {}
        # tool params
        params = {}
        # Workflow inputs
        i_input = workflow.json['inputs'].keys()[0]

        # Handle workflow main input file
        # before creating the workspace etc.
        uploaded_file = request.FILES.get("file") or request.POST.get("file")
        # We check that a file has been given
        if not uploaded_file:
            context = self.get_context_data(object=self.object)
            context['fileerror'] = "No input file given"
            workflow.delete(gi)
            return render(request, self.template_name, context)
        # Then we check input file format
        nseq=0
        length=0
        try:
            tmp_file, uploadfile_name, nseq, length = self.process_file_to_upload(
                uploaded_file)
        except WorkflowInputFileFormatError as e:
            context = self.get_context_data(object=self.object)
            context['fileerror'] = str(e)
            workflow.delete(gi)
            return render(request, self.template_name, context)

        # We check form validity
        if not self.check_form_validity(request, context):
            workflow.delete(gi)
            return self.get(request, *args, **kwargs)

        # We create an history (local and on galaxy)
        wksph = create_history(
            self.request, name="NGPhylogeny Analyse - " + workflow.name)
        # we send the file to galaxy
        output = gi.tools.upload_file(path=tmp_file.name,
                                      file_name=uploadfile_name,
                                      history_id=wksph.history)
        galaxy_file = output.get('outputs')[0].get('id')
        dataset_map[i_input] = {'id': galaxy_file, 'src': 'hda'}

        # We analyze submited forms and upload files to
        # galaxy
        try:
            self.analyze_forms(request, context, workflow, params, gi, wksph, nseq, length)
        except WorkflowInvalidFormError as e:
            # if one form is not valid
            workflow.delete(gi)
            delete_history(wksph.history)
            return self.get(request, *args, **kwargs)
        except WorkflowInputFileFormatError as e:
            context = self.get_context_data(object=self.object)
            context['fileerror'] = str(e)
            workflow.delete(gi)
            return render(request, self.template_name, context)

        # We run the galaxy workflow
        try:
            output = self.request.galaxy.workflows.invoke_workflow(
                workflow_id=workflow.id_galaxy,
                history_id=wksph.history,
                inputs=dataset_map,
                params=params,
                allow_tool_state_corrections=True,
            )

            self.succes_url = reverse_lazy("history_detail", kwargs={
                                           'history_id': wksph.history})
            # Start monitoring (for sending emails)
            initializeworkspacejob.delay(wksph.history)
            wksph.monitored = True
            wksph.save()

            return HttpResponseRedirect(self.succes_url)

        except Exception:
            delete_history(wksph.history)
            raise
        finally:
            # delete the workflow copy of oneclick workflow when
            # the workflow has been run
            workflow.delete(gi)


def valid_fasta(fasta_file):
    # Check uploaded file or pasted content
    nbseq = 0
    length = 0
    for r in SeqIO.parse(fasta_file, "fasta"):
        tlen = len(r.seq)
        length = tlen if tlen > length else length
        nbseq += 1
    return (nbseq, length)
