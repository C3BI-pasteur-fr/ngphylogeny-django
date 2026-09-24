from django.utils.decorators import method_decorator
from django.utils.text import slugify
from django.shortcuts import render
from django.views.generic import View
from django.views.generic.detail import SingleObjectMixin
from django.core.files.uploadedfile import InMemoryUploadedFile, TemporaryUploadedFile
from django.urls import reverse_lazy
from django.http import HttpResponseRedirect

import tempfile

from galaxy.decorator import connection_galaxy
from tools.models import Tool
from tools.models import ToolFieldWhiteList
from tools.forms import ToolForm
from workspace.views import create_history, delete_history
from workflows.views.generic import (
    GALAXY_UNREACHABLE_EXCEPTIONS, WorkflowListView,
    galaxy_unavailable_response)
from workflows.exceptions import WorkflowInvalidFormError
from workflows.models import Workflow
from workflows.exceptions import WorkflowInputFileFormatError
from blast.models import BlastRun
from workspace.tasks import initializeworkspacejob

from bioblend.galaxy.tools.inputs import inputs

from utils import biofile

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
        try:
            context = self.get_context_data(object=self.object)
        except GALAXY_UNREACHABLE_EXCEPTIONS:
            return galaxy_unavailable_response(request)
        return render(request, self.template_name, context)


    def compatible_inputs(self, extensions, files):
        """
        Returns a list of files [{ext:,history:,name:}] from self.session_files
        that are compatible with the given list of extenstions
        """
        outlist = []
        if files:
            for key, sf in files.items():
                if sf.get('ext') in extensions:
                    outlist.append(sf)
        return outlist
    
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

        if self.request.session.get('files'):
            # We get the first tool (after input_tool) of the workflow
            first_tool_id = self.object.json['steps'].get('1').get('tool_id')
            # Then we get the acceptable extentions of this tool
            extensions = gi.tools.show_tool(first_tool_id,io_details=True).get('inputs')[0].get('extensions')
            compatibleinputs = self.compatible_inputs(extensions,self.request.session.get('files'))
            context['compatibleinputs'] = compatibleinputs
            
        #gi.attrfield.get('extensions')
        #compatible_inputs = compatibleInputs()
        
            
        for t in self.object.detail:
            tools.append(t[1])
            context['tool_list'].append(slugify(t[1].name))

        context['form_list'] = form_list(tools, self.request)

        return context

    # Copies and checks format of input
    # file to upload on galaxy server
    def process_file_to_upload(self, file_to_upload, uploadfile_name):
        if isinstance(file_to_upload, InMemoryUploadedFile) or isinstance(file_to_upload, TemporaryUploadedFile):
            tmp_file = tempfile.NamedTemporaryFile()
            for chunk in file_to_upload.chunks():
                tmp_file.write(chunk)
            tmp_file.flush()
        else:
            # Reached for pasted text (request.POST.get("file") is a plain
            # str, unlike the uploaded-file case above) and BlastRun.to_fasta()
            # results - NamedTemporaryFile() defaults to binary mode, and
            # tmp_file.write(file_to_upload) here used to crash outright
            # under Python 3 ("a bytes-like object is required, not 'str'")
            # for any str input - only caught by actually pasting text
            # through the live A La Carte / advanced form, not by any
            # existing test.
            tmp_file = tempfile.NamedTemporaryFile()
            data = file_to_upload
            if isinstance(data, str):
                data = data.encode('utf-8')
            tmp_file.write(data)
            tmp_file.flush()

        # Rewrite sequence ids to something every downstream Galaxy tool
        # in the pipeline (MAFFT, PhyML/PhyML-SMS, newick_utilities'
        # nw_display, ...) will tokenize identically - see
        # biofile.sanitize_fasta_content's own docstring for the real
        # production bug (a non-breaking space survived alignment/tree
        # building untouched, then broke the Newick Display step) this
        # is here to prevent from recurring.
        tmp_file.seek(0)
        sanitized = biofile.sanitize_fasta_content(tmp_file.read())
        tmp_file.seek(0)
        tmp_file.truncate()
        tmp_file.write(sanitized)
        tmp_file.flush()

        # Check that input file is Fasta and is not empty
        # open() in binary mode, not text mode: valid_fasta() branches
        # on isinstance(raw, bytes) to decode with errors='replace' -
        # but a plain open(tmp_file.name) (text mode, the default)
        # already tries to decode as UTF-8 *inside* .read() itself,
        # before valid_fasta() ever gets a chance to handle it, and
        # raises UnicodeDecodeError uncaught for any non-UTF-8 upload
        # (e.g. a real UTF-16 fasta file saved from Windows Notepad/
        # Excel - 0xFF as the very first byte is a UTF-16LE BOM).
        nseq, length, seqaa = biofile.valid_fasta(open(tmp_file.name, 'rb'))
        if nseq < 4 :
            raise WorkflowInputFileFormatError(
                "Input data is malformed or contain less than 4 sequences"
            )
        return tmp_file, uploadfile_name, nseq, length, seqaa

    def check_form_validity(self, request, context):
        for form in context['form_list']:
            if not form.is_valid():
                return False
        return True

    def analyze_forms(self, request, context, workflow, params, gi, wksph, nseq, length, seqaa):
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
                    if fields.get(key,"") == 'bootstrap|do_bootstrap' and value == 'true':
                        boot = True
                    tool_inputs.set_param(fields.get(key), value)
            if not boot:
                nboot = 0
            if not t.can_run_on_data(nseq, length, nboot, seqaa):
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
                        file_name=uploaded_file.name.encode('ascii','ignore').decode('ascii'),
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
        if wf.category == 'automaker' or wf.category == 'duplicated':
            workflow = wf
        else:
            workflow = wf.duplicate(gi)

        workflow.fetch_details(gi, self.restricted_toolset)
        workflow.save()
        
        context = self.get_context_data(object=self.object)

        # input file
        dataset_map = {}
        # tool params
        params = {}
        # Workflow inputs
        i_input = list(workflow.json['inputs'].keys())[0]

        # Handle workflow main input file
        # before creating the workspace etc.
        uploaded_file = request.FILES.get("file") or request.POST.get("file")
        blastrun = request.POST.get("blastrun")
        galaxy_file = request.POST.get("galaxyfile")
        # Then we check input file format
        nseq=0
        length=0
        seqaa=False
        if galaxy_file == "--":
            # We check that a file has been given
            if blastrun != "--":
                b = BlastRun.objects.get(pk=blastrun)
                uploaded_file = b.to_fasta()
                upload_filename= "Blast_%s_%s" % (b.query_id,str(blastrun))
            elif not uploaded_file:
                context = self.get_context_data(object=self.object)
                context['fileerror'] = "No input file given"
                workflow.delete_from_galaxy(gi)
                return render(request, self.template_name, context)
            elif isinstance(uploaded_file, InMemoryUploadedFile) or isinstance(uploaded_file, TemporaryUploadedFile):
                upload_filename = uploaded_file.name.encode('ascii','ignore').decode('ascii')
            else:
                upload_filename = "uploaded_content"
                
            try:
                tmp_file, uploadfile_name, nseq, length, seqaa = self.process_file_to_upload(
                    uploaded_file,
                    upload_filename,
                )
            except WorkflowInputFileFormatError as e:
                context = self.get_context_data(object=self.object)
                context['fileerror'] = str(e)
                workflow.delete_from_galaxy(gi)
                return render(request, self.template_name, context)

        # We check form validity
        if not self.check_form_validity(request, context):
            workflow.delete_from_galaxy(gi)
            return self.get(request, *args, **kwargs)

        # We create an history (local and on galaxy)
        cat = workflow.category
        if cat == "":
            cat = "Advanced"
        wksph = create_history(
            self.request,
            name="NGPhylogeny Analyse - " + workflow.description,
            wf_category=cat,
            wf_steps=workflow.tooldesc,
        )

        if galaxy_file == "--":
            # we send the file to galaxy
            output = gi.tools.upload_file(path=tmp_file.name,
                                          file_name=uploadfile_name,
                                          history_id=wksph.history)
            galaxy_file = output.get('outputs')[0].get('id')

        dataset_map[i_input] = {'id': galaxy_file, 'src': 'hda'}

        # We analyze submited forms and upload files to
        # galaxy
        try:
            self.analyze_forms(request, context, workflow, params, gi, wksph, nseq, length, seqaa)
        except WorkflowInvalidFormError as e:
            # if one form is not valid
            workflow.delete_from_galaxy(gi)
            delete_history(request, wksph.history)
            return self.get(request, *args, **kwargs)
        except WorkflowInputFileFormatError as e:
            context = self.get_context_data(object=self.object)
            context['fileerror'] = str(e)
            workflow.delete_from_galaxy(gi)
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
            wksph.workflow = workflow
            wksph.save()

            return HttpResponseRedirect(self.succes_url)

        except Exception:
            # Real production bug: this used to call delete_history
            # (workspace/views.py, @connection_galaxy-decorated,
            # signature (request, history_id)) as delete_history(
            # wksph.history) - a single positional string arg, which
            # bound the history id string to the *request* parameter
            # instead. connection_galaxy's wrapper then crashed with
            # AttributeError: 'str' object has no attribute 'session' on
            # request.session.get(...) - caught by the decorator's own
            # broad except Exception (logged, HttpResponseGone
            # returned), so the AttributeError never surfaced directly,
            # but the cleanup silently never ran either: every failed
            # Advanced-workflow submission (e.g. a tool parameter Galaxy
            # rejects, like the real "randstart" out-of-range case that
            # surfaced this) left its just-created WorkspaceHistory row
            # and local session state uncleaned. tools/views.py's own
            # call site already used the correct (request, history_id)
            # form - matched here.
            delete_history(request, wksph.history)
            raise
        #finally:
            # delete the workflow copy of oneclick workflow when
            # the workflow has been run
            # workflow.delete(gi)


