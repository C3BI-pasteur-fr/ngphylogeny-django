from django.utils.decorators import method_decorator
from django.utils.text import slugify
from django.shortcuts import render
from django.views.generic import View
from django.views.generic.detail import SingleObjectMixin
from django.core.files.uploadedfile import InMemoryUploadedFile

import tempfile

from data.views import UploadView
from galaxy.decorator import connection_galaxy
from tools.models import Tool

from workspace.views import create_history, delete_history

from workflows.forms import tool_form_factory
from workflows.views.generic import WorkflowWizard, WorkflowListView
from workflows.views.viewmixing import WorkflowDuplicateMixin, WorkflowDetailMixin

from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponseRedirect

WORKFLOW_ADV_FLAG = "wadv"


def form_class_list(tools):
    """
    Create ToolForm classes on the fly to be used by WizardForm
    :param gi:
    :param tools:
    :return: list af Class form
    """
    tool_form_class = []
    for tool in tools:
        tool_form_class.append(tool_form_factory(tool))
    return tool_form_class


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowAdvancedListView(WorkflowListView):
    """
        Workflow Advanced ListView
    """
    template_name = 'workflows/workflows_advanced_list.html'
    restricted_toolset = Tool.objects.filter(toolflag__name=WORKFLOW_ADV_FLAG)


class WorkflowAdvancedFormView(WorkflowDetailMixin, SingleObjectMixin):
    """
        Generic Workflow Advanced class-based view:
        - Cut selected workflow into multiple tool forms
    """
    # object = workflow

    template_name = 'workflows/workflows_advanced_form.html'
    context_object_name = "workflow"
    restricted_toolset = Tool.objects.filter(toolflag__name=WORKFLOW_ADV_FLAG)

    def get_context_data(self, **kwargs):

        if not self.object:
            self.object = self.get_object()
        self.fetch_workflow_detail(self.object)
        context = super(WorkflowAdvancedFormView, self).get_context_data(**kwargs)
        context['workflow_list'] = [self.object, ]

        context['tool_list'] = []
        tools = []

        for t in self.object.detail:
            tools.append(t[1])
            context['tool_list'].append(slugify(t[1].name))

        context['form_list'] = form_class_list(tools)

        return context


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowAdvancedSinglePageView(WorkflowDuplicateMixin, WorkflowAdvancedFormView, View):
    template_name = 'workflows/workflows_adv_singlepage_form.html'

    def get(self, request, *args, **kwargs):
        context = self.get_context_data(object=self.object)

        return render(request, self.template_name, context)
    
    
    def post(self, request, *args, **kwargs):

        gi = request.galaxy
        message = ""

        workflow = self.get_workflow()
        context = self.get_context_data(object=self.object)
        # input file
        dataset_map = {}

        # tool params
        params = {}

        i_input = workflow.json['inputs'].keys()[0]
        steps = workflow.json['steps']
        step_id = u'0'

        history_id = create_history(self.request, name="NGPhylogeny Analyse - " + workflow.name)

        for form in context['form_list']:

            tool_form = form(data=request.POST, prefix=form.prefix)

            if not tool_form.is_valid():
                # if form is not valid
                return self.get(request, *args, **kwargs)
            else:

                print (tool_form.cleaned_data)
                tool_inputs = inputs()

                # mapping between form id to galaxy params names
                fields = tool_form.fields_ids_mapping
                inputs_data = set(tool_form.input_file_ids)

                # set the Galaxy parameter (name, value)
                for key, value in tool_form.cleaned_data.items():
                    if not key in inputs_data:
                        tool_inputs.set_param(fields.get(key), value)

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
                        outputs = gi.tools.upload_file(path=tmp_file.name, file_name=uploaded_file.name,
                                                       history_id=history_id)
                        file_id = outputs.get('outputs')[0].get('id')

                        tool_inputs.set_dataset_param(fields.get(inputfile), file_id)

                    else:
                        # else paste content
                        content = tool_form.cleaned_data.get(inputfile)
                        if content:
                            tmp_file = tempfile.NamedTemporaryFile()
                            tmp_file.write(content)
                            tmp_file.flush()

                            # send file to galaxy
                            input_fieldname = tool_form.fields_ids_mapping.get(inputfile)
                            outputs = gi.tools.upload_file(path=tmp_file.name,
                                                           file_name=input_fieldname + " pasted_sequence",
                                                           history_id=history_id)
                            file_id = outputs.get('outputs')[0].get('id')
                            tool_inputs.set_dataset_param(fields.get(inputfile),
                                                          file_id)

                # workflow step
                # get from which step the tools are used
                for i, step in steps.items():
                    if getattr(tool_form, 'tool_id', 'null') == step.get('tool_id'):
                        step_id = i
                        break

                # convert inputs to dict
                params[step_id] = tool_inputs.to_dict()
                if not params[step_id]:
                    del params[step_id]
                print (params)

        # main input file worklfow
        uploaded_file = request.FILES.get("file") or request.POST.get("file")

        if isinstance(uploaded_file, InMemoryUploadedFile):

            tmp_file = tempfile.NamedTemporaryFile()
            for chunk in uploaded_file.chunks():
                tmp_file.write(chunk)
            tmp_file.flush()
            uploadfile_name = str(uploaded_file.name)

        else:
            uploadfile_name = "Upload File"
            tmp_file = tempfile.NamedTemporaryFile()
            tmp_file.write(uploaded_file)
            tmp_file.flush()

        if uploaded_file:
            # send file to galaxy

            output = gi.tools.upload_file(path=tmp_file.name,
                                          file_name=uploadfile_name,
                                          history_id=history_id)

            galaxy_file = output.get('outputs')[0].get('id')
            dataset_map[i_input] = {'id': galaxy_file, 'src': 'hda'}

        try:
            output = self.request.galaxy.workflows.invoke_workflow(workflow_id=workflow.id_galaxy,
                                                                   history_id=history_id,
                                                                   inputs=dataset_map,
                                                                   params=params,
                                                                   allow_tool_state_corrections=True,
                                                                   )

            self.succes_url = reverse_lazy("history_detail", kwargs={'history_id': history_id})
            return HttpResponseRedirect(self.succes_url)

        except Exception:
            delete_history(history_id)
            raise
        finally:
            # delete the workflow copy of oneclick workflow when the workflow has been run
            print (self.clean_copy())
