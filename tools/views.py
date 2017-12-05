import json
import tempfile

from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.views.generic import DetailView, ListView

from galaxy.decorator import connection_galaxy
from workspace.views import create_history, delete_history
from .forms import ToolForm
from .models import Tool, ToolFieldWhiteList, ToolFlag


class ToolListView(ListView):
    """

    """
    CATEGORY = ['algn',
                'clean',
                'tree',
                'visu',
                ]

    def get_queryset(self):
        CATEGORY = ToolFlag.objects.filter(rank=0).values_list('name', flat=True) or self.CATEGORY
        return Tool.objects.filter(galaxy_server__current=True, visible=True).filter(toolflag__name__in=CATEGORY)


class ToolDetailView(DetailView):
    queryset = Tool.objects.filter(galaxy_server__current=True, visible=True)


@connection_galaxy
def tool_exec_view(request, pk, store_output=None):
    """
    :param request:
    :param pk:
    :param store_output:
    :return:
    """

    gi = request.galaxy
    message = ""

    tool_obj = get_object_or_404(Tool, pk=pk)

    toolfield, created = ToolFieldWhiteList.objects.get_or_create(tool=tool_obj, context='t')
    tool_visible_field = toolfield.saved_params

    tool_inputs_details = gi.tools.show_tool(tool_id=tool_obj.id_galaxy, io_details='true')
    tool_form = ToolForm(tool_params=tool_inputs_details['inputs'], tool_id=pk, whitelist=tool_visible_field,
                         data=request.POST or None)

    if request.method == 'POST':

        if tool_form.is_valid():

            history_id = create_history(request, name="Analyse with " + tool_obj.name)

            tool_inputs = inputs()
            for key, value in request.POST.items():
                tool_inputs.set_param(tool_form.fields_ids_mapping.get(key), value)

            if request.FILES:
                for input_file_id in request.FILES:

                    uploaded_file = request.FILES.get(input_file_id)
                    tmp_file = tempfile.NamedTemporaryFile()
                    for chunk in uploaded_file.chunks():
                        tmp_file.write(chunk)
                    tmp_file.flush()

                    # send file to galaxy
                    outputs = gi.tools.upload_file(path=tmp_file.name, file_name=uploaded_file.name,
                                                   history_id=history_id)
                    file_id = outputs.get('outputs')[0].get('id')
                    tool_inputs.set_dataset_param(tool_form.fields_ids_mapping.get(input_file_id.strip('[]')), file_id)

            else:
                # else paste content
                for input_file_id in set(tool_form.input_file_ids):

                    if input_file_id in request.POST.keys():
                        content = request.POST.get(input_file_id)
                        if content:
                            tmp_file = tempfile.NamedTemporaryFile()
                            tmp_file.write(content)
                            tmp_file.flush()

                            # send file to galaxy
                            outputs = gi.tools.upload_file(path=tmp_file.name, file_name="pasted_sequence",
                                                           history_id=history_id)
                            file_id = outputs.get('outputs')[0].get('id')
                            tool_inputs.set_dataset_param(tool_form.fields_ids_mapping.get(input_file_id.strip('[]')),
                                                          file_id)

            try:

                tool_outputs = gi.tools.run_tool(history_id=history_id,
                                                 tool_id=tool_obj.id_galaxy,
                                                 tool_inputs=tool_inputs)

                if store_output:
                    request.session['output'] = tool_outputs

                return redirect(reverse_lazy("history_detail", kwargs={'history_id': history_id}))

            except Exception as e:
                import ast
                from django.forms import ValidationError

                message = ast.literal_eval(e.body)
                reverse_dict_field = {v: k for k, v in tool_form.fields_ids_mapping.items()}

                for field, err_msg in message.get('err_data').items():
                    tool_form.add_error(reverse_dict_field.get(field),
                                        ValidationError(err_msg, code='invalid'))

                delete_history(request, history_id)

    context = {"toolform": tool_form,
               "tool": tool_obj,
               "message": message,
               }

    return render(request, 'tools/tool_form.html', context)


@connection_galaxy
def get_tool_name(request):
    """
    Find the name of tool from an id_tool (used by AJAX)
    """
    context = dict()

    if request.POST:
        gi = request.galaxy
        toolid = request.POST.get('tool_id')

        if toolid:
            tool = gi.tools.get_tools(tool_id=toolid)[0]
            context.update({'tool_id': toolid, 'name': tool.get('name')})

    return HttpResponse(json.dumps(context), content_type='application/json')
