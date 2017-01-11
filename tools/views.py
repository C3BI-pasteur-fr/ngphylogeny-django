import tempfile

from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy

from django.shortcuts import render, get_object_or_404, redirect
from django.views.generic import DetailView, ListView

from account.decorator import connection_galaxy
from .forms import ToolForm
from .models import Tool
from workspace.views import create_history


class ToolListView(ListView):
    queryset = Tool.objects.filter(galaxy_server__galaxyconf__active=True)


class ToolDetailView(DetailView):
    queryset = Tool.objects.filter(galaxy_server__galaxyconf__active=True)


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
    tool_inputs_details = gi.tools.show_tool(tool_id=tool_obj.id_galaxy, io_details='true')
    print tool_inputs_details
    tool_form = ToolForm(tool_params=tool_inputs_details['inputs'], tool_id=pk, data=request.POST or None)

    if request.method == 'POST':

        if tool_form.is_valid():

            history_id = create_history(request)

            tool_inputs = inputs()
            for key, value in request.POST.items():
                tool_inputs.set_param(tool_form.fieds_ids_mapping.get(key), value)

            if request.FILES:
                for input_file_id in request.FILES:

                    uploaded_file = request.FILES.get(input_file_id)
                    tmp_file = tempfile.NamedTemporaryFile()
                    for chunk in uploaded_file.chunks():
                        tmp_file.write(chunk)
                    tmp_file.flush()

                    # send file to galaxy
                    outputs = gi.tools.upload_file(tmp_file.name, history_id, file_name=uploaded_file.name)
                    file_id = outputs.get('outputs')[0].get('id')
                    tool_inputs.set_dataset_param(tool_form.fieds_ids_mapping.get(input_file_id.strip('[]')), file_id)

            else:
                for input_file_id in tool_form.input_file_ids:
                    if input_file_id in request.POST.keys():
                        tool_inputs.set_dataset_param(
                                                      tool_form.fieds_ids_mapping.get(input_file_id.strip('[]')),
                                                      request.POST.get(input_file_id))

            try:

                tool_outputs = gi.tools.run_tool(history_id=history_id,
                                                 tool_id=tool_obj.id_galaxy,
                                                 tool_inputs=tool_inputs)
                print "#4", tool_outputs
                if store_output:
                    request.session['output'] = tool_outputs

                return redirect(reverse_lazy("history_detail", kwargs={'history_id': history_id}))

            except Exception as e:
                import ast
                from django.forms import ValidationError

                message = ast.literal_eval(e.body)
                reverse_dict_field = {v: k for k, v in tool_form.fieds_ids_mapping.items()}

                for field, err_msg in message.get('err_data').items():
                    tool_form.add_error(reverse_dict_field.get(field),
                                        ValidationError(err_msg, code='invalid'))


    context = {"toolform": tool_form,
               "tool": tool_obj,
               "message":message,
               }

    return render(request, 'tools/tool_form.html', context)


@connection_galaxy
def tool_exec(request, tool_form, store_output=None):

        gi = request.galaxy
        if tool_form.is_valid():

            history_id = create_history(request)

            tool_inputs = inputs()
            for key, value in request.POST.items():
                tool_inputs.set_param(tool_form.fieds_ids_mapping.get(key), value)

            if request.FILES:
                for input_file_id in request.FILES:

                    uploaded_file = request.FILES.get(input_file_id)
                    tmp_file = tempfile.NamedTemporaryFile()
                    for chunk in uploaded_file.chunks():
                        tmp_file.write(chunk)
                    tmp_file.flush()

                    # send file to galaxy
                    outputs = gi.tools.upload_file(tmp_file.name, history_id, file_name=uploaded_file.name)
                    file_id = outputs.get('outputs')[0].get('id')
                    tool_inputs.set_dataset_param(tool_form.fieds_ids_mapping.get(input_file_id.strip('[]')), file_id)

            else:
                for input_file_id in tool_form.input_file_ids:
                    if input_file_id in request.POST.keys():
                        tool_inputs.set_dataset_param(
                                                      tool_form.fieds_ids_mapping.get(input_file_id.strip('[]')),
                                                      request.POST.get(input_file_id))

            try:
                tool_outputs = gi.tools.run_tool(history_id=history_id,
                                                 tool_id=tool_form.tool_id_galaxy,
                                                 tool_inputs=tool_inputs)

                if store_output:
                    request.session['output'] = tool_outputs



            except Exception as e:
                import ast
                from django.forms import ValidationError

                message = ast.literal_eval(e.body)
                reverse_dict_field = {v: k for k, v in tool_form.fieds_ids_mapping.items()}

                for field, err_msg in message.get('err_data').items():
                    tool_form.add_error(reverse_dict_field.get(field),
                                        ValidationError(err_msg, code='invalid'))
            return tool_outputs