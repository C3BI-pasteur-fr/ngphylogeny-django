import tempfile

from bioblend.galaxy.tools.inputs import inputs
from django.core.urlresolvers import reverse_lazy
from django.http import JsonResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.views.generic import DetailView, ListView

from account.decorator import connection_galaxy
from .forms import ToolForm
from .models import Tool


class ToolListView(ListView):
    model = Tool


class ToolDetailView(DetailView):
    model = Tool


class ToolJSONView(ToolDetailView):
    def get(self, request, *args, **kwargs):
        tool = self.get_object()
        jsondata = {'id': tool.id, 'name': tool.name,
                    'id_galaxy': tool.id_galaxy,
                    'inputs_data': list(str(i) for i in tool.input_data),
                    'outputs_data': list(str(o) for o in tool.output_data),
                    'tool_linkable': list(str(t) for t in tool.compatible_tool),
                    }
        return JsonResponse(jsondata)


@connection_galaxy
def tool_form_view(request, pk):
    tool_obj = get_object_or_404(Tool, pk=pk)
    data = request.galaxy.tools.show_tool(tool_id=tool_obj.id_galaxy, io_details='true')

    context = {"toolform": ToolForm(data['inputs']),
               "tool": tool_obj}

    return render(request, 'tools/tool_form.html', context)


@connection_galaxy
def tool_exec_view(request, pk):
    gi = request.galaxy
    message = ""

    tool_obj = get_object_or_404(Tool, pk=pk)
    tool_inputs_details = gi.tools.show_tool(tool_id=tool_obj.id_galaxy, io_details='true')
    tool_form = ToolForm(tool_inputs_details['inputs'], request.POST or None)

    if request.method == 'POST':

        if tool_form.is_valid():

            try:
                history_id = gi.histories.get_most_recently_used_history().get('id')

            except:
                # Create a new galaxy history and delete older if the user is not authenticated
                history_id = gi.histories.create_history().get('id')

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

            try:
                tool_outputs = gi.tools.run_tool(history_id=history_id,
                                                 tool_id=tool_obj.id_galaxy,
                                                 tool_inputs=tool_inputs)

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
               }

    return render(request, 'tools/tool_form.html', context)
