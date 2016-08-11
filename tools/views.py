from django.http import JsonResponse
from django.shortcuts import render ,get_object_or_404
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
