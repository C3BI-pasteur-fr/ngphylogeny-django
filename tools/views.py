from django.http import JsonResponse
from django.views.generic import DetailView, ListView

from models import Tool

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
