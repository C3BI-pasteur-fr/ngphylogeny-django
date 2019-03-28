import random
import string

from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator

from galaxy.decorator import connection_galaxy
from tools.models import ToolFlag, Tool
from workflows.models import Workflow
from workflows.models import WorkflowGalaxyFactory

from workflows.views.wkadvanced import WorkflowAdvancedFormView

WORKFLOW_MAKER_FLAG = 'wmake'


@connection_galaxy
def workflows_alacarte_build(request):
    """
        Workflow maker display
         - GET: Display a list of tools regroup by category
         - POST:
            o Retrieve selected tools
            o Created Galaxy json Workflow file
            o Import workflow into Galaxy
            :return Galaxy id workflow [wk_id]
    """
    WORKFLOW_STATIC_STEPS = [{"step": 0, "category": ['algn', ], "group": []},
                             {"step": 1, "category": ['clean', ], "group": []},
                             {"step": 2, "category": ['tree', 'model'], "group": []},
                             {"step": 3, "category": ['visu', ], "group": []},
                             ]

    for step in WORKFLOW_STATIC_STEPS:
        flags = ToolFlag.objects.filter(name__in=step.get('category'))
        for flag in flags:
            step['group'].append({'flag': flag,
                                  'tools': Tool.objects.filter(
                                      galaxy_server=request.galaxy_server,
                                      visible=True,
                                      toolflag=flag
                                  )
                                 .filter(toolflag__name=WORKFLOW_MAKER_FLAG)
                                  })

    if request.method == 'POST':
        tools = []
        previous_tool = ""
        form_is_valid = True
        for step in WORKFLOW_STATIC_STEPS:
            category = step.get('category')
            for cat in category:
                select_tool_pk = request.POST.get(cat, None)
                if select_tool_pk:
                    break
            else:
                select_tool_pk = None
            # TODO use data input/output format to invalidate form
            if (("visu" in category) and
                    select_tool_pk and
                    (not previous_tool or not tools)):
                form_is_valid = False
                messages.add_message(
                    request,
                    messages.ERROR,
                    "This combination of tools is not allowed")
                break
            if select_tool_pk:
                tools.append(select_tool_pk)
            previous_tool = select_tool_pk
        if form_is_valid:
            dict_tools = Tool.objects.in_bulk(tools)
            # sort tool
            list_tool = [dict_tools.get(int(t)) for t in tools]
            if list_tool:
                gi = request.galaxy
                server = request.galaxy_server
                # create temp history to build workflow
                history_id = gi.histories.create_history(name="temp").get("id")
                # build workflow object
                description = request.POST.get('wkname')
                name = "Workflow["+"".join(random.sample(string.ascii_letters + string.digits, 8))+"]"
                if not description:
                    description = name
                category='automaker'
                slug=name
                wkg = WorkflowGalaxyFactory(category, name, description)
                if wkg.build(gi, list_tool, history_id):
                    if wkg.valid:
                        wk = wkg.to_ngworkflow(server, gi)
                        # remove temp history
                        gi.histories.delete_history(history_id, purge=True)
                        return redirect(reverse(
                            "workflow_maker_form", args=[wk.id_galaxy]))
                gi.histories.delete_history(history_id, purge=True)
    context = {"workflow": WORKFLOW_STATIC_STEPS}
    return render(request, 'workflows/workflows_alacarte.html', context)


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowMakerView(WorkflowAdvancedFormView):
    """
    Workflow form with the list of tools and launch workflow
    """
    template_name = "workflows/workflows_maker_form.html"
    restricted_toolset = Tool.objects.filter(
        toolflag__name=WORKFLOW_MAKER_FLAG)

    def get_workflow(self, detail=True):

        if detail:
            self.copy_workflow = self.get_workflow_detail()
        else:
            self.copy_workflow = self.get_object()

        return self.copy_workflow

    def get_object(self, queryset=None, detail=True):
        # load workflow
        wk_json = self.request.galaxy.workflows.show_workflow(
            workflow_id=self.kwargs['id'])
        wkname = wk_json.get('name')

        # create workflow
        wk_obj = Workflow.objects.filter(id_galaxy=self.kwargs['id']).exclude(category='base').first()
        # add galaxy json information
        wk_obj.json = wk_json
        wk_obj.save()
        return wk_obj
