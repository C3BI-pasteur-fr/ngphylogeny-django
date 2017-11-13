import random
import string

from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import View

from galaxy.decorator import connection_galaxy
from tools.models import ToolFlag, Tool
from workflows.models import Workflow
from workflows.models import WorkflowGalaxyFactory
from workflows.views.generic import WorkflowWizard, UploadView, DetailView
from workflows.views.viewmixing import WorkflowDeleteWorkingCopyMixin
from workflows.views.wkadvanced import WorkflowAdvancedFormView
from workspace.views import create_history


@connection_galaxy
def workflows_alacarte_build(request):
    """
        Workflow maker display
         - GET: Display a list of tools regroup by category
         - POST:
            o Retrieve selected tools
            o Created Galaxy json Workflow file
            o Import into Galaxy
            :return wizard_form_wiew
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
                                  'tools': Tool.objects.filter(galaxy_server=request.galaxy_server,
                                                               toolflag=flag,
                                                               visible=True).filter(toolflag__name='wmake')
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
            if ("visu" in category) and select_tool_pk and (not previous_tool or not tools):
                form_is_valid = False
                messages.add_message(request, messages.ERROR, "This combination of tools is not allowed")
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
                history_id = create_history(request)

                # build workflow object
                wkg = WorkflowGalaxyFactory(list_tool, gi, history_id)
                wkg.name = request.POST.get('wkname')

                if not wkg.name:
                    wkg.name = "Workflow_generated_"
                    wkg.name += ''.join(random.sample(string.ascii_letters + string.digits, 8))

                # create the JSON of generated workflow
                # create new entry into Galaxy
                wkgi = gi.workflows.import_workflow_json(wkg.to_json())
                wk_id = wkgi.get('id')

                return redirect(reverse("workflow_maker_form", args=[wk_id]))
                # return HttpResponse(json.dumps(wkg.to_json()), content_type='application/json')

    context = {"workflow": WORKFLOW_STATIC_STEPS}

    return render(request, 'workflows/workflows_alacarte.html', context)


@method_decorator(connection_galaxy, name="dispatch")
class WorkflowMakerView(WorkflowAdvancedFormView, DetailView):
    """
    Workflow form with the list of tools and launch workflow
    """
    restrict_toolset = Tool.objects.filter(toolflag__name='wmake')

    def get_object(self, queryset=None, detail=True):
        # load workflow
        wk_json = self.request.galaxy.workflows.show_workflow(workflow_id=self.kwargs['id'])
        wkname = wk_json.get('name')

        # create workflow
        wk_obj = Workflow(galaxy_server=self.request.galaxy_server,
                          id_galaxy=self.kwargs['id'],
                          name=wkname,
                          category='automaker',
                          description="Auto Generated Workflow",
                          slug=wkname)

        # add galaxy json information
        wk_obj.json = wk_json
        return wk_obj


class WorkflowMarkerWizardView(WorkflowDeleteWorkingCopyMixin, WorkflowWizard, WorkflowMakerView):
    """
        Workflow maker Wizard Form
    """

    template_name = 'workflows/workflows_maker_form.html'

    def done(self, *args, **kwargs):
        render_wizard = super(WorkflowMarkerWizardView, self).done(*args, **kwargs)

        if self.succes_url:
            # delete the workflow when the workflow has been run
            self.delete_workflow(self.get_object().id_galaxy)
        return render_wizard


class WorkflowsMarkerRedirectView(WorkflowMakerView, View):
    """
        Redirect to WizardformView build with selected tools
    """

    def get(self, request, *args, **kwargs):
        context = self.get_context_data()
        context['form_list'] = [UploadView.form_class] + context['form_list']

        return WorkflowMarkerWizardView.as_view(form_list=context['form_list'])(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        return self.get(request, *args, **kwargs)
