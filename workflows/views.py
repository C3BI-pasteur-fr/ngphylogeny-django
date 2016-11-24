from django.shortcuts import render, redirect
from django.utils.decorators import method_decorator
from django.views.generic import TemplateView
from formtools.wizard.views import SessionWizardView
from django.core.urlresolvers import reverse_lazy

from account.decorator import connection_galaxy
from tools.forms import ToolForm
from tools.models import Tool,ToolFlag
from tools.views import tool_exec


@method_decorator(connection_galaxy, name="dispatch")
class HistoryDetailView(TemplateView):
    """
        Display Galaxy like history information
    """
    template_name = 'workspace/history.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(HistoryDetailView, self).get_context_data(**kwargs)

        gi = self.request.galaxy

        history_content = gi.histories.show_history(context['history_id'], contents=True)

        context['history_info'] = gi.histories.show_history(context['history_id'])
        context['history_content'] = history_content

        return context


def form_class_list(gi, ids_tools):
    """
    Create ToolForm classes on the fly to be used by WizardForm
    :param gi:
    :param ids_tools:
    :return: list af Class form
    """
    tools = Tool.objects.filter(pk__in=ids_tools)
    tools_inputs_details = []
    for tool in tools:
        tool_inputs_details = gi.tools.show_tool(tool_id=tool.id_galaxy, io_details='true')
        tools_inputs_details.append(
            type(str(tool.name) + 'Form',
                 (ToolForm,),
                 {'tool_params': tool_inputs_details.get('inputs'),'tool_id': tool.id_galaxy }
                 )
        )
    return tools_inputs_details


@connection_galaxy
def workflow_form(request):

    tools = Tool.objects.all()
    form_list = form_class_list(request.galaxy, tools)
    selected_tools = request.session.get('selected_tools')

    if selected_tools:
         form_list = form_class_list(request.galaxy, selected_tools)

    return WorkflowWizard.as_view(form_list=form_list)(request)


class WorkflowWizard(SessionWizardView):

    template_name = 'workflows/workflows_form.html'

    def done(self, form_list, **kwargs):
        for form in form_list:
            #launch tools
            output = tool_exec(self.request, form)
            print output

        return redirect(reverse_lazy("history_detail", kwargs={'history_id': output }))



@connection_galaxy
def workflows_advanced_mode_build(request):

    WORKFLOW_ADVANCED_MODE = [{"step": 0, "category": 'algn', "group": ""},
                              {"step": 1, "category": 'clean', "group": ""},
                              {"step": 2, "category": 'tree', "group": ""},
                              {"step": 3, "category": 'visu', "group": ""},
                              ]

    for step in WORKFLOW_ADVANCED_MODE:
        step['group'] = ToolFlag.objects.get(name=step.get('category'))

    if request.method == 'POST':

        tools = []
        for step in WORKFLOW_ADVANCED_MODE:
            tools.append(request.POST.get(step.get('category')))

        request.session['selected_tools'] = tools
        return redirect('workflow_form')

    context = {"workflow": WORKFLOW_ADVANCED_MODE}
    return render(request, 'workflows/workflows_advanced.html', context)



def workflows_oneclick_mode_build(request):

    pass

def workflows_alacarte_mode_build(request):

    WORKFLOW_ADVANCED_MODE = [{"step": 0, "category": 'algn', "group": ""},
                              {"step": 1, "category": 'clean', "group": ""},
                              {"step": 2, "category": 'tree', "group": ""},
                              {"step": 3, "category": 'visu', "group": ""},
                              ]

    for step in WORKFLOW_ADVANCED_MODE:
        step['group'] = ToolFlag.objects.get(name=step.get('category'))
    context = {"workflow": WORKFLOW_ADVANCED_MODE}

    return render(request, 'workflows/workflows_alacarte.html', context)
