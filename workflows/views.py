from django.views.generic import TemplateView
from django.utils.decorators import method_decorator
from account.decorator import connection_galaxy
from django.forms import formset_factory
from tools.forms import ToolForm
from tools.models import Tool

from django.shortcuts import render, redirect
from formtools.wizard.views import SessionWizardView


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


class WorflowsWizard(SessionWizardView):
    template_name = 'workflows/workflows_form.html'
    def done(self, form_list, **kwargs):
        return redirect('home')


@connection_galaxy
def workflows_form(request):

    def form_class_list(request):
        gi = request.galaxy
        tools = Tool.objects.all()
        tools_inputs_details = []
        for tool in tools:
            tool_inputs_details = gi.tools.show_tool(tool_id=tool.id_galaxy, io_details='true')
            tools_inputs_details.append(type(str(tool.name)+'Form', (ToolForm,), {'tool_params':tool_inputs_details.get('inputs')} ))
        return tools_inputs_details

    form_list = form_class_list(request)

    return WorflowsWizard.as_view(form_list=form_list)(request)



def workflows_build(request):

    tools = Tool.objects.all()

    context = {"tool_list": tools}

    return render(request, 'workflows/workflows_build.html', context )







