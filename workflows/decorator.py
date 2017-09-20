


from tools.forms import tool_form_factory

from models import Workflow, WorkflowStepInformation
from galaxy.models import GalaxyUser, Server


def form_class_list(galaxy_server, tools):
    """
    Create ToolForm classes on the fly to be used by WizardForm
    :param gi:
    :param tools:
    :return: list af Class form
    """
    tool_form_class = []
    for tool in tools:
        tool_form_class.append(tool_form_factory(tool, galaxy_server=galaxy_server))
    return tool_form_class



def initformlist(view_function, *args, **kwargs):
    """Initiating wizard form list"""

    def wrapper(request, *args, **kwargs):

        wk = Workflow.objects.filter().first()
        galaxy_server = Server.objects.get(current=True)
        gi = GalaxyUser.objects.get(galaxy_server__current=True, anonymous=True).get_galaxy_instance()
        workflow_json = gi.workflows.show_workflow(workflow_id=wk.id_galaxy)

        # parse galaxy workflows json information
        wk.detail = WorkflowStepInformation(workflow_json).sorted_tool_list

        tools = [ dict(id_galaxy=t[1].id_galaxy,name=t[1].name) for t in wk.detail]

        request.session['tool_list'] = tools

        return view_function(request, *args, **kwargs)

    return wrapper