from __future__ import unicode_literals

import collections
import json

from django.db import models

from galaxy.models import Server
from tools.models import Tool, ToolOutputData


class Workflow(models.Model):
    """
    Galaxy Workflow information
    """

    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE, null=True, blank=True)
    id_galaxy = models.CharField(max_length=250, unique=True)
    name = models.CharField(max_length=100, unique=True)
    category = models.CharField(max_length=100, blank=True)
    version = models.CharField(max_length=10, blank=True)
    description = models.CharField(max_length=250)
    slug = models.SlugField(max_length=100, unique=True)
    rank = models.IntegerField(default=999, help_text="Workflows order")

    class Meta:
        ordering = ["rank", ]


class WorkflowStepInformation(object):
    """
    Parse Galaxy Workflow json:

        :steps_tooldict: dictionary , {'step_id' : { 'tool_idgalaxy':..,
                                                      'annotation': ..,
                                                     'params':.. }
                                        }
        :sorted_tool_list: list of tuple, [(step_id, queryset.tool)..]
    """

    tool_queryset = Tool.objects.filter()

    def get_tools(self):
        # get known tools
        return self.tool_queryset.filter(id_galaxy__in=self.toolset).prefetch_related('toolflag_set')

    def update_dict_tools(self):
        # remove unknown tools
        for nbstep, step in self.steps_tooldict.items():

            for tool in self.get_tools():
                if step.get('tool_idgalaxy') == tool.id_galaxy:
                    step['tool'] = tool
                    break
            else:
                del self.steps_tooldict[nbstep]

    def __init__(self, workflow_json, tools=None):
        """
        :param workflow_json: Galaxy workflow json
        :param tools <queryset>: limits the list of available tools
        """
        if isinstance(tools, type(self.tool_queryset)) and tools:
            self.tool_queryset = tools

        self.workflow_json = workflow_json
        self.steps_tooldict = {}
        self.sorted_tool_list = []
        self.toolset = set()

        # parse galaxy workflow information
        for step_id, step in self.workflow_json.get('steps').items():
            if step.get('tool_id'):
                self.toolset.update([step.get('tool_id')])
                self.steps_tooldict[step_id] = {'tool': None,
                                                'tool_idgalaxy': step.get('tool_id'),
                                                'annotation': step.get('annotation'),
                                                'params': step.get('tool_inputs')
                                                }

        # update dict with known tool and remove steps with unknown tools
        self.update_dict_tools()

        # sort steps
        ord_step = collections.OrderedDict.fromkeys(sorted(self.steps_tooldict.keys()))
        ord_step.update(self.steps_tooldict)

        self.steps_tooldict = ord_step
        self.sorted_tool_list = list((k, v.get('tool')) for k, v in ord_step.iteritems() if v and 'tool' in v)


class WorkflowGalaxyFactory(object):
    """
     Take an ordered list of <tool> to construct Galaxy like workflow
     - simple linear workflow maker
    """

    def __init__(self, list_tools, gi=None, history_id=""):

        # self.uuid= ""
        # self.format-version = ""
        self.a_galaxy_workflow = "true"
        self.annotation = ""
        self.name = "Imported workflow"
        self.steps = dict()

        if gi and history_id:
            self.set_steps(list_tools, gi, history_id)

    def set_steps(self, list_tool, gi, history_id):

        # Input data first step
        self.steps["0"] = WorkflowToolInformation()
        self.steps["0"].id = 0

        # iterate on list of tool
        # contruct workflow steps chain
        for tool in list_tool:

            step = len(self.steps)
            self.steps[str(step)] = WorkflowToolInformation(tool, gi, history_id)
            self.steps[str(step)].set_id(step)

            for inputdata in tool.toolinputdata_set.filter():

                # TODO: improve algorithm with recursive function
                # output compatible from previous step
                previous_step = self.steps.get(str(step - 1), None)
                if previous_step:

                    if previous_step.type == "data_input":

                        self.steps.get(str(step)).set_input_connections(
                            stepid=step - 1,
                            input_name=inputdata.name,
                            output_name="output"
                        )

                    else:
                        compatible_outputs = ToolOutputData.objects.filter(compatible_inputs=inputdata,
                                                                           tool__id_galaxy=previous_step.tool_id
                                                                           ).values_list("name", flat=True)

                        if compatible_outputs:
                            # create steps
                            self.steps[str(step)] = WorkflowToolInformation(tool, gi, history_id)
                            self.steps[str(step)].set_id(step)

                            for o in compatible_outputs:
                                self.steps.get(str(step)).set_input_connections(
                                    stepid=step - 1,
                                    input_name=inputdata.name,
                                    output_name=o
                                )
                        else:
                            # Search format conversion tool
                            # Edam : operation_0335
                            compatible_output = ToolOutputData.objects.filter(compatible_inputs=inputdata,
                                                                              tool__toolflag__name='conv'
                                                                              ).first()
                            if compatible_output:
                                conv_tool = compatible_output.tool

                                # move step...
                                self.steps[str(step + 1)] = self.steps.pop(str(step))
                                self.steps[str(step + 1)].set_id(step + 1)
                                self.steps[str(step + 1)].set_input_connections(
                                    stepid=step,
                                    input_name=inputdata.name,
                                    output_name=compatible_output.name
                                )

                                # ... to insert convertion tool step
                                self.steps[str(step)] = WorkflowToolInformation(conv_tool, gi, history_id)
                                self.steps[str(step)].set_id(step)

                                conv_tool_input = conv_tool.toolinputdata_set.filter().first()
                                output_name = ToolOutputData.objects.filter(compatible_inputs=conv_tool_input,
                                                                            tool__id_galaxy=previous_step.tool_id
                                                                            )

                                self.steps[str(step)].set_input_connections(
                                    stepid=step - 1,
                                    input_name=conv_tool_input.name,
                                    output_name=output_name.first().name
                                )

    def __repr__(self):
        return str(self.__dict__)

    def to_json(self):
        import ast
        return ast.literal_eval(str(self))


class WorkflowToolInformation(object):
    """
    Take a Tool to make Galaxy workflow formated tool:
    tool <object>
    url: <server_galaxy>/api/tools
    history_id: galaxy history_id
    """

    def __init__(self, tool=None, gi=None, history_id=''):

        self.annotation = ""
        self.content_id = None

        self.label = None
        self.input_connections = {}
        self.inputs = []
        self.outputs = []
        self.position = {'left': 300, 'top': 300}
        self.tool_errors = None
        self.tool_id = None
        self.tool_version = None
        self.uuid = "None"
        self.workflow_outputs = []

        if tool:
            self.content_id = tool.id_galaxy
            self.name = tool.name
            self.post_job_actions = {}
            self.tool_id = tool.id_galaxy
            self.tool_version = tool.version
            self.tool_shed_repository = {}
            self.type = "tool"
            self.set_inputs(tool)
            self.set_outputs(tool)

            if gi and history_id:
                self.set_tool_state(tool, gi, history_id)

        else:
            self.inputs.append(
                {
                    'description': "input dataset",
                    'name': "Input Dataset"
                })
            self.name = "Input dataset"
            self.tool_state = "{\"name\": \"Input Dataset\"}"
            self.type = "data_input"

    def set_id(self, id):
        self.id = id
        self.position = {'left': 300 + (200 * int(id)), 'top': 300}

    def set_inputs(self, tool):

        for i in tool.toolinputdata_set.filter():
            self.inputs.append(
                {
                    'description': 'runtime parameter for tool ' + self.name,
                    'name': i.name
                }
            )

    def set_outputs(self, tool):

        for o in tool.tooloutputdata_set.filter():
            self.outputs.append(
                {
                    'name': o.name,
                    'type': o.extension
                }
            )

    def set_input_connections(self, input_name, stepid="", output_name=""):

        self.input_connections.update(
            {input_name:
                {
                    'id': stepid,
                    'output_name': output_name
                }
            }
        )

    def set_tool_state(self, tool, gi, history_id):

        tool_build = gi.make_get_request(url=gi.tools.url + '/' + tool.id_galaxy + '/build',
                                         params=dict(history_id=history_id))

        tool_state = tool_build.json()['state_inputs']
        tool_state.update(__page__=0)
        self.tool_state = json.dumps(tool_state)

    def __repr__(self):
        return str(self.__dict__)
