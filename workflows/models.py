from __future__ import unicode_literals

import collections
import json

from django.db import models

from galaxy.models import Server
from tools.models import Tool, ToolOutputData
from datetime import datetime


class Workflow(models.Model):
    """
    Galaxy Workflow information
    """

    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE)
    id_galaxy = models.CharField(max_length=250, unique=True)
    name = models.CharField(max_length=100)
    category = models.CharField(max_length=100, blank=True)
    version = models.CharField(max_length=10, blank=True)
    description = models.CharField(max_length=250)
    slug = models.SlugField(max_length=100, unique=True)
    rank = models.IntegerField(default=999, help_text="Workflows order")
    # Date added
    date = models.DateTimeField(default=datetime.now, blank=True)
    deleted = models.BooleanField(default=False)
    tooldesc = models.CharField(max_length=250,blank=True, default="")
    # Json representation of the workflow
    json =  None
    # Details about Steps (WorkflowStepInformation.sorted_tool_list)
    detail = None
    
    def fetch_details(self, galaxyinstance, toolset=None):
        self.json = galaxyinstance.workflows.show_workflow(
            workflow_id=self.id_galaxy)
        self.detail = WorkflowStepInformation(
            self.json,
            tools=toolset).sorted_tool_list
        self.tooldesc = ""
        i=0
        for step, tool in self.detail:
            if i>0:
                self.tooldesc+=','
            self.tooldesc+= "%s.%s" % (step, tool.name)
            i+=1

    def duplicate(self, galaxyinstance):
        """
        Returns a copy of this workflow.
        It copies the workflow locally on the django db
        and remotely on the galaxy server
        """
        workflow_copy = None
        try:
            # Try to import directly the publicly shared workflow
            # as a new workflow
            wf_import = galaxyinstance.workflows.import_shared_workflow(
                self.id_galaxy)
        except Exception:
            # Otherwise try to export it as dict and
            # import it again as a new workflow
            wk_cp = galaxyinstance.workflows.export_workflow_dict(
                self.id_galaxy)
            wf_import = galaxyinstance.workflows.import_workflow_dict(wk_cp)

        # makes the copy
        workflow_copy = Workflow(galaxy_server=self.galaxy_server,
                                 id_galaxy=wf_import.get('id'),
                                 name=self.name,
                                 category='duplicated',
                                 description=self.description,
                                 slug=wf_import.get('id')+"_"+self.name+"_copy")
        return workflow_copy

    def delete_from_galaxy(self, galaxyinstance):
        """
        Deletes the given workflow from galaxy server
        """
        msg = galaxyinstance.workflows.delete_workflow(
            workflow_id=self.id_galaxy)
        return msg

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

    tool_queryset = Tool.objects.all()

    def get_tools(self):
        # get known tools
        return self.tool_queryset.filter(
            id_galaxy__in=self.toolset).prefetch_related('toolflag_set')

    def update_dict_tools(self):

        _tools = self.get_tools()

        # remove unknown tools
        for nbstep, step in self.steps_tooldict.items():
            for tool in _tools:
                if tool.id_galaxy in step.get('tool_idgalaxy'):
                    step['tool'] = tool
                    break
            else:
                del self.steps_tooldict[nbstep]

    def __init__(self, workflow_json, tools=None):
        """
        :param workflow_json: Galaxy workflow json
        :param tools <queryset>: limit the list of available tools
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
                self.steps_tooldict[step_id] = {
                    'tool': None,
                    'tool_idgalaxy': step.get('tool_id'),
                    'annotation': step.get('annotation'),
                    'params': step.get('tool_inputs')
                }

        # update dict with known tool and remove steps with unknown tools
        self.update_dict_tools()

        # sort steps
        ord_step = collections.OrderedDict.fromkeys(
            sorted(self.steps_tooldict.keys()))
        ord_step.update(self.steps_tooldict)

        self.steps_tooldict = ord_step
        self.sorted_tool_list = list(
            (k, v.get('tool')) for k, v in ord_step.iteritems() if v and 'tool' in v)


class WorkflowGalaxyFactory(object):
    """
     Take an ordered list of <tool> to construct Galaxy like workflow
     - simple linear workflow maker
     - if no link exists between two steps, self.valid = False
    """
    CONVERSION_TOOL_FLAG = 'conve'

    def __init__(self, category, name, description):

        self.a_galaxy_workflow = "true"
        self.annotation = ""
        self.name = name
        self.steps = dict()
        self.valid = False
        self.description = description
        self.id_galaxy = ""
        self.slug = self.name
        self.category = category
        
    def build(self, galaxy_instance, list_tools, history_id):
        self.set_steps(galaxy_instance, list_tools, history_id)
        if self.valid:
            print self.to_json()
            wkgi = galaxy_instance.workflows.import_workflow_json(self.to_json())
            wk_id = wkgi.get('id')
            self.id_galaxy = wk_id
        return self.valid

    def to_ngworkflow(self, galaxy_server, galaxy_instance):
        wk_obj = None
        if self.valid:
            wk_obj = Workflow(galaxy_server=galaxy_server,
                              id_galaxy=self.id_galaxy,
                              name=self.name,
                              category=self.category,
                              description=self.description,
                              slug=self.slug)
            wk_obj.fetch_details(galaxy_instance)
            wk_obj.save()
        return wk_obj
    
    def set_steps(self, galaxy_instance, list_tool, history_id):
        self.valid = True
        step = 0
        # Input data : First step
        wfi = WorkflowToolInformation()
        wfi.id = step
        self.steps[str(step)] = wfi
        # iterate on list of tool
        # construct workflow steps chain
        for tool in list_tool:
            step += 1
            previous_step = wfi
            wfi = WorkflowToolInformation(tool, galaxy_instance, history_id)
            wfi.set_id(step)
            self.steps[str(step)] = wfi
            tmpvalid = False
            # link output compatible from previous step
            for inputdata in tool.toolinputdata_set.all():
                if previous_step:
                    # i.e first step : We link input data only to
                    # fields marked as galaxy_input_data
                    if (previous_step.type == "data_input" and
                            inputdata.galaxy_input_data):
                        wfi.set_input_connections(
                            stepid=previous_step.id,
                            input_name=inputdata.name,
                            output_name="output"
                        )
                    if previous_step.type == 'tool':
                        # find outputs of tool in step-1 compatible
                        # with current tool input
                        compatible_outputs = ToolOutputData.objects.filter(
                            compatible_inputs=inputdata,
                            # current tool input
                            tool__id_galaxy=previous_step.tool_id
                            # previous step output
                        )
                        if compatible_outputs:
                            tmpvalid = True
                            # set_input_connections
                            for o in compatible_outputs:
                                wfi.set_input_connections(
                                    stepid=previous_step.id,
                                    input_name=inputdata.name,
                                    output_name=o.name
                                )
                            i_extensions = inputdata.get_extensions()
                            if (i_extensions and
                                    not (o.extension in i_extensions)):
                                first_ext = "".join(i_extensions[:1])
                                previous_step.convert_action(o.name, first_ext)
            if not tmpvalid and step > 1:
                self.valid = False

    def __repr__(self):
        return str(self.__dict__)

    def to_json(self):
        import ast
        print(str(self))
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
        self.post_job_actions = {}
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
        for o in tool.tooloutputdata_set.all():
            self.outputs.append(
                {
                    'name': o.name,
                    'type': o.extension
                }
            )

    def set_input_connections(self, input_name, stepid, output_name):
        self.input_connections.update(
            {input_name:
                {
                    'id': stepid,
                    'output_name': output_name
                }
             }
        )

    def set_tool_state(self, tool, gi, history_id):

        tool_build = gi.make_get_request(
            url=gi.tools.url + '/' + tool.id_galaxy + '/build',
            params=dict(history_id=history_id))

        tool_state = tool_build.json()['state_inputs']
        tool_state.update(__page__=0)
        self.tool_state = json.dumps(tool_state)

    def convert_action(self, output, newtype):
        # convert previous output in the right input format
        self.post_job_actions.update(
            {
                "ChangeDatatypeAction" + output: {
                    "action_arguments": {
                        "newtype": newtype
                    },
                    "action_type": "ChangeDatatypeAction",
                    "output_name": output
                }
            }
        )

    def __repr__(self):
        return str(self.__dict__)
