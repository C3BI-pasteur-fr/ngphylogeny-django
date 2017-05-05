from __future__ import unicode_literals

from django.db import models

from galaxy.models import Server
from tools.models import Tool


class Workflow(models.Model):
    """
    Galaxy Workflow informations
    """
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE, null=True, blank=True)
    id_galaxy = models.CharField(max_length=250, unique=True)
    name = models.CharField(max_length=100, unique=True)
    category = models.CharField(max_length=100 , blank=True)
    version = models.CharField(max_length=10, blank=True)
    description = models.CharField(max_length=250)
    slug = models.SlugField(max_length=100,unique=True)


class WorkflowStepInformation(object):
    """
    Parse Galaxy Workflow json:
        :graph: {step_input: [step_output, step_output], ...}
        :params: {step_id: param_dict}
        :tool_list: list of tuple, [(step_id, queryset.tool)..]
        :sorted_tool_list: list of tuple, [(step_id, queryset.tool)..]
    """

    @staticmethod
    def __graph_search(graph, node):
        """Breadth-first search (BFS) is an algorithm
           :return: list of node
        """
        step_visited = []
        next_steps = [node]
        while next_steps:
            current_step = next_steps.pop(0)
            if not (current_step in step_visited):
                step_visited.append(current_step)
                next_steps.extend(graph.get(current_step, []))
        return step_visited


    @staticmethod
    def __sort_tools( sorted_step, unsorted_tool_list):
        """
        :param sorted_step: list of steps
        :param unsorted_tool_list: list of tuple (step, tool)
        :return: list of tuple sorted (step, tools)
        """

        sorted_tool_list = []
        for e, step in enumerate(sorted_step):
            for y, x in unsorted_tool_list:
                if step == y:
                    if x:
                        sorted_tool_list.append((e, x[0]))
        return sorted_tool_list

    def __init__(self, workflow_json):
        """
        :param workflow_json: Galaxy workflow json
        """
        self.workflow_json = workflow_json
        self.graph = {}
        self.params = {}
        self.tool_list = []
        self.sorted_tool_list =[]

        # get known tools
        query = Tool.objects.filter()

        # parse galaxy workflow information
        for step_id, step in self.workflow_json.get('steps').items():
            self.params[step_id] = step.get('tool_inputs')
            self.tool_list.append([step_id, query.filter(id_galaxy=step.get('tool_id'))])

            for input, step_output in step.get("input_steps", {}).items():
                self.graph.setdefault(str(step_output.get('source_step')), []).append(str(step_id))


        first_step = self.workflow_json['inputs'].keys()[0]
        sorted_step = self.__graph_search( self.graph,first_step )
        self.sorted_tool_list = self.__sort_tools(sorted_step, self.tool_list)