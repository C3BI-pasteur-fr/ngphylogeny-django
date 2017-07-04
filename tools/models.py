import ast
import urllib

import requests
from django.db import models

from forms import ToolForm
from galaxy.models import Server

TOOL_TAGS = (
    ('', '----'),
    ('algn', 'Multiple Alignment'),
    ('clean', 'Alignment curation'),
    ('tree', 'Construction of phylogenetic tree'),
    ('visu', 'Visualisation of phylogenetic tree')
)



class Tool(models.Model):
    """

    """

    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE, null=True, blank=True)
    id_galaxy = models.CharField(max_length=250)
    toolshed = models.CharField(max_length=100, blank=True)
    name = models.CharField(max_length=100, blank=True)
    version = models.CharField(max_length=20, blank=True)
    description = models.CharField(max_length=250)
    toolshed_revision= models.CharField(max_length=250, blank=True)
    visible = models.BooleanField(default=True, help_text="display this tool on the user web interface")

    @property
    def toolflags(self):
        return ",".join(self.toolflag_set.all().values_list('verbose_name',flat=True))

    class Meta:
        unique_together = (('galaxy_server','id_galaxy'),)


    @classmethod
    def import_tools_from_url(cls, galaxy_url, query="phylogeny" ):
        try:
            galaxy_server = Server.objects.get(url=galaxy_url)
        except:
            raise Exception("NGPhylogeny server is not properly configured,"
                            "Please ensure that the Galaxy server is correctly set up")
        return cls.import_tools(galaxy_server, query)

    @classmethod
    def import_tools(cls, galaxy_server, query="phylogeny"):

        params = urllib.urlencode({'q': query}, True)
        tools_url = '%s/%s/%s/?%s' % (galaxy_server.url, 'api', 'tools', params)
        connection = requests.get(tools_url)
        tools_ids = []
        tools_created = []
        if connection.status_code == 200:
            tools_ids = connection.json()

            for id_tool in tools_ids:
                params = urllib.urlencode({'io_details': "true"}, True)
                tool_url = '%s/%s/%s/%s/?%s' % (galaxy_server.url, 'api', 'tools', id_tool, params)
                tool_info_request = requests.get(tool_url)
                tool_info = tool_info_request.json()

                t, created = Tool.objects.get_or_create(id_galaxy=id_tool, galaxy_server=galaxy_server)

                if created:
                    tools_created.append(t)

                    # save tool informations
                    t.name = tool_info.get('name')
                    t.description = tool_info.get('description')
                    t.version = tool_info.get('version')

                    toolshed = tool_info.get('tool_shed_repository')
                    if toolshed:
                        t.toolshed = toolshed.get('tool_shed')
                        t.toolshed_revision = toolshed.get('changeset_revision')

                    t.save()

                    # save inputs
                    inputs_tools = tool_info.get('inputs')
                    for input_d in inputs_tools:

                        if input_d.get('type') == 'data':

                            edam = input_d.get('edam')
                            ed_format = ""
                            if edam:
                                ed_format = edam.get('edam_formats')

                            input_obj = ToolInputData(name=input_d.get('name'),
                                                      edam_formats=ed_format,
                                                      extensions=input_d.get('extensions'),
                                                      tool=t
                                                      )
                            input_obj.save()

                    # save outputs
                    outputs_tools = tool_info.get('outputs')
                    for output_d in outputs_tools:
                        output_obj = ToolOutputData(name=output_d.get('name'),
                                                    edam_format=output_d.get('edam_format'),
                                                    extension=output_d.get('format'),
                                                    tool=t
                                                    )
                        output_obj.save()

            return tools_created

    @property
    def compatible_tool(self):
        """
            return the list of tools that have inputs compatible with the outputs of the tool
            to create workflow.

            ie: tool.output.format == next_tool.input.format
        """
        outputs = ToolOutputData.objects.filter(tool=self)
        outputs_formats = outputs.values_list('edam_formats', flat=True)

        tools_data_input = ToolInputData.objects.filter(tool=self)
        tools_compatible = []

        from ast import literal_eval
        for data_input in tools_data_input:

            list_edam_formats = literal_eval(data_input.edam_formats)
            for ed in list_edam_formats:
                if ed in set(outputs_formats):
                    tools_compatible.append(data_input.tool.pk)

        return Tool.objects.filter(pk__in=tools_compatible)

    def form_class(self, galaxy_server):
        """
        toolform class factory
        :param galaxy_server:
        :return: Form class
        """

        params = urllib.urlencode({'io_details': "true"}, True)
        tool_url = '%s/%s/%s/%s/?%s' % (galaxy_server.url, 'api', 'tools', self.id_galaxy, params)
        tool_info_request = requests.get(tool_url)
        tool_inputs_details = tool_info_request.json()

        return type(str(self.name) + 'Form',
                        (ToolForm,),
                         {'tool_params': tool_inputs_details.get('inputs'),
                        'tool_id': self.id_galaxy
                        }
                    )

    def __unicode__(self):
        return self.name


class ToolData(models.Model):
    """
        Generic model for Input/output data
    """

    name = models.CharField(max_length=250)
    edam_data = models.CharField(max_length=250, null=True, blank=True)

    class Meta:
        abstract = True


class ToolInputData(ToolData):
    """
    Tool input data information extract from Galaxy
    """
    edam_formats = models.CharField(max_length=250, null=True, blank=True)
    extensions = models.CharField(max_length=100)

    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)


    def get_extensions(self):
        return  ast.literal_eval(self.extensions)


    def search_compatible_outputs(self):
        """
        Search compatible data produced by others tools
        :return: list of <ToolOutputData>
        """
        l_ext = self.get_extensions()
        return ToolOutputData.objects.filter(extension__in=l_ext)



    def __unicode__(self):
        return "%s %s" % (self.tool, self.extensions)

    class Meta:
        verbose_name_plural = "Tool input data"



class ToolOutputData(ToolData):
    """
    Tool output data extract from Galaxy
    """
    edam_format = models.CharField(max_length=250, null=True, blank=True)
    extension = models.CharField(max_length=100)

    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)
    compatible_inputs = models.ManyToManyField(ToolInputData, blank=True)

    def search_compatible_inputs(self):
        return ToolInputData.objects(extensions__contains=self.extension)


    def __unicode__(self):
        return "%s %s" % (self.tool, self.extension)

    class Meta:
        verbose_name_plural = "Tool output data"


class ToolFlag(models.Model):
    """
    Help to Organize tools by category
    """

    name = models.CharField(max_length=5, unique=True)
    verbose_name = models.CharField(max_length=250, unique=True)
    tool = models.ManyToManyField(Tool)
    rank = models.IntegerField(default=999, help_text="order flag")

    def __unicode__(self):
        return self.verbose_name


    class Meta:
        ordering = ["rank",]







