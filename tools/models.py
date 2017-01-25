import urllib
import requests
from django.db import models

from account.models import GalaxyServer

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

    galaxy_server = models.ForeignKey(GalaxyServer, null=True, blank=True)
    id_galaxy = models.CharField(max_length=250, unique=True)
    toolshed = models.CharField(max_length=100, blank=True)
    name = models.CharField(max_length=100, blank=True)
    version = models.CharField(max_length=20, blank=True)
    description = models.CharField(max_length=250)

    @classmethod
    def import_tools_from_url(cls, galaxy_url, query="phylogeny"):
        try:
            galaxy_server = GalaxyServer.objects.get(url=galaxy_url)
        except:
            raise Exception("NGPhylogeny server is not properly configured,"
                            "Please ensure that the Galaxy server is correctly set up")
        return cls.import_tools(galaxy_server)

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

                    if "toolshed" in t.id_galaxy:
                        t.toolshed = t.id_galaxy.split('/')[0]
                    t.save()

                    # save inputs
                    inputs_tools = tool_info.get('inputs')
                    for input_d in inputs_tools:

                        if input_d.get('type') == 'data':

                            edam = input_d.get('edam')
                            ed_format = ""
                            if edam:
                                ed_format = edam.get('edam_formats')

                            input_obj = ToolDataInput(name=input_d.get('name'),
                                                      edam_formats=ed_format,
                                                      extensions=input_d.get('extensions'),
                                                      tool=t
                                                      )
                            input_obj.save()

                    # save outputs
                    outputs_tools = tool_info.get('outputs')
                    for output_d in outputs_tools:
                        output_obj = ToolDataOutput(name=output_d.get('name'),
                                                    edam_formats=output_d.get('edam_format'),
                                                    extensions=output_d.get('format'),
                                                    type="o",
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
        outputs = ToolDataOutput.objects.filter(tool=self)
        outputs_formats = outputs.values_list('edam_formats', flat=True)

        tools_data_input = ToolDataInput.objects.filter(tool=self)
        tools_compatible = []

        from ast import literal_eval
        for data_input in tools_data_input:

            list_edam_formats = literal_eval(data_input.edam_formats)
            for ed in list_edam_formats:
                if ed in set(outputs_formats):
                    tools_compatible.append(data_input.tool.pk)

        return Tool.objects.filter(pk__in=tools_compatible)

    def __unicode__(self):
        return self.name


class ToolData(models.Model):
    """
        Generic model for Input/output data
    """
    DATA_TYPE_CHOICES = (
        ('i', "input"),
        ('o', "output"),
    )
    name = models.CharField(max_length=250)
    format = models.CharField(max_length=100, blank=True)
    extensions = models.CharField(max_length=100)
    edam_formats = models.CharField(max_length=250, null=True)
    type = models.CharField(max_length=1, choices=DATA_TYPE_CHOICES, default='i')
    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)

    class Meta:
        verbose_name_plural = "Tool Data"

    def __unicode__(self):
        return "%s_%s" % (self.name, self.tool)


class ToolDataOutputManager(models.Manager):
    def get_queryset(self):
        return super(ToolDataOutputManager, self).get_queryset().filter(
            type='o')


class ToolDataInputManager(models.Manager):
    def get_queryset(self):
        return super(ToolDataInputManager, self).get_queryset().filter(
            type='i')


class ToolDataInput(ToolData):
    objects = ToolDataInputManager()

    class Meta:
        proxy = True
        default_related_name = 'data_input'


class ToolDataOutput(ToolData):
    objects = ToolDataOutputManager()

    class Meta:
        proxy = True
        default_related_name = 'data_output'


class ToolFlag(models.Model):
    name = models.CharField(max_length=5, unique=True)
    verbose_name = models.CharField(max_length=250, unique=True)
    tool = models.ManyToManyField(Tool)

    def __unicode__(self):
        return self.verbose_name


class ToolInterconnections(models.Model):
    data = models.ForeignKey(ToolDataOutput, related_name='data_output_producted')

    @property
    def producted_by(self):
        return self.tool_output.tool

    @property
    def compatible_with(self):
        return self.next_tool_input.tool

    input_field = models.ForeignKey(ToolDataInput, related_name='data_input_compatible')

    class Meta:
        verbose_name_plural = "Tool Interconnections"
