import ast

import requests
from django.db import models

from data.models import ExampleFile
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
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE)
    id_galaxy = models.CharField(max_length=250)
    toolshed = models.CharField(max_length=100, null=True, blank=True)
    name = models.CharField(max_length=100, blank=True)
    version = models.CharField(max_length=20, blank=True)
    description = models.CharField(max_length=250)
    toolshed_revision = models.CharField(max_length=250, null=True, blank=True)
    visible = models.BooleanField(default=True, help_text="Display this tool on the user web interface")
    oneclick = models.BooleanField(default=False)

    @property
    def toolflags(self):
        return ",".join(self.toolflag_set.all().values_list('verbose_name', flat=True))

    @property
    def get_params_detail(self):
        return self.fetch_tool_json().get('inputs')

    class Meta:
        unique_together = (('galaxy_server', 'id_galaxy'),)

    @classmethod
    def import_tools_from_url(cls, galaxy_url, query="phylogeny"):
        try:
            galaxy_server = Server.objects.get(url=galaxy_url)
        except:
            raise Exception("NGPhylogeny server is not properly configured,"
                            "Please ensure that the Galaxy server is correctly set up")
        return cls.import_tools(galaxy_server, query)

    def fetch_tool_json(self):
        tool_url = '%s/%s/%s/%s' % (self.galaxy_server.url, 'api', 'tools', self.id_galaxy)
        tool_info_request = requests.get(tool_url, params={'io_details': "true"})
        return tool_info_request.json()

    def import_tool(self):
        """
        Extract tool information from Galaxy server
        """
        tool_info = self.fetch_tool_json()

        # save tool information
        self.name = tool_info.get('name')
        self.description = tool_info.get('description')
        self.version = tool_info.get('version')

        toolshed = tool_info.get('tool_shed_repository')
        if toolshed:
            self.toolshed = toolshed.get('tool_shed')
            self.toolshed_revision = toolshed.get('changeset_revision')

        self.save()

        # save inputs
        inputs_tools = tool_info.get('inputs')
        inputs_list = []
        for input_d in inputs_tools:

            if input_d.get('type') == 'data':

                edam = input_d.get('edam')
                ed_format = ""
                if edam:
                    ed_format = edam.get('edam_formats')

                input_obj, created_input = ToolInputData.objects.get_or_create(name=input_d.get('name'),
                                                                               tool=self
                                                                               )
                if created_input:
                    input_obj.edam_formats = ed_format
                    input_obj.extensions = input_d.get('extensions')
                    input_obj.save()

            else:
                # remove input data from whitelist for workflows
                inputs_list.append(input_d.get('name'))

        # save outputs
        outputs_tools = tool_info.get('outputs')
        for output_d in outputs_tools:
            output_obj, created_output = ToolOutputData.objects.get_or_create(name=output_d.get('name'),
                                                                              tool=self
                                                                              )

            if created_output:
                output_obj.edam_format = output_d.get('edam_format')
                output_obj.extension = output_d.get('format')
                output_obj.save()

        # pre-compute whitelist for workflows
        w, created_w = ToolFieldWhiteList.objects.get_or_create(tool_id=self.id, context='w')
        if created_w:
            w._params = ",".join(inputs_list)
            w.save()

    @classmethod
    def import_tools(cls, galaxy_server, tools=None, query="phylogeny", force=False):

        tools_ids = tools or []
        tools_created = []

        if not tools:
            # fetch tools list
            tools_url = '%s/%s/%s/' % (galaxy_server.url, 'api', 'tools')
            connection = requests.get(tools_url, params={'q': query})
            if connection.status_code == 200:
                tools_ids = connection.json() or []

        for id_tool in tools_ids:

            t, created = Tool.objects.get_or_create(id_galaxy=id_tool, galaxy_server=galaxy_server)

            if created:
                t.import_tool()
                tools_created.append(t)

            elif force:
                t.import_tool()
                tools_created.append(t)

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
        unique_together = ['tool', 'name']


class ToolInputData(ToolData):
    """
    Tool input data information extract from Galaxy
    """
    edam_formats = models.CharField(max_length=250, null=True, blank=True)
    extensions = models.CharField(max_length=100)
    examplefile = models.ForeignKey(ExampleFile, null=True, blank=True)

    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)

    def get_extensions(self):
        return ast.literal_eval(self.extensions)

    def search_compatible_outputs(self, ignore=None):
        """
        Search compatible data produced by others tools
        :param ignore : list of ignored extension
        :return: list of <ToolOutputData>
        """

        l_ext = self.get_extensions()
        l_ext_filtered = [ext for ext in l_ext if ext not in ([] or ignore)]
        print l_ext_filtered
        galaxy_server = self.tool.galaxy_server
        return ToolOutputData.objects.filter(extension__in=l_ext_filtered,
                                             tool__galaxy_server=galaxy_server)

    def __unicode__(self):
        return "%s | %s: %s" % (self.tool, self.name, self.extensions)

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
        return ToolInputData.objects.filter(extensions__contains=self.extension)

    def __unicode__(self):
        return "%s | %s: %s" % (self.tool, self.name, self.extension)

    class Meta:
        verbose_name_plural = "Tool output data"


class ToolFlag(models.Model):
    """
    Flag tools by category
    """

    name = models.CharField(max_length=5, unique=True)
    verbose_name = models.CharField(max_length=250, unique=True)
    tool = models.ManyToManyField(Tool)
    rank = models.IntegerField(default=999, help_text="flags order")

    def __unicode__(self):
        return self.verbose_name

    class Meta:
        ordering = ['rank', 'verbose_name']


class ToolFieldWhiteList(models.Model):
    """
    Models use to filter Galaxy render form
    """
    CONTEXT_CHOICES = (('w', 'Workflows'),
                       ('aw', 'Advanced Workflows'),
                       ('t', 'Tool Forms'),
                       ('at', 'Advanced Tool Forms'),
                       )

    DICT_CONTEXT_CHOICES = dict(CONTEXT_CHOICES)

    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)
    context = models.CharField(max_length=2, choices=CONTEXT_CHOICES, help_text='In which context use the whitelist')
    _params = models.TextField(max_length=1000, blank=True, default="")

    @property
    def saved_params(self):
        if self._params:
            return self._params.split(',')
        else:
            return []

    @property
    def get_params(self):
        l = []
        for p in self.tool.get_params_detail:
            l.append(p.get('name'))
            cond = p.get('test_param')
            if cond:
                l.append(cond.get('name'))

        return l

    @property
    def get_json_params(self):
        return self.tool.get_params_detail

    def __unicode__(self):
        return "%s ,%s" % (self.tool.name, self.DICT_CONTEXT_CHOICES.get(self.context))

    class Meta:
        unique_together = ('tool', 'context')
