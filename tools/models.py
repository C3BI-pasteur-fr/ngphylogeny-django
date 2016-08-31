from __future__ import unicode_literals
from django.db import models


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
    id_galaxy = models.CharField(max_length=250, unique=True)
    name = models.CharField(max_length=100)
    description = models.CharField(max_length=250)

    @property
    def input_data(self):
        """return the list of input of the tool"""
        return self.tooldata_set.filter(type='i')


    @property
    def output_data(self):
        """return the list of output of the tool"""
        return self.tooldata_set.filter(type='o')


    @property
    def compatible_tool(self):
        """
            return the list of tools that have inputs compatible with the outputs of the tool
            to create workflow.

            ie: tool.output.format == next_tool.input.format
        """
        outputs = self.output_data
        outputs_formats = outputs.values_list('edam_formats', flat=True)

        tools_data_input = ToolData.objects.filter(type="i").select_related()
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

    """
    DATA_TYPE_CHOICES = (
        ('i', "input"),
        ('o', "output"),
    )
    name = models.CharField(max_length=250)
    format = models.CharField(max_length=25, blank=True)
    extensions = models.CharField(max_length=25)
    edam_formats = models.CharField(max_length=250, null=True)
    type = models.CharField(max_length=1, choices=DATA_TYPE_CHOICES, default='i')
    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)

    class Meta:
        verbose_name_plural = "Tool_data"

    def __unicode__(self):
        return "%s_%s" % (self.name, self.tool)


class ToolFlag(models.Model):

    name = models.CharField(max_length=5, unique=True)
    verbose_name = models.CharField(max_length=250, unique=True)
    tool = models.ManyToManyField(Tool)

    def __unicode__(self):
        return self.verbose_name