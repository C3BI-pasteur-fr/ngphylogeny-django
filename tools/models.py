from __future__ import unicode_literals

from django.db import models

# Create your models here.

class Tool(models.Model):
    """

    """
    id_galaxy = models.CharField(max_length=250, unique=True )
    name = models.CharField(max_length=250)
    #params =

    @property
    def input_data(self):
        "return the list of input of the tool"
        return ToolData.objects.filter(tool=self.id, type='i')

    @property
    def output_data(self):
        "return the list of output of the tool"
        return ToolData.objects.filter(tool=self.id, type='o')

    @property
    def compatible_tool(self):
        """
            return the list of tools that have inputs compatible with the outputs of the tool
            to create workflow.

            ie: tool.output.format = next_tool.input.format
        """
        outputs = self.output_data
        outputs_formats = set([o.edam_formats for o in outputs])

        return Tool.objects.filter(tooldata__edam_formats__in=outputs_formats, tooldata__type="i")


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
    type = models.CharField(max_length=1,choices=DATA_TYPE_CHOICES, default='i')
    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)

    class Meta:
        verbose_name_plural = "Tool_data"


    def __unicode__(self):
        return "%s_%s"%(self.name, self.tool)