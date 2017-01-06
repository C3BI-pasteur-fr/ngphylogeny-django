from django import forms
from django.contrib import admin

# Register your models here.
from .models import *


class CustomMultipleChoiceField(forms.ModelMultipleChoiceField):

    def __init__(self, *args, **kwargs):
        super(CustomMultipleChoiceField, self).__init__(queryset=None,*args, **kwargs)
        self.required = False
        self.widget = forms.SelectMultiple(attrs={'disabled': 'disabled'})


class ToolForm(forms.ModelForm, forms.Form ):

    compatible_tool = CustomMultipleChoiceField()

    class Meta:
        model = Tool
        fields = ('name', 'id_galaxy')
        readonly_fields = ('compatible_tool')

    def __init__(self, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)
        self.fields['compatible_tool'].queryset = kwargs.get('instance').compatible_tool


class ToolDataOutputInline(admin.TabularInline):

    model = ToolDataOutput
    extra = 0


class ToolDataInputInline(admin.TabularInline):

    model = ToolDataInput
    extra = 0


class ToolFlagInline(admin.TabularInline):

    model = ToolFlag.tool.through
    extra = 0


class ToolInterConnectionInline(admin.TabularInline):

    model = ToolData
    extra = 0


class ToolAdmin(admin.ModelAdmin):

    fields = ('galaxy_server','toolshed','name','version','id_galaxy',)
    inlines = [
        ToolDataInputInline,
        ToolDataOutputInline,
        ToolFlagInline
    ]


admin.site.register(Tool, ToolAdmin)
admin.site.register(ToolData)
admin.site.register(ToolFlag)
admin.site.register(ToolInterConnections)
