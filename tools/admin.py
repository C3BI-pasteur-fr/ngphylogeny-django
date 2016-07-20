from django.contrib import admin
from django import forms
# Register your models here.
from .models import *


class CustomMultipleChoiceField(forms.ModelMultipleChoiceField):
    def __init__(self, *args, **kwargs):
        super(CustomMultipleChoiceField, self).__init__(queryset=None,*args, **kwargs)
        self.required = False
        self.widget = forms.SelectMultiple(attrs={'disabled': 'disabled'})

class ToolForm(forms.ModelForm, forms.Form ):

    input_data = CustomMultipleChoiceField()
    output_data = CustomMultipleChoiceField()
    compatible_tool = CustomMultipleChoiceField()

    class Meta:
        model = Tool
        fields = ('name', 'id_galaxy' )
        readonly_fields = ('input_data', 'output_data', 'compatible_tool')

    def __init__(self, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)
        self.fields['input_data'].queryset = kwargs.get('instance').input_data
        self.fields['output_data'].queryset = kwargs.get('instance').output_data
        self.fields['compatible_tool'].queryset = kwargs.get('instance').compatible_tool


class ToolDataInline(admin.TabularInline):
    model = ToolData
    extra = 0


class ToolAdmin(admin.ModelAdmin):

    form = ToolForm
    fields = ('name', 'id_galaxy', 'input_data','output_data','compatible_tool')
    inlines = [
                ToolDataInline,
             ]


admin.site.register(Tool, ToolAdmin)
admin.site.register(ToolData)