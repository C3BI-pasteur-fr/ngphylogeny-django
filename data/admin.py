from django import forms
from django.contrib import admin
from django.contrib.admin.widgets import FilteredSelectMultiple

from tools.models import ToolInputData
from .models import ExampleFile


class ToolInputDataForm(forms.ModelForm):
    tool = forms.ModelMultipleChoiceField(queryset=ToolInputData.objects.all(),
                                          widget=FilteredSelectMultiple(verbose_name=ToolInputData._meta.verbose_name,
                                                                        is_stacked=False))


class ExampleFileAdmin(admin.ModelAdmin):
    """
    """
    form = ToolInputDataForm
    list_display = ('name', 'ext')
    fields = ['name', 'ext', 'type', 'upload', 'comment', 'tool']

    def get_form(self, request, obj=None, **kwargs):

        params_choices = ''
        if obj:
            params_choices = obj.toolinputdata_set.values_list('pk', flat=True)
        form = super(ExampleFileAdmin, self).get_form(request, obj, **kwargs)
        form.base_fields['tool'].initial = params_choices
        return form

    def save_model(self, request, obj, form, change):

        obj.save()
        obj.toolinputdata_set.clear()
        tools = form.cleaned_data.get('tool')
        for t in tools:
            t.examplefile = obj
            t.save()

        super(ExampleFileAdmin, self).save_model(request, obj, form, change)


admin.site.register(ExampleFile, ExampleFileAdmin)
