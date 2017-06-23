from django import forms
from django.contrib import admin
from django.contrib import messages

from .models import Tool, ToolInputData, ToolOutputData, ToolFlag


class CustomMultipleChoiceField(forms.ModelMultipleChoiceField):
    def __init__(self, *args, **kwargs):
        super(CustomMultipleChoiceField, self).__init__(queryset=None, *args, **kwargs)
        self.required = False
        self.widget = forms.SelectMultiple(attrs={'disabled': 'disabled'})


class ToolForm(forms.ModelForm, forms.Form):
    compatible_tool = CustomMultipleChoiceField()

    class Meta:
        model = Tool
        fields = ('name', 'id_galaxy')
        readonly_fields = ('compatible_tool')

    def __init__(self, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)
        self.fields['compatible_tool'].queryset = kwargs.get('instance').compatible_tool


class ToolOutputDataInline(admin.TabularInline):
    model = ToolOutputData
    extra = 0


class ToolInputDataInline(admin.TabularInline):
    model = ToolInputData
    extra = 0


class ToolFlagInline(admin.TabularInline):
    model = ToolFlag.tool.through
    extra = 0


class ToolOutputDataAdmin(admin.ModelAdmin):
    pass


class ToolInputDataAdmin(admin.ModelAdmin):
    pass


class ToolAdmin(admin.ModelAdmin):
    """

    """

    list_display = ['name', 'galaxy_server', 'toolshed', 'version', 'visible', 'toolflags']
    list_filter = ['galaxy_server', 'toolshed']
    fields = ('galaxy_server', 'toolshed', 'name', 'version', 'id_galaxy',)
    inlines = [
        ToolInputDataInline,
        ToolOutputDataInline,
        ToolFlagInline
    ]
    actions = ['search_more_tools_from_this_tool_server', 'hide_tools', 'display_tools']

    def search_more_tools_from_this_tool_server(self, request, queryset):
        """
            find more tools based on selected tool galaxy server
        """
        tools_list = []
        for tool in queryset:
            tools_list.extend(Tool.import_tools(tool.galaxy_server))
        self.message_user(request, "%s successfully imported new tools." % (len(tools_list)))

    def hide_tools(self, request, queryset):
        """
        hide selected tool from  user website interface
        """
        queryset.select_for_update().update(visible=False)

    def display_tools(self, request, queryset):
        """
        display selected tool from user website interface
        """
        queryset.select_for_update().update(visible=True)

    def get_actions(self, request):
        actions = super(ToolAdmin, self).get_actions(request)
        """add flag assignation in action Tool"""
        for flag in ToolFlag.objects.all():
            action = make_actions_addflag(flag)
            actions[action.__name__] = (action,
                                        action.__name__,
                                        action.short_description)

        return actions


def make_actions_addflag(flag):
    """ Create functions to add a flag on selected tools
    """

    def assign_flag(modeladmin, request, queryset):
        for tool in queryset:
            tool.toolflag_set.add(flag)
            messages.info(request, "Flag {0} assigned to {1}".format(flag.verbose_name,
                                                                     tool.name
                                                                     ))

    assign_flag.short_description = 'Add {0}'.format(flag.verbose_name)
    assign_flag.__name__ = 'add_{0}_flag'.format(flag.name)
    return assign_flag


admin.site.register(Tool, ToolAdmin)
admin.site.register(ToolInputData, ToolInputDataAdmin)
admin.site.register(ToolOutputData, ToolOutputDataAdmin)
admin.site.register(ToolFlag)
