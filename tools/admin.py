from django.contrib import admin
from django.contrib import messages

from .forms import ToolFieldWhiteListForm
from .forms import ToolForm
from .models import Tool, ToolInputData, ToolOutputData, ToolFlag, ToolFieldWhiteList


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


class ToolOutputDataInline(admin.TabularInline):
    model = ToolOutputData
    extra = 0


class ToolInputDataInline(admin.TabularInline):
    model = ToolInputData
    extra = 0


class ToolFlagInline(admin.TabularInline):
    model = ToolFlag.tool.through
    extra = 0


class ToolAdmin(admin.ModelAdmin):
    """

    """

    list_display = ['name', 'galaxy_server', 'toolshed', 'version', 'visible', 'toolflags']
    list_filter = ['galaxy_server', 'toolshed']
    fields = ('name', 'version', 'oneclick', 'galaxy_server', 'toolshed', 'id_galaxy',)
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


class ToolInputDataAdmin(admin.ModelAdmin):
    pass


class ToolOutputDataAdmin(admin.ModelAdmin):
    filter_horizontal = ('compatible_inputs',)


class ToolInputOutputLinkAdmin(admin.ModelAdmin):
    list_display = ['pk', 'tooloutputdata', 'toolinputdata']


class ToolFlagAdmin(admin.ModelAdmin):
    filter_horizontal = ('tool',)


class ToolFieldWhiteListAdmin(admin.ModelAdmin):
    """
    Custom admin add view to previous tool form fields according to whitelist
    """
    change_form_template = "tools/admin/custom_change_form.html"
    list_display = ['tool', 'context']
    form = ToolFieldWhiteListForm
    fields = ['tool', 'context', '_params', ]

    def change_view(self, request, object_id, form_url='', extra_context=None):
        obj = self.get_object(request, object_id)
        extra_context = extra_context or {}
        extra_context['preview_toolform'] = ToolForm(obj.get_json_params, obj.tool.id_galaxy, whitelist=obj._params)
        extra_context['tool'] = obj.tool.name
        return super(ToolFieldWhiteListAdmin, self).change_view(request, object_id, form_url,
                                                                extra_context=extra_context, )

    def get_form(self, request, obj=None, **kwargs):
        params_choices = [('', '--All--')]
        if obj:
            obj._params = obj.saved_params
            params_choices += list(zip(obj.get_params, obj.get_params))

        form = super(ToolFieldWhiteListAdmin, self).get_form(request, obj, **kwargs)
        form.base_fields['_params'].choices = params_choices

        return form


admin.site.register(Tool, ToolAdmin)
admin.site.register(ToolInputData)
admin.site.register(ToolOutputData, ToolOutputDataAdmin)
admin.site.register(ToolFlag, ToolFlagAdmin)
admin.site.register(ToolOutputData.compatible_inputs.through, ToolInputOutputLinkAdmin)
admin.site.register(ToolFieldWhiteList, ToolFieldWhiteListAdmin)
