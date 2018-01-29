import requests
from django.contrib import admin
from django.contrib import messages
from django.core.validators import ValidationError

from tools.models import Tool
from .models import Server, GalaxyUser


class GalaxyServerAdmin(admin.ModelAdmin):
    """

    """
    list_display = ('name', 'url', 'current')
    list_display_links = ('name', 'url')
    actions = ['activate_configuration', 'import_new_tools']

    def save_model(self, request, obj, form, change):
        # get automatically Galaxy version if it not set
        if not obj.version:
            try:
                r = requests.get(obj.url + '/api/version')
                obj.version = r.json().get('version_major')
            except Exception as e:
                raise ValidationError(message="Please set the version of Galaxy, or make sure that URL is working")

        super(GalaxyServerAdmin, self).save_model(request, obj, form, change)

    def import_new_tools(self, request, queryset):
        """
        Import tools description from Galaxy server
        """
        for galaxy_server in queryset:
            tool_import_report = Tool.import_tools(galaxy_server)

            if tool_import_report['new']:
                self.message_user(request, "%s successfully imported new tools." % (len(tool_import_report['new'])))
            else:
                self.message_user(request, "No new tools imported.",
                                  level=messages.WARNING
                                  )

    def activate_configuration(self, request, queryset):
        """
        Activate the default Galaxy Server to run analyses
        """
        if len(queryset) > 1:
            self.message_user(request, "Please select only one available configuration", level=messages.WARNING)

        elif queryset:
            galaxy_server = queryset[0]
            galaxy_server.current = True
            galaxy_server.save()

            self.message_user(request,
                              "%s configuration is now activated and usable by the application" % (galaxy_server.name),
                              level=messages.SUCCESS)


class GalaxyUserAdmin(admin.ModelAdmin):
    """
    """
    list_display = ('user', 'galaxy_server', 'anonymous')
    list_filter = ('user', 'galaxy_server', 'anonymous')


admin.site.register(Server, GalaxyServerAdmin)
admin.site.register(GalaxyUser, GalaxyUserAdmin)
