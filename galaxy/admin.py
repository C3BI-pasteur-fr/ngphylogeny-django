from django.contrib import admin
from django.contrib import messages

from tools.models import Tool
from .models import *


class GalaxyServerAdmin(admin.ModelAdmin):
    """

    """
    list_display = ('name', 'url', 'current')
    list_display_links = ('name', 'url')
    actions = ['activate_configuration','import_new_tools' ]

    def import_new_tools(self, request, queryset):
        """
        Import tools description from Galaxy server
        """
        for galaxy_server in queryset:
            tools_list = Tool.import_tools(galaxy_server)
            self.message_user(request, "%s successfully imported new tools." % (len(tools_list)))

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

            self.message_user(request, "%s configuration is now activated and usable by the application" %(galaxy_server.name),
                              level=messages.SUCCESS )


class GalaxyUserAdmin(admin.ModelAdmin):
    """

    """
    list_display = ('user', 'galaxy_server', 'anonymous')
    list_filter = ('user', 'galaxy_server', 'anonymous')


admin.site.register(Server, GalaxyServerAdmin)
admin.site.register(GalaxyUser, GalaxyUserAdmin)
