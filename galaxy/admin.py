from django.contrib import admin
from django.contrib import messages

from tools.models import Tool
from .models import *


class GalaxyServerAdmin(admin.ModelAdmin):

    list_display = ('name', 'url', 'current' )
    list_display_links = ('name', 'url')
    actions = ['import_new_tools','activate_configuration']

    def import_new_tools(self, request, queryset):
        """
        Import tools description from Galaxy server
        """
        for galaxy_server in queryset:

            tools_list = Tool.import_tools(galaxy_server)
            self.message_user(request, "%s successfully imported new tools." %(len(tools_list)))

    def activate_configuration(self, request, queryset):

        if len(queryset)>1:
            self.message_user(request, "Please select only one available configuration",level=messages.WARNING)
        elif queryset:
            galaxy_conf = queryset[0]
            galaxy_conf.active=True
            galaxy_conf.save()

            self.message_user(request, "%s configuration is now activated" %(galaxy_conf.name))





class GalaxyUserAdmin(admin.ModelAdmin):

    list_display = ('user','galaxy_server')
    list_filter = ('user','galaxy_server')

admin.site.register(Server, GalaxyServerAdmin)
admin.site.register(GalaxyUser, GalaxyUserAdmin)
