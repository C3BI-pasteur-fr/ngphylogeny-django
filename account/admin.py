from django.contrib import admin
from django.contrib import messages
# Register your models here.
from .models import *
from tools.models import Tool

class GalaxyConfAdmin(admin.ModelAdmin):

    list_display = ('name', 'galaxy_server', 'active')
    actions = ['activate_configuration']

    def activate_configuration(self, request, queryset):

        if len(queryset)>1:
            self.message_user(request, "Please select only one available configuration",level=messages.WARNING)
        elif queryset:
            galaxy_conf = queryset[0]
            galaxy_conf.active=True
            galaxy_conf.save()

            self.message_user(request, "%s configuration is now activated" %(galaxy_conf.name))


class GalaxyServerAdmin(admin.ModelAdmin):

    list_display = ('name', 'url')
    list_display_links = ('name', 'url')
    actions = ['import_new_tools']

    def import_new_tools(self, request, queryset):
        """
        Import tools description from Galaxy server
        """
        for galaxy_server in queryset:

            tools_list = Tool.import_tools(galaxy_server)
            self.message_user(request, "%s successfully imported new tools." %(len(tools_list)))


class GalaxyUserAdmin(admin.ModelAdmin):

    list_display = ('user','galaxy_server')
    list_filter = ('user','galaxy_server')

admin.site.register(GalaxyServer, GalaxyServerAdmin)
admin.site.register(GalaxyConf, GalaxyConfAdmin)
admin.site.register(GalaxyUser, GalaxyUserAdmin )
