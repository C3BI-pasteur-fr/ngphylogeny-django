from django.contrib import admin
from django.contrib import messages
# Register your models here.
from .models import *

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

admin.site.register(GalaxyServer)
admin.site.register(GalaxyConf, GalaxyConfAdmin)
admin.site.register(GalaxyUser)
