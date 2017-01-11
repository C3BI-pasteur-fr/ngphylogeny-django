from django.contrib import admin

# Register your models here.
from .models import *

class GalaxyConfAdmin(admin.ModelAdmin):

    list_display = ('name', 'galaxy_server', 'active')


admin.site.register(GalaxyServer)
admin.site.register(GalaxyConf, GalaxyConfAdmin)
admin.site.register(GalaxyUser)
