from django.contrib import admin

# Register your models here.
from .models import *


admin.site.register(GalaxyServer)
admin.site.register(GalaxyConf)
admin.site.register(GalaxyUser)
