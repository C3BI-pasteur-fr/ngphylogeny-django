from django.contrib import admin

# Register your models here.
# Register your models here.
from .models import *

class WorkspaceAdmin(admin.ModelAdmin):

    list_display = ('user', 'galaxy_server',)

admin.site.register(WorkspaceHistory,WorkspaceAdmin)
