from django.contrib import admin

from .models import WorkspaceHistory


class WorkspaceAdmin(admin.ModelAdmin):

    list_display = ('name','history','user','created_date', 'galaxy_server',)

admin.site.register(WorkspaceHistory,WorkspaceAdmin)
