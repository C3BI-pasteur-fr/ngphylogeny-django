from django.contrib import admin
from models import Workflow
# Register your models here.


class WorkflowsAdmin(admin.ModelAdmin):
    prepopulated_fields = {"slug": ("name",)}


admin.site.register(Workflow,WorkflowsAdmin)