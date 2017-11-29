from django.contrib import admin

from .models import ExampleFile


class ExampleFileAdmin(admin.ModelAdmin):
    """
    """
    list_display = ('name', 'ext',)


admin.site.register(ExampleFile, ExampleFileAdmin)
