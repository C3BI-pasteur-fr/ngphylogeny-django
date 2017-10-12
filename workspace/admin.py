from django.contrib import admin

from .models import WorkspaceHistory


class WorkspaceAdmin(admin.ModelAdmin):
    change_form_template = "workspace/admin/custom_change_form.html"
    list_display = ('name', 'history', 'user', 'created_date', 'galaxy_server',)

    def change_view(self, request, object_id, form_url='', extra_context=None):
        obj = self.get_object(request, object_id)
        extra_context = extra_context or {}
        gu = obj.get_galaxy_user()
        if gu:
            gi = gu.get_galaxy_instance
            gi.nocache = True
            extra_context['history_content'] = gi.histories.show_history(obj.history, contents=True)

        return super(WorkspaceAdmin, self).change_view(request, object_id, form_url,
                                                       extra_context=extra_context, )


admin.site.register(WorkspaceHistory, WorkspaceAdmin)
