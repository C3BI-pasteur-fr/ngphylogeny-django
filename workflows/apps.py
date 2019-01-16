from __future__ import unicode_literals

from django.apps import AppConfig


class WorkflowsConfig(AppConfig):
    name = 'workflows'

    def ready(self):
        import workspace.signals
