from __future__ import unicode_literals

from django.apps import AppConfig


class GalaxyConfig(AppConfig):
    name = 'galaxy'

    def ready(self):
        import signals
