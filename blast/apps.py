# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.apps import AppConfig


class BlastConfig(AppConfig):
    name = 'blast'
    def ready(self):
        super(BlastConfig,self).ready()
        pass
