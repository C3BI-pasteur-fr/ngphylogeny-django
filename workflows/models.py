from __future__ import unicode_literals

from django.db import models

from galaxy.models import Server


class Workflow(models.Model):
    """
    Galaxy Workflow informations
    """
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE, null=True, blank=True)
    id_galaxy = models.CharField(max_length=250, unique=True)
    name = models.CharField(max_length=100, unique=True)
    category = models.CharField(max_length=100 , blank=True)
    version = models.CharField(max_length=10, blank=True)
    description = models.CharField(max_length=250)
    slug = models.SlugField(max_length=100,unique=True)

