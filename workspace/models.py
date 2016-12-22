from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from account.models import GalaxyServer
# Create your models here.


class WorkspaceHistory(models.Model):

    history = models.CharField(max_length=20,unique=True)
    name = models.CharField(max_length=100)
    created_date = models.DateField(auto_now=True)
    user = models.ForeignKey(User, null=True)
    galaxy_server = models.ForeignKey(GalaxyServer)

    class Meta:

        unique_together = (("history", "galaxy_server"),)
