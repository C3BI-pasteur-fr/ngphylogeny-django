from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from django.conf import settings


class GalaxyUser(models.Model):
    """
    Model to save User APIkey from Galaxy
    """

    user = models.OneToOneField(User, primary_key=True)
    api_key = models.CharField(max_length=100, blank=True)

    def __unicode__(self):
        return self.user.username


