from __future__ import unicode_literals

import requests
from django.contrib.auth.models import User
from django.db import models
from django.utils.functional import cached_property

from .galaxylib import GalaxyInstance


class Server(models.Model):
    """
    Galaxy Server Information
    """
    name = models.CharField(max_length=80)
    url = models.URLField(max_length=254, unique=True, null=False)
    synopsis = models.CharField(max_length=254)
    description = models.TextField(max_length=500)
    contact = models.CharField(max_length=254)
    version = models.CharField(max_length=9, blank=True)
    current = models.BooleanField(default=False,
                                  help_text="Galaxy server currently used to run jobs. "
                                            "Only one server can be used at the same time"
                                  )

    def save(self, *args, **kwargs):

        if self._state.adding:
            connection = requests.get(self.url + '/api/version')
            if connection.status_code == 200:
                self.version = connection.json().get('version_major')

        # only one server can be used at the same time
        if self.current is True:
            Server.objects.all().update(current=False)

        super(Server, self).save(*args, **kwargs)

    def __unicode__(self):
        return "%s %s" % (self.name, self.url)


class GalaxyUser(models.Model):
    """
        Model to save User APIkey associated with a Galaxy server
    """
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE)
    api_key = models.CharField(max_length=100, blank=True)
    anonymous = models.BooleanField(default=False,
                                    verbose_name="Share this user api key with all anonymous users ",
                                    help_text="Important: User api key will be used to make all Galaxy requests "
                                              "for NGPhylogeny Anonymous User"
                                    )

    def save(self, *args, **kwargs):
        # only one Anonymous user can be used for one Galaxy server the same time
        if self.anonymous is True:
            GalaxyUser.objects.filter(galaxy_server=self.galaxy_server).update(anonymous=False)
        super(GalaxyUser, self).save(*args, **kwargs)

    @cached_property
    def get_galaxy_instance(self):
        """
            :return bioblend Galaxy instance object
        """
        if self.api_key:
            return GalaxyInstance(url=self.galaxy_server.url, key=self.api_key)
        else:
            raise ValueError("API key must be set")

    def __unicode__(self):
        return "%s %s" % (self.user.username, self.api_key)

    class Meta:
        unique_together = ['user', 'galaxy_server']
