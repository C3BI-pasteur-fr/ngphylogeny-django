from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User

from django.core.exceptions import ValidationError
from galaxylib import GalaxyInstance


class GalaxyServer(models.Model):
    """
    Galaxy Server Informations
    """
    name = models.CharField(max_length=80)
    url = models.URLField(max_length=254, unique=True, null=False)
    synopsis = models.CharField(max_length=254)
    description = models.TextField(max_length=500)
    contact = models.CharField(max_length=254)
    version = models.CharField(max_length=9)

    def __str__(self):

        return "%s %s" %(self.name, self.url)


class GalaxyConf(models.Model):
    """
    Model to setup Galaxy Configuration and manage activation of Galaxyserver
    :galaxy_server: in GalaxyConf is used to Galaxy requests when status active is true.
    """
    name = models.CharField(max_length=80,)
    galaxy_server = models.ForeignKey(GalaxyServer)
    anonymous_user_api_access = models.BooleanField(default=False, verbose_name="Enable Anonymous api?", help_text="If True, Anonymous user must be set")
    anonymous_user = models.OneToOneField('GalaxyUser', blank=True, null=True, verbose_name="Global Anonymous User",
                                          help_text="This User is used to make Galaxy requests for NGPhylogeny Anonymous Users")
    active = models.BooleanField(default=False, verbose_name="Active Configuration?")

    def clean(self,*args, **kwargs):
        """make anonymous_user as required field them anonymous_user_api_access is true"""
        if self.anonymous_user_api_access == True and (not self.anonymous_user):
            raise ValidationError({'anonymous_user': 'This field is required.',
                                   })
        if self.anonymous_user:
            if not (self.anonymous_user.galaxy_server == self.galaxy_server):
                raise ValidationError({'anonymous_user': 'Anonymous user must have the same Galaxy Server',
                                       })
        super(GalaxyConf,self).clean(*args, **kwargs)

    def save(self, *args, ** kwargs):
        if self.active is True:
            """One configuration can be actived"""
            deactivated_config = GalaxyConf.objects.all()
            deactivated_config.update(active=False)
        super(GalaxyConf, self).save(*args, **kwargs)

    def __unicode__(self):
        return self.name


class GalaxyUser(models.Model):
    """
    Model to save User APIkey from Galaxy
    """
    user = models.ForeignKey(User)
    galaxy_server = models.ForeignKey(GalaxyServer)
    api_key = models.CharField(max_length=100, blank=True)

    def get_galaxy_instance(self):
        """:return bioblend Galaxy instance object
        """
        if self.api_key:
            return GalaxyInstance(url=self.galaxy_server.url, key=self.api_key)
        else:
            raise "API key must be set"

    def __unicode__(self):
        return self.user.username

    class Meta:

        unique_together = ['user','galaxy_server']
