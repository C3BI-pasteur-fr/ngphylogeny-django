from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from django.conf import settings


class GalaxyServer(models.Model):
    """
    Galaxy Server Informations
    """

    url = models.URLField(max_length=254, unique=True, null=False)
    synopsis = models.CharField(max_length=254)
    description = models.TextField(max_length=500)
    contact = models.CharField(max_length=254)
    version = models.CharField(max_length=9)

    def __str__(self):

        return self.url


class GalaxyConf(models.Model):
    """
    Configuration of Galaxy
    """

    active = models.BooleanField(default=False, verbose_name="Active Configuration?")
    galaxy_anonymous_api_access = models.BooleanField(default=False, verbose_name="Anonymous api is activated?")
    galaxy_server = models.ForeignKey(GalaxyServer)
    global_user_name = models.CharField(max_length=50, verbose_name="User used to connect Anonymous Users")
    global_api_key = models.CharField(max_length=50, verbose_name="Anonymous User shared api key")

    def save(self, *args, ** kwargs):
        if self.active is True:
            config = self.objects.all()
            config.update(active=False)
        if self.galaxy_anonymous_api_access is True:
            if not (self.global_user_name and self.global_api_key):
                return

        super(GalaxyConf, self).save(*args, **kwargs)


class GalaxyUser(models.Model):
    """
    Model to save User APIkey from Galaxy
    """

    galaxy_server = models.ForeignKey(GalaxyServer)
    user = models.OneToOneField(User, primary_key=True)
    api_key = models.CharField(max_length=100, blank=True)

    def __unicode__(self):
        return self.user.username


