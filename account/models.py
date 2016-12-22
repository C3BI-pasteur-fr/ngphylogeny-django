from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from django.conf import settings

from django.core.exceptions import ValidationError

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

    galaxy_server = models.ForeignKey(GalaxyServer)
    galaxy_anonymous_api_access = models.BooleanField(default=False, verbose_name="Enable Anonymous api?", help_text="If True")
    global_user_name = models.CharField(max_length=50, null=True, blank=True, verbose_name="User used to connect Anonymous Users")
    global_api_key = models.CharField(max_length=50, null=True, blank=True, verbose_name="Anonymous User shared api key")
    active = models.BooleanField(default=False, verbose_name="Active Configuration?")

    def clean(self,*args, **kwargs):

        if self.galaxy_anonymous_api_access == True and ((not self.global_api_key) or (not self.global_user_name)):
            raise ValidationError({'global_api_key': 'This field is required.',
                                   'global_user_name': 'This field is required.'
                                   })
        super(GalaxyConf,self).clean(*args, **kwargs)

    def save(self, *args, ** kwargs):
        if self.active is True:
            deactivated_config = GalaxyConf.objects.all()
            deactivated_config.update(active=False)
        super(GalaxyConf, self).save(*args, **kwargs)


class GalaxyUser(models.Model):
    """
    Model to save User APIkey from Galaxy
    """

    user = models.OneToOneField(User, primary_key=True)
    galaxy_server = models.ForeignKey(GalaxyServer)
    api_key = models.CharField(max_length=100, blank=True)

    def __unicode__(self):
        return self.user.username


