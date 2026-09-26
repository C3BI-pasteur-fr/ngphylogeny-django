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

    def __str__(self):
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

    # bioblend's own default (timeout=None) waits forever - a request
    # handler blocked on one Galaxy call ties up one of a fixed, small
    # number of worker slots (uwsgi's own "--processes 4 --threads 2" in
    # the k8s deployment - see manifest.yaml) until Galaxy answers or
    # uwsgi's own 120s harakiri kills that worker outright. Bound well
    # under that so a stalled/degraded Galaxy frees the worker back up
    # quickly instead - real production incident: a galaxy.pasteur.fr
    # 502 spell (workspace.views.get_dataset_toolprovenance et al., now
    # polled every 10s per dataset by the history detail page's step
    # chain/table - see CLAUDE.md) piled up requests on every open
    # history page until enough worker slots were stuck that the
    # liveness probe (a plain GET /status, no special exemption from
    # needing a free worker) itself couldn't get served in time - 3
    # missed probes and Kubernetes restarts the pod (observed: 5
    # restarts, exit code 137/SIGKILL).
    GALAXY_REQUEST_TIMEOUT = 30

    @cached_property
    def get_galaxy_instance(self):
        """
            :return bioblend Galaxy instance object
        """
        if self.api_key:
            # bioblend.galaxy.GalaxyInstance.__init__ (unlike the lower-
            # level GalaxyClient it wraps) doesn't actually accept/forward
            # a timeout kwarg at all - has to be set as a plain attribute
            # afterward instead (make_get_request/make_post_request just
            # read self.timeout off the instance on each call).
            gi = GalaxyInstance(url=self.galaxy_server.url, key=self.api_key)
            gi.timeout = self.GALAXY_REQUEST_TIMEOUT
            return gi
        else:
            raise ValueError("API key must be set")

    def __str__(self):
        return "%s %s" % (self.user.username, self.api_key)

    class Meta:
        unique_together = ['user', 'galaxy_server']
