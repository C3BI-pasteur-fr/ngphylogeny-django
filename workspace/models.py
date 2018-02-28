from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models

from galaxy.models import Server, GalaxyUser

import uuid
class WorkspaceHistory(models.Model):
    """
    Galaxy history information
    """
    # uuid = models.UUIDField(unique=True, default=uuid.uuid4, editable=False)
    history = models.CharField(max_length=20)
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE)
    name = models.CharField(max_length=100)
    created_date = models.DateTimeField(auto_now=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)

    def get_galaxy_user(self):
        if self.user:
            return GalaxyUser.objects.get(user=self.user, galaxy_server=self.galaxy_server)

    def rename(self):
        """Rename history galaxy"""
        gu = self.get_galaxy_user()
        if gu:
            gi = gu.get_galaxy_instance
            gi.histories.update_history(history_id=self.history, name=self.name)

    def save(self, *args, **kwargs):
        if self.name:
            self.rename()
        super(WorkspaceHistory, self).save(*args, **kwargs)

    class Meta:
        verbose_name_plural = "Workspace histories"
        unique_together = (("history", "galaxy_server"),)


