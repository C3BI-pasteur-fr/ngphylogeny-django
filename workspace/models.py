from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import pre_delete

from account.models import GalaxyServer, GalaxyUser


class WorkspaceHistory(models.Model):
    """
    Galaxy history information
    """
    history = models.CharField(max_length=20,unique=True)
    name = models.CharField(max_length=100)
    created_date = models.DateField(auto_now=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    galaxy_server = models.ForeignKey(GalaxyServer, on_delete=models.CASCADE )

    def get_galaxy_user(self):
        if self.user:
            return GalaxyUser.objects.get(user=self.user, galaxy_server=self.galaxy_server)

    def rename(self):
        """Rename history galaxy"""
        gu = self.get_galaxy_user()
        if gu:
            gi = gu.get_galaxy_instance()
            gi.histories.update_history(history_id=self.history, name=self.name)

    def save(self, *args, **kwargs):
        if self.name:
            self.rename()
        super(WorkspaceHistory, self).save(*args, **kwargs)


    class Meta:
        verbose_name_plural = "Workspace histories"
        unique_together = (("history", "galaxy_server"),)


def send_delete_galaxy_history(sender, instance, using, **kwargs):
    """remove history from db and from Galaxy serveur"""

    gu = instance.get_galaxy_user()
    if gu:
        gi = gu.get_galaxy_instance()
        gi.histories.delete_history(instance.history, purge=True)

pre_delete.connect(send_delete_galaxy_history, sender=WorkspaceHistory)