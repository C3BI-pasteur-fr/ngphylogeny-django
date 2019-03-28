from __future__ import unicode_literals

from django.contrib.auth.models import User
from django.db import models

from galaxy.models import Server, GalaxyUser
from workflows.models import Workflow


class WorkspaceHistory(models.Model):
    """
    Galaxy history information
    """
    # uuid = models.UUIDField(unique=True, default=uuid.uuid4, editable=False)
    history = models.CharField(max_length=20)
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE)
    # The potential workflow that has been executed in the workspace
    workflow = models.ForeignKey(Workflow, null=True)
    name = models.CharField(max_length=100)
    created_date = models.DateTimeField(auto_now=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    email = models.CharField(max_length=100)
    monitored = models.BooleanField(default=False)
    finished = models.BooleanField(default=False)
    source_ip=models.CharField(max_length=20, default="")
    # Where the workflow comes from : Oneclick, Advanced or ALaCarte
    # OneClick|ALaCarte|SingleTool
    workflow_category = models.CharField(max_length=100, default="")
    # Workflow Steps PhyML, etc.
    workflow_steps = models.CharField(max_length=100, default="")
    # If the galaxy history has been deleted
    # We keep the workspace on the django side but
    # describe it as deleted (it won't appear anymore
    # on workspace history)
    deleted = models.BooleanField(default=False)
    # history Json coming from galaxy server: stored in the database
    history_content_json = models.TextField(default="{}")
    history_info_json = models.TextField(default="{}")
    # Deserialized json, not stored in the django database
    history_content = None
    history_info = None

    def get_galaxy_user(self):
        if self.user:
            return GalaxyUser.objects.get(user=self.user,
                                          galaxy_server=self.galaxy_server)

    def rename(self):
        """Rename history galaxy"""
        gu = self.get_galaxy_user()
        if gu:
            gi = gu.get_galaxy_instance
            gi.histories.update_history(history_id=self.history,
                                        name=self.name)

    def save(self, *args, **kwargs):
        if self.name:
            self.rename()
        super(WorkspaceHistory, self).save(*args, **kwargs)

    class Meta:
        verbose_name_plural = "Workspace histories"
        unique_together = (("history", "galaxy_server"),)
