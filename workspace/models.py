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
    workflow = models.ForeignKey(Workflow, null=True, on_delete=models.SET_NULL)
    name = models.CharField(max_length=100)
    # db_index=True: workspace/reports.py's daily-report queries
    # (gather_last_7_days, gather_period_totals) all filter/order/group
    # by this column - with no index at all, real production data
    # (701,957 rows) measured ~3.7-4s *each* for these (EXPLAIN ANALYZE,
    # 2026-09-23), a full sequential scan/external sort every time,
    # dominating the report's total generation time far more than the
    # matplotlib rendering this file used to (and still does) blame.
    created_date = models.DateTimeField(auto_now_add=True, db_index=True)
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
        indexes = [
            # workspace.views.running_jobs_view's WorkspaceHistory query
            # filters on exactly these three columns - with no index at
            # all, real production data (701,949 rows, only 9 actually
            # matching) measured a ~2.8s parallel sequential scan for
            # this one query alone (EXPLAIN ANALYZE, 2026-09-23),
            # directly explaining that page's slow load. A partial
            # index (only the still-running rows) stays tiny forever -
            # it grows with how many jobs are currently running, not
            # with the whole table - rather than a full index over
            # every historical row.
            models.Index(
                fields=['created_date'],
                condition=models.Q(
                    monitored=True, finished=False, deleted=False),
                name='wsph_running_jobs_idx',
            ),
        ]
