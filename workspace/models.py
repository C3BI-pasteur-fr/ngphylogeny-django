from __future__ import unicode_literals

import json

from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone
from django.utils.functional import cached_property

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

    # workspace.tasks.deleteoldgalaxyhistory's own daily cleanup cutoff -
    # defined here (not just as a local constant in tasks.py) so the
    # "days until deletion" estimate shown on the workspace/account pages
    # (days_until_deletion below) can't drift out of sync with the value
    # that cleanup task actually uses. Reads settings.
    # NGPHYLO_WORKSPACE_RETENTION_DAYS (see settings/base.py) rather than
    # a bare literal, so a deployment can configure this via the
    # WORKSPACE_RETENTION_DAYS GitLab CI/CD variable / NGPHYLO_WORKSPACE_
    # RETENTION_DAYS env var without a code change - 14 if unset, this
    # project's original hardcoded value. Safe to read at class-body
    # (i.e. import) time: Django only imports app models after settings
    # are fully configured, same as GalaxyUser.GALAXY_REQUEST_TIMEOUT's
    # own class-level constant elsewhere in this codebase.
    RETENTION_DAYS = settings.NGPHYLO_WORKSPACE_RETENTION_DAYS

    @property
    def days_until_deletion(self):
        """
        Days left before workspace.tasks.deleteoldgalaxyhistory's daily
        2am cleanup removes this history - same RETENTION_DAYS cutoff,
        same created_date reference point that task itself filters on.
        That task only actually acts on finished=True rows, but every
        history gets there eventually (either normally, or forced via
        the 24h stale-job cancellation - see that task's own docstring),
        and the cutoff is always measured from created_date regardless,
        so this stays a meaningful estimate even for a still-running
        history. Clamped to 0 rather than going negative once past the
        cutoff (deletion is a daily batch job, not instantaneous - a
        history can briefly sit at "0 days left" for up to a day before
        that job actually runs, or longer still if a transient Galaxy
        failure left it deleted=False for a retry - see that task's own
        notes on this).
        """
        days_left = self.RETENTION_DAYS - (timezone.now() - self.created_date).days
        return max(0, days_left)

    @property
    def type_label(self):
        """
        Human-readable workflow_category label ('duplicated' ->
        'Advanced', etc.) for the Workspace/account history tables' own
        "Type" column - the exact same mapping workspace.views.
        running_jobs_view and workspace.reports already use
        (workspace.reports.CATEGORY_LABELS, via its own _category_label
        helper). Local import, not a module-level one: workspace.reports
        itself imports WorkspaceHistory, so importing reports at this
        module's top level would be circular.
        """
        from .reports import _category_label
        return _category_label(self.workflow_category)

    @cached_property
    def _history_content_steps(self):
        """
        The dict-shaped entries of history_content_json - same tolerant
        parsing convention as running_jobs_view/WorkspaceHistoryObjectMixin
        (a malformed/non-list value degrades to "0 steps" rather than
        crashing whatever's rendering it). Cached per-instance since both
        steps_done and steps_total below read this on every access -
        rendering a table row needs both, not just one.
        """
        try:
            content = json.loads(self.history_content_json or '[]')
        except (ValueError, TypeError):
            content = []
        return [f for f in content if isinstance(f, dict)]

    @property
    def steps_done(self):
        return sum(
            1 for f in self._history_content_steps
            if 'ok' in (f.get('state') or ''))

    @property
    def steps_total(self):
        return len(self._history_content_steps)

    def get_galaxy_user(self):
        """
        Every authenticated visitor now authenticates to Galaxy through
        the same single, shared key regardless of which NGPhylogeny
        account owns this history (see galaxy.decorator.
        connection_galaxy) - nothing creates a personal, per-user
        GalaxyUser row anymore, so an authenticated owner's own row
        almost never exists. Falls back to the shared anonymous
        GalaxyUser instead of raising GalaxyUser.DoesNotExist (a plain
        .get() used to do exactly that) - rename() below calls this on
        every single save(), so that used to mean every save() of a
        history owned by a real account would crash outright once this
        fell back to the shared key elsewhere. A personal row is still
        preferred first if one happens to exist (e.g. from before this
        change), for backwards compatibility.

        Deliberately still returns None (no fallback at all) for a
        history with no self.user - i.e. a session-only, not-logged-in
        submission - preserving this method's pre-existing behavior for
        that case exactly (the previous version had no else branch
        here at all, so rename() below already silently no-ops for
        every anonymous-visitor history and always has; that's
        unrelated to the shared-key change and not something to alter
        as a side effect of it).

        Only prefers a personal row when it actually has an api_key -
        caught live against the real local dev DB: it still had two
        personal GalaxyUser rows with a blank api_key, artifacts of the
        old (now-removed) connection_galaxy code path that used to
        get_or_create() one for every authenticated visitor regardless
        of whether they'd ever set a key. Without this check, `if gu:`
        alone treats that empty-key row as "found" and returns it
        straight away, never reaching the shared-key fallback below -
        exactly the ValueError('API key must be set') this whole
        change exists to prevent, just for a different reason (a stale
        row instead of no row at all).
        """
        if not self.user:
            return None
        gu = GalaxyUser.objects.filter(
            user=self.user, galaxy_server=self.galaxy_server).first()
        if gu and gu.api_key:
            return gu
        return GalaxyUser.objects.filter(
            anonymous=True, galaxy_server=self.galaxy_server).first()

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
            # workspace/reports.py's gather_all_time() groups every row
            # (all-time, by definition) by exactly these three columns -
            # created_date's own index above doesn't help this one at
            # all, since this query never touches created_date. Real
            # production data measured this at ~2.9-3.2s (EXPLAIN
            # ANALYZE, 2026-09-23), a plain Seq Scan reading the whole
            # (wide - history_content_json/history_info_json are TEXT)
            # heap for every one of 701,970 rows just to project out 3
            # narrow columns. A covering index lets Postgres satisfy
            # this via an Index Only Scan instead - the same fix already
            # measured cutting gather_period_totals() from ~3.7s to
            # ~1.3s for the same reason (its own created_date index).
            # Not partial (unlike wsph_running_jobs_idx above): this
            # query genuinely needs every row, all-time.
            models.Index(
                fields=['workflow_category', 'workflow_steps', 'workflow'],
                name='wsph_alltime_category_idx',
            ),
        ]
