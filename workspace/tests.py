from datetime import timedelta
from unittest.mock import Mock, patch

from django.test import TestCase
from django.utils import timezone

from galaxy.models import Server
from workflows.models import Workflow
from workspace.models import WorkspaceHistory
from workspace.tasks import deleteoldgalaxyhistory


class DeleteOldGalaxyHistoryTest(TestCase):
    """
    Regression tests for deleteoldgalaxyhistory(): it used to mark a
    WorkspaceHistory (and its associated Workflow) as deleted=True
    unconditionally, even when the actual Galaxy-side
    delete_history()/delete_workflow() call silently failed (both
    deletegalaxyhistory() and deletegalaxyworkflow() swallow their own
    exceptions and just log a warning). Since every future run only
    looks at deleted=False rows, a transient Galaxy failure at cleanup
    time made the row look "cleaned up" forever - the data could still
    exist on Galaxy with no way for this task to ever notice or retry.
    Only caught by code review, not by any existing test.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    def _make_old_history(self, workflow=None):
        h = WorkspaceHistory.objects.create(
            history='hist1', name='test', email='', monitored=True,
            finished=True, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server, workflow=workflow,
            history_content_json='{"some": "content"}',
            history_info_json='{"more": "content"}')
        # created_date is auto_now_add - .save()/.create() always stamp it
        # with "now", so backdate it via a queryset update (bypasses
        # auto_now_add) to simulate an old, finished history.
        WorkspaceHistory.objects.filter(pk=h.pk).update(
            created_date=timezone.now() - timedelta(days=15))
        h.refresh_from_db()
        return h

    def test_marks_deleted_and_clears_json_on_success(self):
        h = self._make_old_history()
        with patch('workspace.tasks.deletegalaxyhistory', return_value=True):
            deleteoldgalaxyhistory()
        h.refresh_from_db()
        self.assertTrue(h.deleted)
        self.assertEqual(h.history_content_json, "")
        self.assertEqual(h.history_info_json, "")

    def test_leaves_deleted_false_when_galaxy_history_delete_fails(self):
        h = self._make_old_history()
        with patch('workspace.tasks.deletegalaxyhistory', return_value=False):
            deleteoldgalaxyhistory()
        h.refresh_from_db()
        self.assertFalse(h.deleted)
        # Not touched either - still there to retry against next run.
        self.assertEqual(h.history_content_json, '{"some": "content"}')

    def test_leaves_deleted_false_when_associated_workflow_delete_fails(self):
        wf = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='wfid1',
            name='PhyML OneClick', category='duplicated',
            description='PhyML OneClick',
            slug='wfid1_PhyML OneClick_copy')
        h = self._make_old_history(workflow=wf)
        with patch('workspace.tasks.deletegalaxyworkflow',
                   return_value=False), \
             patch('workspace.tasks.deletegalaxyhistory',
                   return_value=True):
            deleteoldgalaxyhistory()
        h.refresh_from_db()
        wf.refresh_from_db()
        self.assertFalse(h.deleted)
        self.assertFalse(wf.deleted)

    def test_one_failing_history_does_not_block_the_others(self):
        h_fail = self._make_old_history()
        h_fail.history = 'hist-fail'
        h_fail.save()
        h_ok = self._make_old_history()
        h_ok.history = 'hist-ok'
        h_ok.save()
        # Re-backdate: .save() above didn't touch created_date (not
        # auto_now), so both are still old - just confirm the setup.
        h_fail.refresh_from_db()
        h_ok.refresh_from_db()

        def fake_delete_history(historyid):
            return historyid != 'hist-fail'

        with patch('workspace.tasks.deletegalaxyhistory',
                   side_effect=fake_delete_history):
            deleteoldgalaxyhistory()

        h_fail.refresh_from_db()
        h_ok.refresh_from_db()
        self.assertFalse(h_fail.deleted)
        self.assertTrue(h_ok.deleted)
