from datetime import timedelta
from unittest.mock import Mock, patch

from bioblend.galaxy.client import ConnectionError
from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from galaxy.models import GalaxyUser, Server
from tools.models import Tool
from workflows.models import Workflow, WorkflowStepInformation
from workflows.tasks import deleteoldgalaxyworkflows
from workflows.views.wkadvanced import WorkflowAdvancedFormView
from workspace.models import WorkspaceHistory


class WorkflowStepInformationTest(TestCase):
    """
    Regression test: WorkflowStepInformation.update_dict_tools() drops
    steps whose tool isn't locally known by deleting from
    self.steps_tooldict while iterating it. That relied on Python 2's
    dict.items() returning a list snapshot; under Python 3 (a live
    view), it raised "RuntimeError: dictionary changed size during
    iteration" for any real Galaxy workflow containing a step this
    app's Tool table doesn't have a matching row for - only caught by
    testing against real (restored production) data, not by any of
    the existing unit tests, which never exercised the "unknown tool"
    branch.
    """

    def test_step_with_unknown_tool_is_dropped_without_error(self):
        workflow_json = {
            "steps": {
                "0": {"tool_id": None},
                "1": {
                    "tool_id": "toolshed.example.org/repos/x/y/unknown_tool/1.0",
                    "annotation": "",
                    "tool_inputs": {},
                },
            }
        }
        # No Tool rows exist in the test DB, so the step's tool is
        # necessarily unknown - exercising the deletion branch.
        self.assertEqual(Tool.objects.count(), 0)
        info = WorkflowStepInformation(workflow_json)
        self.assertEqual(info.steps_tooldict, {})
        self.assertEqual(info.sorted_tool_list, [])


class ImportWorkflowsCommandTest(TestCase):
    """
    Regression test: import_workflows() used to build a new Workflow row
    for every matching Galaxy workflow on every run, with no dedup. The
    Galaxy image used in deployment re-imports its bundled workflows (with
    a fresh id_galaxy each time) on every restart, so re-running this
    command against the same server crashed with "duplicate key value
    violates unique constraint workflows_workflow_slug_key" instead of
    refreshing the existing row - only caught by running the actual
    docker-compose deployment's init step twice against a real Galaxy
    server.
    """

    def setUp(self):
        user = User.objects.create_user('admin')
        # Server.save() itself pings <url>/api/version on creation.
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)

    @staticmethod
    def _mock_workflow_list(id_galaxy):
        return Mock(status_code=200, json=lambda: [
            {'id': id_galaxy, 'name': 'PhyML OneClick'},
        ])

    def test_rerunning_against_the_same_server_updates_not_duplicates(self):
        with patch('tools.management.commands.importworkflows.requests.get',
                   return_value=self._mock_workflow_list('galaxyid1')):
            call_command('importworkflows', galaxyurl=self.server.url,
                         wfnamefile='wfnames.txt')

        self.assertEqual(Workflow.objects.count(), 1)
        self.assertEqual(Workflow.objects.get().id_galaxy, 'galaxyid1')

        # Simulate Galaxy re-importing its bundled workflow on restart:
        # same name, fresh id_galaxy.
        with patch('tools.management.commands.importworkflows.requests.get',
                   return_value=self._mock_workflow_list('galaxyid2')):
            call_command('importworkflows', galaxyurl=self.server.url,
                         wfnamefile='wfnames.txt')

        self.assertEqual(Workflow.objects.count(), 1)
        self.assertEqual(Workflow.objects.get().id_galaxy, 'galaxyid2')

    def test_wfids_imports_directly_by_id(self):
        """
        --wfids fetches specific known workflow ids directly
        (/api/workflows/{id}, one GET per id) instead of listing/filtering
        the server's whole workflow collection - precise and fast once
        the target Galaxy's stable base-workflow ids are already known
        (e.g. from a database dump), unlike the oldest-N heuristic
        import_workflows() uses.
        """
        mock_response = Mock(status_code=200, json=lambda: {
            'id': 'known-id', 'name': 'FastME OneClick'})
        with patch('tools.management.commands.importworkflows.requests.get',
                   return_value=mock_response):
            call_command('importworkflows', galaxyurl=self.server.url,
                         wfids='known-id')

        self.assertEqual(Workflow.objects.count(), 1)
        wf = Workflow.objects.get()
        self.assertEqual(wf.id_galaxy, 'known-id')
        self.assertEqual(wf.category, 'base')

    def test_ignores_per_run_duplicates_sharing_the_same_name(self):
        """
        Regression test: Workflow.duplicate() gives every real user run its
        own Galaxy-side copy of the base workflow (same name, fresh
        Galaxy id, tracked locally as its own category='duplicated' row -
        see workflows/models.py). Galaxy's /api/workflows/ list returns
        those alongside the one true base workflow, all sharing the exact
        same name and therefore the exact same slug. Re-running
        importworkflows used to try to collapse every one of them into
        the single 'base' row via update_or_create(slug=...), which
        crashed with "duplicate key value violates unique constraint
        workflows_workflow_id_galaxy_key" as soon as one of those already
        belonged to an existing 'duplicated' row - only caught by
        redeploying against a real Galaxy server that had accumulated
        real workflow runs.
        """
        base = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='base-id',
            name='PhyML OneClick', category='base',
            description='PhyML OneClick', slug='phyml-oneclick')
        Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='dup-id',
            name='PhyML OneClick', category='duplicated',
            description='PhyML OneClick',
            slug='dup-id_PhyML OneClick_copy')

        mock_response = Mock(status_code=200, json=lambda: [
            {'id': 'dup-id', 'name': 'PhyML OneClick'},
            {'id': 'base-id', 'name': 'PhyML OneClick'},
        ])
        with patch('tools.management.commands.importworkflows.requests.get',
                   return_value=mock_response):
            call_command('importworkflows', galaxyurl=self.server.url,
                         wfnamefile='wfnames.txt')

        self.assertEqual(Workflow.objects.count(), 2)
        base.refresh_from_db()
        self.assertEqual(base.id_galaxy, 'base-id')


class DeleteOldGalaxyWorkflowsTest(TestCase):
    """
    Regression tests for deleteoldgalaxyworkflows(): it used to call
    Galaxy's delete_workflow() directly and unconditionally follow up with
    w.delete() (a hard delete of the Django row) regardless of whether the
    Galaxy call actually succeeded - delete_workflow() failures only
    surfaced as an exception caught by the loop's own try/except, which
    logged a warning and then aborted the *entire* remaining batch for
    that run (rather than skipping just the one problem row), and any
    workflow processed before the failure was already gone from Django
    even on runs where the Galaxy-side delete silently didn't take
    effect. Only caught by code review, not by any existing test.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    def _make_orphan_workflow(self, id_galaxy):
        return Workflow.objects.create(
            galaxy_server=self.server, id_galaxy=id_galaxy,
            name='PhyML OneClick', category='duplicated',
            description='PhyML OneClick',
            slug='%s_PhyML OneClick_copy' % id_galaxy,
            # Older than deleteoldgalaxyworkflows()'s 7-day cutoff.
            date=timezone.now() - timedelta(days=8))

    def test_deletes_orphaned_workflow_on_success(self):
        wf = self._make_orphan_workflow('orphan-ok')
        with patch('workflows.tasks.deletegalaxyworkflow',
                   return_value=True):
            deleteoldgalaxyworkflows()
        self.assertFalse(Workflow.objects.filter(pk=wf.pk).exists())

    def test_keeps_workflow_row_when_galaxy_delete_fails(self):
        wf = self._make_orphan_workflow('orphan-fail')
        with patch('workflows.tasks.deletegalaxyworkflow',
                   return_value=False):
            deleteoldgalaxyworkflows()
        self.assertTrue(Workflow.objects.filter(pk=wf.pk).exists())

    def test_one_failure_does_not_block_the_rest_of_the_batch(self):
        wf_fail = self._make_orphan_workflow('orphan-fail')
        wf_ok = self._make_orphan_workflow('orphan-ok')

        def fake_delete(id_galaxy):
            return id_galaxy != 'orphan-fail'

        with patch('workflows.tasks.deletegalaxyworkflow',
                   side_effect=fake_delete):
            deleteoldgalaxyworkflows()

        self.assertTrue(Workflow.objects.filter(pk=wf_fail.pk).exists())
        self.assertFalse(Workflow.objects.filter(pk=wf_ok.pk).exists())

    def test_workflow_associated_with_a_history_is_also_deleted(self):
        """
        deleteoldgalaxyworkflows() no longer skips a workflow just because
        a WorkspaceHistory references it - real usage left hundreds of
        thousands of old, actually-run duplicated workflows never cleaned
        up by either this task (originally never-run-only) or
        deleteoldgalaxyhistory's own narrower per-history cleanup. A
        workflow's Galaxy definition is safe to drop independently of its
        history's own (much longer) retention window, since the history's
        actual data lives separately from the workflow "recipe" that
        launched it.
        """
        wf = self._make_orphan_workflow('in-use')
        history = WorkspaceHistory.objects.create(
            history='hist1', name='test', email='', monitored=True,
            finished=False, source_ip='127.0.0.1',
            workflow_category='OneClick', workflow_steps='',
            galaxy_server=self.server, workflow=wf)
        with patch('workflows.tasks.deletegalaxyworkflow',
                   return_value=True) as mock_delete:
            deleteoldgalaxyworkflows()
        mock_delete.assert_called_once_with('in-use')
        self.assertFalse(Workflow.objects.filter(pk=wf.pk).exists())
        # w.delete() SET_NULLs any WorkspaceHistory.workflow FK pointing
        # here - deleteoldgalaxyhistory relies on seeing workflow=None to
        # know there's nothing left to delete if it processes this same
        # history afterwards.
        history.refresh_from_db()
        self.assertIsNone(history.workflow)


class ProcessFileToUploadTest(TestCase):
    """
    Regression test: process_file_to_upload()'s non-uploaded-file branch
    (request.POST.get("file") for pasted text, or a BlastRun.to_fasta()
    result) used to crash with "TypeError: a bytes-like object is required,
    not 'str'" - it wrote a plain str straight into a NamedTemporaryFile(),
    which defaults to binary mode. Only caught by actually pasting text
    through the live A La Carte/advanced form (POST /workflows/wkmake/<id>),
    not by any existing test.
    """

    def test_pasted_text_str_input_does_not_crash(self):
        view = WorkflowAdvancedFormView()
        fasta = ">s1\nACGT\n>s2\nACGT\n>s3\nACGT\n>s4\nACGT\n"
        tmp_file, name, nseq, length, seqaa = view.process_file_to_upload(
            fasta, "pasted.fasta")
        self.assertEqual(nseq, 4)
        self.assertEqual(length, 4)


class GalaxyUnavailableTest(TestCase):
    """
    Regression test: a Galaxy outage (e.g. Galaxy's own maintenance
    mode, returning a 503) used to 500 every workflow list/form page -
    WorkflowListView.workflow_list, WorkflowFormView.get_context_data()
    and WorkflowAdvancedFormView.get_context_data() all called
    bioblend's show_workflow() with nothing catching a connection
    failure. Real traceback this reproduces: bioblend.ConnectionError
    ("Site Maintenance") raised straight out of workflows/models.py's
    fetch_details(), from a plain GET /workflows/advanced/. All three
    pages should now render workflows/galaxy_unavailable.html with a
    503 instead of Django's generic 500 page.
    """

    def setUp(self):
        user = User.objects.create_user('admin')
        # Server.save() itself pings <url>/api/version on creation.
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)
        self.workflow = Workflow.objects.create(
            galaxy_server=self.server, id_galaxy='galaxyid1',
            name='PhyML OneClick', category='base', description='d',
            slug='phyml-oneclick')

    @staticmethod
    def _galaxy_unreachable():
        return patch(
            'bioblend.galaxy.workflows.WorkflowClient.show_workflow',
            side_effect=ConnectionError(
                "GET: error 503: b'Site Maintenance'", status_code=503))

    def test_oneclick_list_view_returns_503_not_500(self):
        with self._galaxy_unreachable():
            response = self.client.get(reverse('workflow_oneclick_list'))
        self.assertEqual(response.status_code, 503)
        self.assertTemplateUsed(
            response, 'workflows/galaxy_unavailable.html')

    def test_advanced_list_view_returns_503_not_500(self):
        with self._galaxy_unreachable():
            response = self.client.get(reverse('workflows_advanced'))
        self.assertEqual(response.status_code, 503)
        self.assertTemplateUsed(
            response, 'workflows/galaxy_unavailable.html')

    def test_oneclick_form_view_returns_503_not_500(self):
        with self._galaxy_unreachable():
            response = self.client.get(reverse(
                'workflow_oneclick_form',
                kwargs={'slug': self.workflow.slug}))
        self.assertEqual(response.status_code, 503)
        self.assertTemplateUsed(
            response, 'workflows/galaxy_unavailable.html')

    def test_advanced_form_view_returns_503_not_500(self):
        with self._galaxy_unreachable():
            response = self.client.get(reverse(
                'workflows_advanced_fullsteps',
                kwargs={'slug': self.workflow.slug}))
        self.assertEqual(response.status_code, 503)
        self.assertTemplateUsed(
            response, 'workflows/galaxy_unavailable.html')
