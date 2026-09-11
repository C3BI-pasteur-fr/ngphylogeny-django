from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase

from galaxy.models import GalaxyUser, Server
from tools.models import Tool
from workflows.models import Workflow, WorkflowStepInformation
from workflows.views.wkadvanced import WorkflowAdvancedFormView


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
