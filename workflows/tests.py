from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase

from galaxy.models import GalaxyUser, Server
from tools.models import Tool
from workflows.models import Workflow, WorkflowStepInformation


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
