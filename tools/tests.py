from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase
from django.urls import reverse

from galaxy.models import GalaxyUser, Server
from tools.models import Tool, ToolFlag


class ToolCanRunOnDataTest(TestCase):
    """
    Pure boundary-logic tests for Tool.can_run_on_data(), built on an
    unsaved Tool instance so no DB write / Galaxy API call is triggered
    (Tool.save() fetches tool JSON from a live Galaxy server).
    """

    def make_tool(self, **limits):
        return Tool(**limits)

    def test_no_limits_always_allowed(self):
        tool = self.make_tool()
        self.assertTrue(tool.can_run_on_data(nseq=100000, length=100000, nboot=1000, seqaa=False))

    def test_max_nbseq_enforced(self):
        tool = self.make_tool(max_nbseq=10)
        self.assertTrue(tool.can_run_on_data(nseq=10, length=1, nboot=0, seqaa=False))
        self.assertFalse(tool.can_run_on_data(nseq=11, length=1, nboot=0, seqaa=False))

    def test_max_nbseq_scaled_for_amino_acids(self):
        tool = self.make_tool(max_nbseq=10, aa_scale_factor=2)
        self.assertTrue(tool.can_run_on_data(nseq=5, length=1, nboot=0, seqaa=True))
        self.assertFalse(tool.can_run_on_data(nseq=6, length=1, nboot=0, seqaa=True))

    def test_max_boot_enforced(self):
        tool = self.make_tool(max_boot=100)
        self.assertTrue(tool.can_run_on_data(nseq=1, length=1, nboot=100, seqaa=False))
        self.assertFalse(tool.can_run_on_data(nseq=1, length=1, nboot=101, seqaa=False))

    def test_string_representation(self):
        tool = self.make_tool(name="PhyML", version="3.1")
        self.assertEqual(str(tool), "PhyML - 3.1")


class ToolListViewGroupingTest(TestCase):
    """
    Regression test for templates/tools/tool_list.html's
    {% regroup tool_list|dictsort:"first_flag.verbose_name" by
    first_flag %}: it used to read toolflag_set.first.verbose_name
    directly, relying on dictsort auto-calling .first(). Django 3.1
    hardened dictsort's variable resolver to never call methods (to
    stop sort keys from triggering side-effecting methods), so that
    stopped working - dictsort silently returned "" and the tools
    page rendered as if there were zero tools, for any tool that
    actually has a flag (i.e. always, in real data). Only caught by
    testing against real (restored production) data; Tool.save() and
    Server.save() both make live Galaxy HTTP calls, which is why this
    uses bulk_create() to build fixtures without touching either.
    """

    def setUp(self):
        server = Server(name="Test Galaxy", url="http://example.invalid",
                        current=True)
        Server.objects.bulk_create([server])
        self.server = Server.objects.get(url="http://example.invalid")

        tool = Tool(galaxy_server=self.server, id_galaxy="test_tool",
                    name="PhyML", version="3.1", description="")
        Tool.objects.bulk_create([tool])
        self.tool = Tool.objects.get(id_galaxy="test_tool")

        flag = ToolFlag(name="tree", verbose_name="Tree Inference", rank=0)
        ToolFlag.objects.bulk_create([flag])
        self.flag = ToolFlag.objects.get(name="tree")
        self.tool.toolflag_set.add(self.flag)

    def test_tool_and_its_group_appear_on_the_page(self):
        response = self.client.get(reverse('tools'))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, self.tool.name)
        self.assertContains(response, self.flag.verbose_name)


class AddGalaxyKeyCommandTest(TestCase):
    """
    Regression test: addgalaxykey used to build a new GalaxyUser row with a
    plain GalaxyUser(...).save() every run. Re-running it for a user/server
    pair that already has one - e.g. docker-compose's init service running
    again on an existing DB - crashed with "duplicate key value violates
    unique constraint
    galaxy_galaxyuser_user_id_galaxy_server_id_53f2ea9d_uniq" instead of
    refreshing the key. Only caught by actually redeploying against a
    server that already had one imported.
    """

    def setUp(self):
        self.user = User.objects.create_user('admin')
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)

    def test_rerunning_for_the_same_user_and_server_updates_not_duplicates(self):
        call_command('addgalaxykey', user='admin',
                     galaxyurl=self.server.url, galaxykey='key1')
        self.assertEqual(GalaxyUser.objects.count(), 1)
        self.assertEqual(GalaxyUser.objects.get().api_key, 'key1')

        # Simulate re-running the same setup step later (e.g. a redeployed
        # Galaxy with a freshly generated admin key).
        call_command('addgalaxykey', user='admin',
                     galaxyurl=self.server.url, galaxykey='key2')
        self.assertEqual(GalaxyUser.objects.count(), 1)
        self.assertEqual(GalaxyUser.objects.get().api_key, 'key2')
