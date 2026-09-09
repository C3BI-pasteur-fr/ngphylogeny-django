from django.test import TestCase
from django.urls import reverse

from galaxy.models import Server
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
