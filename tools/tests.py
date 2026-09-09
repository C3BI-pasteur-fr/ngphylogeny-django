from django.test import TestCase

from tools.models import Tool


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
