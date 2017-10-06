from workflows.models import Workflow, WorkflowStepInformation


class WorkflowDetailMixin(object):
    """
        GET Workflow
    """
    model = Workflow
    object = None

    def fetch_workflow_detail(self, workflow):

        # add galaxy workflows json information
        workflow.json = self.request.galaxy.workflows.show_workflow(workflow_id=workflow.id_galaxy)

        # add workflow sorted tool list information
        workflow.detail = WorkflowStepInformation(workflow.json, tools=self.restrict_toolset).sorted_tool_list
        return workflow

    def get_workflow_detail(self):
        """
            Get Workflow with Galaxy detail
        """
        wk_obj = self.get_object()
        self.fetch_workflow_detail(wk_obj)
        return wk_obj

    def get_workflow(self, detail=True):
        """
            Get Workflow
        """
        # add more information load from Galaxy
        if detail:
            wk_obj = self.get_workflow_detail()
        else:
            wk_obj = self.get_object()

        return wk_obj


class WorkflowDeleteWorkingCopyMixin(object):
    """
    Clean working copy to Galaxy
    """

    def delete_workflow(self, workflow_id):
        wk = Workflow.objects.filter(id_galaxy=workflow_id)
        msg = ''
        # Verify workflow doesnt exist in db
        if not wk:
            msg = self.request.galaxy.workflows.delete_workflow(workflow_id=workflow_id)
        print msg
        return msg


class WorkflowDuplicateMixin(WorkflowDeleteWorkingCopyMixin):
    """
        get_workflow return a working copy from workflow
    """
    copy_workflow = None

    def get_workflow(self):
        """
        Create a working copy of workflow in Galaxy:  import workflow oneclick shared
        """
        if self.copy_workflow is None:
            wk = super(WorkflowDuplicateMixin, self).get_workflow()
            gi = self.request.galaxy
            try:
                # import published shared workflow
                wf_import = gi.workflows.import_shared_workflow(wk.id_galaxy)
                print
            except:
                # make a working copy workflow
                wk_cp = gi.workflows.export_workflow_dict(wk.id_galaxy)
                wf_import = gi.workflows.import_workflow_dict(wk_cp)

            if wf_import:
                # replace workflow with working copy of original workflow shared
                wk.pk = None
                wk.id_galaxy = wf_import.get('id')
                # return obj without pk and new galaxy id
                self.copy_workflow = wk

        return self.copy_workflow

    def clean_copy(self):
        """
        Clean working copy to Galaxy
        """
        if self.copy_workflow:
            return self.delete_workflow(self.copy_workflow.id_galaxy)
