from __future__ import absolute_import

import logging

from celery import shared_task
from datetime import timedelta

from celery.utils.log import get_task_logger

from workflows.models import Workflow
from galaxy.decorator import galaxy_connection
from django.db import transaction
from django.utils import timezone

def flush_transaction():
    transaction.commit()
    
logger = get_task_logger(__name__)

@shared_task
def deletegalaxyworkflow(workflow_galaxyid):
    """
    Deletes a workflow on the Galaxy server.

    Returns True if the deletion actually succeeded, False otherwise - so
    callers can decide whether it's safe to mark their own bookkeeping as
    deleted, rather than assuming success just because this didn't raise.
    """
    try:
        galaxycon = galaxy_connection()
        galaxycon.nocache = True
        galaxycon.workflows.delete_workflow(workflow_galaxyid)
        return True
    except Exception as e:
        logging.warning("Problem while deleting workflow: %s" % (e))
        return False


# Every day at 2am (same schedule as workspace.tasks.deleteoldgalaxyhistory's
# own per-history workflow cleanup, see below), remove every non-base
# workflow older than 7 days from Galaxy - regardless of whether it was
# ever actually run. This used to only catch workflows with zero
# associated WorkspaceHistory (created but never run) on a 1-day cutoff -
# anything that WAS run stayed in Galaxy forever unless its own specific
# linked WorkspaceHistory happened to independently satisfy
# deleteoldgalaxyhistory's narrower conditions (finished=True, 14-day
# cutoff, workflow FK actually set) at cleanup time. Real usage over years
# left hundreds of thousands of duplicated-category rows never touched by
# either task - see the one-off cleanup this mirrors,
# scripts/cleanup_old_galaxy_workflows.sh.
#
# Deleting a workflow definition doesn't touch its associated History's
# actual data (datasets/job outputs live in the History, not the
# Workflow) - so there's no need to wait for that history's own 14-day
# retention window before dropping the workflow "recipe" that launched
# it. w.delete() also SET_NULLs any WorkspaceHistory.workflow FK pointing
# here (see that field's on_delete), so if deleteoldgalaxyhistory
# processes the same history afterwards it correctly sees workflow=None
# and skips re-deleting anything already gone; if it runs first instead,
# it already marks the Workflow row deleted=True itself, so this task's
# own deleted=False filter skips it in turn - safe regardless of which of
# the two tasks Celery happens to run first.
@shared_task
def deleteoldgalaxyworkflows():
    logger.info("Start old workflow deletion task")
    datecutoff = timezone.now() - timedelta(days=7)
    for w in Workflow.objects.exclude(category='base').filter(date__lte=datecutoff).filter(deleted=False):
        try:
            if deletegalaxyworkflow(w.id_galaxy):
                w.delete()
            # else: leave deleted=False - deletegalaxyworkflow already
            # logged why, and this will be retried on the next run
            # instead of silently disappearing from future cleanups.
        except Exception as e:
            logging.warning(
                "Problem while deleting old workflow %s: %s" %
                (w.id_galaxy, e))
    logger.info("Old workflow deletion task finished")
