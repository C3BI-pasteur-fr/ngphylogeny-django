from __future__ import absolute_import

import logging

from celery import shared_task
from datetime import timedelta

from celery.utils.log import get_task_logger

from workflows.models import Workflow
from workspace.models import WorkspaceHistory
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


# Every day at 2am, remove workflows older than 1 day from galaxy and that are not
# associated to a workspace (i.e. have not been executed)
@shared_task
def deleteoldgalaxyworkflows():
    logger.info("Start old workflow deletion task")
    datecutoff = timezone.now() - timedelta(days=1)
    for w in Workflow.objects.exclude(category='base').filter(date__lte=datecutoff).filter(deleted=False):
        try:
            # No associated workspace: this workflow was created (e.g. the
            # user opened a submission form) but never actually run - safe
            # to drop.
            if WorkspaceHistory.objects.filter(workflow=w.id).count() == 0:
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
