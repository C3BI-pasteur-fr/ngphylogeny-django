from __future__ import absolute_import

import logging

from celery import shared_task
from celery.schedules import crontab
from celery.decorators import periodic_task
from datetime import date, timedelta, datetime

from celery.utils.log import get_task_logger

from workflows.models import Workflow
from workspace.models import WorkspaceHistory
from galaxy.decorator import galaxy_connection
from django.db import transaction

def flush_transaction():
    transaction.commit()
    
logger = get_task_logger(__name__)

@shared_task
def deletegalaxyworkflow(workflow_galaxyid):
    try:
        galaxycon = galaxy_connection()
        galaxycon.nocache = True
        galaxycon.workflows.delete_workflow(workflow_galaxyid)
    except Exception as e:
        logging.warning("Problem while deleting workflow: %s" % (e))


# Every day at 2am, remove workflows older than 1 day from galaxy and that are not
# associated to a workspace (i.e. have not been executed)
@periodic_task(run_every=(crontab(hour="02", minute="00", day_of_week="*")))
def deleteoldgalaxyworkflows():
    logger.info("Start old workflow deletion task")
    datecutoff = datetime.now() - timedelta(days=1)
    try:
        for w in Workflow.objects.exclude(category='base').filter(date__lte=datecutoff).filter(deleted=False):
            # No associated workspace
            if WorkspaceHistory.objects.filter(workflow = w.id).count() == 0:
                galaxycon = galaxy_connection()
                galaxycon.nocache = True
                galaxycon.workflows.delete_workflow(w.id_galaxy)
                w.delete()
    except Exception as e:
        logging.warning("Problem while deleting history: %s" % (e))
    logger.info("Old workflow deletion task finished")
