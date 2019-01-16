from __future__ import absolute_import

import logging

from celery.schedules import crontab
from celery.decorators import periodic_task
from datetime import date, timedelta, datetime

from celery.utils.log import get_task_logger

from workflows.models import Workflow
from galaxy.decorator import galaxy_connection
from django.db import transaction

def flush_transaction():
    transaction.commit()
    
logger = get_task_logger(__name__)

# Every day at 2am, remove workflows (wkmaker) older than 1 day from galaxy
@periodic_task(run_every=(crontab(hour="02", minute="00", day_of_week="*")))
def deleteoldgalaxyworkflows():
    logger.info("Start old workflow deletion task")
    datecutoff = datetime.now() - timedelta(seconds=60)
    try:
        for w in Workflow.objects.filter(category='automaker').filter(date__lte=datecutoff).filter(deleted=False):
            galaxycon = galaxy_connection()
            galaxycon.nocache = True
            w.delete(galaxycon)
            w.deleted = True
            w.save()
    except Exception as e:
        logging.warning("Problem while deleting history: %s" % (e))
    logger.info("Old workflow deletion task finished")
