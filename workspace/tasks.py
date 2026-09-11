from __future__ import absolute_import

import time
import json
import re
import logging
import os

from celery import shared_task
from datetime import timedelta

from celery.utils.log import get_task_logger

from smtplib import SMTPException
from workspace.models import WorkspaceHistory
from galaxy.decorator import galaxy_connection
from workflows.tasks import deletegalaxyworkflow

from django.core.mail import send_mail
from django.db import transaction
from django.urls import reverse
from django.core.cache import cache
from django.utils import timezone

LOCK_EXPIRE = 60 * 5 # Lock expires in 5 minutes
LOCK_EXPIRE_SHORT = 9 # Lock expires in 9 seconds


def flush_transaction():
    transaction.commit()
    
logger = get_task_logger(__name__)

@shared_task
def initializeworkspacejob(historyid):
    galaxycon = galaxy_connection()
    galaxycon.nocache = True
    hc = galaxycon.histories.show_history(historyid, contents=True)
    hi = galaxycon.histories.show_history(historyid)
    w = WorkspaceHistory.objects.get(history=historyid)
    w.history_content_json = json.dumps(hc)
    w.history_info_json =  json.dumps(hi)
    w.save()


@shared_task
def launchmonitorworkspaces():
    """
    Celery periodic task that will monitor galaxy workspaces

    It will update content of the workspace django model every minutes

    It will wait for end of execution of all jobs
    and send a mail at the end, if the mail has been
    given by the user.
    """
    for w in WorkspaceHistory.objects.filter(monitored=True, finished=False, deleted=False):
        historyid = w.history
        updateworkspacestatus.delay(historyid)

@shared_task
def updateworkspacestatus(historyid):
    
    ## To be sure that the task is not reexecuted in parallel while
    ## the previous one is still running
    lock_id = "lock_ngphylo_workspacemonitoring_"+historyid
    # We lock this history for 9 seconds, to avoid too frequent refreshs 
    acquire_lock = lambda: cache.add(lock_id, "true", LOCK_EXPIRE_SHORT)
    release_lock = lambda: cache.delete(lock_id)

    if acquire_lock():
        pass
    else:
        return
    
    try:
        galaxycon = galaxy_connection()
        galaxycon.nocache = True
        #print "Monitoring workspace " + historyid
        finished = False
        error = False
        email = None

        hc = galaxycon.histories.show_history(historyid, contents=True)
        hi = galaxycon.histories.show_history(historyid)
        w = WorkspaceHistory.objects.get(history=historyid)
        if w.monitored and not w.finished and not w.deleted:
            w.history_content_json = json.dumps(hc)
            w.history_info_json =  json.dumps(hi)
            w.save()
            if len(hc) > 1:
                finished = True
                for file in hc:
                    if ( 'running' in file.get('state','') or
                         'queued' in file.get('state','') or
                         'new' in file.get('state','')):
                        finished = False
                    if 'error' in file.get('state',''):
                        error = True
                        finished = True
                        break
            if finished:
                w.finished = finished
                logging.warning("history %s finished? %r" % (historyid, w.finished))
                if w and w.email and re.match(r"[^@]+@[^@]+\.[^@]+", w.email):
                    logging.warning("Sending EMail to %s",w.email)
                    try:
                        ngphylohost=os.environ.get('NGPHYLO_HOST')
                        if ngphylohost is None:
                            ngphylohost = "ngphylogeny.fr"
                        citation = "Lemoine F, Correia D, Lefort V, Doppelt-Azeroual O, Mareuil F, Cohen-Boulakia S, Gascuel O\n" \
                                   "NGPhylogeny.fr: new generation phylogenetic services for non-specialists.\n" \
                                   "Nucleic Acids Research 2019 (https://doi.org/10.1093/nar/gkz303).\n"
                        message = "Dear NGPhylogeny user, \n\n"
                        if error:
                            message= message + "Your NGPhylogeny job finished with errors.\n\n"
                        else:
                            message=message + "Your NGPhylogeny job finished successfuly.\n"
                        please = 'Please visit http://%s%s to check results\n\n' % (ngphylohost, reverse('history_detail', kwargs={'history_id':historyid}))
                        message = message + please
                        message = message + "Thank you for using ngphylogeny.fr\n\n"
                        message = message + "NGPhylogeny.fr development team.\n\n"
                        message = message + citation
                        
                        send_mail(
                            'NGPhylogeny.fr results',
                            message,
                            'ngphylogeny@pasteur.fr',
                            [w.email],
                            fail_silently=False,
                        )
                        #print(message)
                    except SMTPException as e:
                        logging.warning("Problem with smtp server : %s" % (e))
                    except Exception as e:
                        logging.warning("Unknown Problem while sending e-mail: %s" % (e))
            w.save()
    except:
        logging.warning('Problem with Galaxy server, will retry later')

    #release_lock()
    
@shared_task
def deletegalaxyhistory(historyid):
    """
    Celery task that will delete an history on the galaxy server in background.

    Returns True if the deletion actually succeeded, False otherwise - so
    callers can decide whether it's safe to mark their own bookkeeping as
    deleted, rather than assuming success just because this didn't raise.
    """
    logging.info("Deleting history %s" % (historyid))
    try:
        galaxycon = galaxy_connection()
        galaxycon.nocache = True
        galaxycon.histories.delete_history(historyid, purge=True)
        return True
    except Exception as e:
        logging.warning("Problem while deleting history: %s" % (e))
        return False


# Every day at 2am, clears analyses older than 14 days
@shared_task
def deleteoldgalaxyhistory():
    logger.info("Start old workspace deletion task")
    datecutoff = timezone.now() - timedelta(days=14)
    for e in WorkspaceHistory.objects.filter(deleted=False).filter(finished=True).filter(created_date__lte=datecutoff):
        try:
            workflow_deleted = True
            if e.workflow is not None:
                workflow_deleted = deletegalaxyworkflow(e.workflow.id_galaxy)
                if workflow_deleted:
                    e.workflow.deleted = True
                    e.workflow.save()

            history_deleted = deletegalaxyhistory(e.history)

            if workflow_deleted and history_deleted:
                e.deleted = True
                # Also drop the cached Galaxy history contents: once the
                # underlying Galaxy history is purged these would just be
                # stale JSON forever, and there's no other cleanup that
                # reclaims this (potentially sizeable) storage.
                e.history_content_json = ""
                e.history_info_json = ""
                e.save()
            else:
                # Leave deleted=False: deletegalaxyworkflow/
                # deletegalaxyhistory already logged why, and this will be
                # picked up again on the next run instead of being marked
                # "cleaned up" while the data may still exist on Galaxy.
                logging.warning(
                    "Could not fully delete workspace history %s "
                    "(workflow_deleted=%r, history_deleted=%r) - will "
                    "retry on the next run" %
                    (e.history, workflow_deleted, history_deleted))
        except Exception as ex:
            logging.warning(
                "Problem while deleting old workspace %s: %s" %
                (e.history, ex))
    logger.info("Old workspace deletion task finished")
