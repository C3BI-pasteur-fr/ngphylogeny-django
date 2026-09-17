from __future__ import absolute_import

import time
import json
import re
import logging

from celery import shared_task
from datetime import timedelta

from celery.utils.log import get_task_logger

from email.mime.image import MIMEImage
from smtplib import SMTPException
from workspace.models import WorkspaceHistory
from galaxy.decorator import galaxy_connection
from workflows.tasks import deletegalaxyworkflow

from django.conf import settings
from django.core.mail import EmailMultiAlternatives
from django.db import transaction
from django.core.cache import cache
from django.utils import timezone

from workspace.emails import send_job_completion_email
from workspace.reports import render_report_email

LOCK_EXPIRE = 60 * 5 # Lock expires in 5 minutes
LOCK_EXPIRE_SHORT = 9 # Lock expires in 9 seconds

# A stuck Galaxy job (a hung cluster node, a tool that never returns) used
# to leave a WorkspaceHistory as finished=False forever - launchmonitorworkspaces
# just keeps polling it, and deleteoldgalaxyhistory only ever looks at
# finished=True rows, so nothing ever cleaned it up either. Same class of
# gap blast.tasks.checkblastruns() was already fixed for
# (PASTEUR_RUN_STALE_AFTER); this is that same fix for regular workflow/
# tool runs. Reference point is WorkspaceHistory.created_date, same
# choice as BlastRun.date there.
WORKFLOW_RUN_STALE_AFTER = timedelta(hours=24)


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

def cancel_stale_jobs(galaxycon, historyid):
    """
    Cancels every job for this history still in a non-terminal state
    (new/queued/running) - used once a run has been going for longer
    than WORKFLOW_RUN_STALE_AFTER. bioblend's cancel_job() deletes that
    *job's own* not-yet-materialized output dataset, but this is a
    per-job operation - every other, already-completed dataset in this
    history (earlier steps that finished fine) is untouched and stays
    downloadable, same as CheckBlastRunsTest's own approach for BLAST.
    One `get_jobs` call, not one per state, to keep this endpoint's own
    Galaxy load down (see CLAUDE.md's step-chain section for why that
    matters here).
    """
    stale_states = ('new', 'queued', 'running')
    for job in galaxycon.jobs.get_jobs(history_id=historyid):
        if job.get('state') in stale_states:
            try:
                galaxycon.jobs.cancel_job(job['id'])
            except Exception as e:
                logging.warning(
                    "Could not cancel stale job %s (history %s): %s",
                    job.get('id'), historyid, e)


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
            stale_states = ('new', 'queued', 'running')
            still_running = any(
                any(s in file.get('state', '') for s in stale_states)
                for file in hc)
            if still_running and timezone.now() - w.created_date > WORKFLOW_RUN_STALE_AFTER:
                logging.warning(
                    "history %s has been running for longer than %s - "
                    "cancelling its still-running/queued Galaxy jobs",
                    historyid, WORKFLOW_RUN_STALE_AFTER)
                cancel_stale_jobs(galaxycon, historyid)
                # Force these to error locally rather than re-fetching to
                # see how Galaxy itself now reports a cancelled job's
                # dataset - a deliberate, known outcome we're causing
                # here, not something to infer indirectly.
                for file in hc:
                    if any(s in file.get('state', '') for s in stale_states):
                        file['state'] = 'error'

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
                        send_job_completion_email(historyid, w.email, error)
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


@shared_task
def send_daily_report():
    """
    Emails the daily workflow-usage report (workspace/reports.py) to
    settings.NGPHYLO_REPORT_RECIPIENTS - see CELERY_BEAT_SCHEDULE in
    settings/base.py for the schedule (8am UTC daily). No-ops (just logs)
    if no recipients are configured, so this is safe to leave enabled on
    deployments that don't want it.
    """
    recipients = settings.NGPHYLO_REPORT_RECIPIENTS
    if not recipients:
        logger.info(
            "NGPHYLO_REPORT_RECIPIENTS not set - skipping daily report")
        return

    try:
        html, images = render_report_email()
        msg = EmailMultiAlternatives(
            'NGPhylogeny.fr - Daily workflow report',
            'This report is only available in HTML - please enable '
            'HTML email to view it.',
            settings.NGPHYLO_REPORT_FROM_EMAIL,
            recipients,
        )
        msg.attach_alternative(html, 'text/html')
        # multipart/related, not the send_mail()-style multipart/mixed
        # default - required for mail clients to treat the images below as
        # inline (cid:-referenced) parts of the HTML rather than as
        # ordinary file attachments alongside it.
        msg.mixed_subtype = 'related'
        for cid, png_bytes in images.items():
            image = MIMEImage(png_bytes, 'png')
            image.add_header('Content-ID', '<%s>' % cid)
            image.add_header('Content-Disposition', 'inline',
                              filename='%s.png' % cid)
            msg.attach(image)
        msg.send(fail_silently=False)
        logger.info("Daily report sent to %s" % (recipients,))
    except SMTPException as e:
        logging.warning("Problem sending daily report: %s" % (e))
    except Exception as e:
        logging.warning("Problem building/sending daily report: %s" % (e))
