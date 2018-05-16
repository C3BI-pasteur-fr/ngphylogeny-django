from __future__ import absolute_import

import time
import json
import re

from celery import shared_task
from workspace.models import WorkspaceHistory
from galaxy.decorator import galaxy_connection
from django.core.mail import send_mail
from django.db import transaction
from django.core.urlresolvers import reverse

def flush_transaction():
    transaction.commit()

@shared_task
def monitorworkspace(historyid):
    """
    Celery task that will monitor galaxy workspace
    It will wait for end of execution of all jobs
    and send a mail at the end, if the mail has been
    given by the user.
    
    this task is launched by workspace.views.HistoryDetailView
    only once, when the WorkspaceHistory is marked as "notmonitored"
    """
    galaxycon = galaxy_connection()
    galaxycon.nocache = True
    #print "Monitoring workspace " + historyid
    finished = False
    error = False
    email = None

    while not finished:
        hi = galaxycon.histories.show_history(historyid, contents=True)
        finished = True
        for file in hi:
            if ( 'running' in file.get('state','') or
                 'queued' in file.get('state','') or
                 'new' in file.get('state','')):
                finished = False
            if 'error' in file.get('state',''):
                error = True
                finished = True
        if not finished:
            time.sleep(2)

    w = WorkspaceHistory.objects.get(history=historyid)
    if w and w.email and re.match(r"[^@]+@[^@]+\.[^@]+", w.email):
        try:
            print "worspace informations"
            print w.history
            print w.email
            print w.name
            message = "Dear NGPhylogeny user, \n\n"
            if error:
                message= message + "Your NGPhylogeny job finished with errors.\n\n"
            else:
                message="Your NGPhylogeny job finished successfuly.\n"
            please = 'Please visit http://%s%s to check results\n\n' % ("ngphylogeny.fr", reverse('history_detail', kwargs={'history_id':historyid}))
            message = message + please
            message = message + "Thank you for using ngphylogeny.fr\n\n"
            message = message + "NGPhylogeny.fr development team.\n"
            send_mail(
                'NGPhylogeny.fr results',
                message,
                'ngphylogeny@pasteur.fr',
                [w.email],
                fail_silently=False,
            )
            print(message)
        except:
            print("Unexpected error:", sys.exc_info()[0])
    
