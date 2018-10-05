from __future__ import absolute_import

from django.core.mail import send_mail
from django.core.urlresolvers import reverse
from smtplib import SMTPException

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

from io import StringIO

import logging
import re

from celery import shared_task
from celery.decorators import periodic_task
from celery.utils.log import get_task_logger
from celery.schedules import crontab

from .models import BlastRun, BlastSubject

from datetime import date, timedelta, datetime


def flush_transaction():
    transaction.commit()

logger = get_task_logger(__name__)

@shared_task
def launchblast(blastrunid, sequence, prog, db, evalue, coverage):
    """
    Celery task that will launch a blast on the public blast server
    """
    logging.info("Blasting %s with %s on %s" % (sequence, prog, db))
    b = BlastRun.objects.get(id=blastrunid)
    try:
        fasta_io = StringIO(sequence)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        #fasta_io.close()
        i=0
        if len(records) == 1:
            b.query_id=records[0].id
            b.query_seq=records[0].seq
            b.evalue=evalue
            b.coverage = coverage
            b.database=db
            b.blastprog = prog
            b.status=BlastRun.RUNNING
            b.save()

            result_handle = NCBIWWW.qblast(prog, db, sequence)
            blast_records = NCBIXML.parse(result_handle)

            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < evalue and (hsp.align_length / blast_record.query_length) >= coverage :
                            s = BlastSubject(subject_id=alignment.title.split(" ")[0],
                                             subject_seq=hsp.sbjct.replace("-",""),
                                             blastrun=b)
                            s.save()
            b.status=BlastRun.FINISHED
        else:
            b.status=BlastRun.ERROR
            b.message = "More than one record in the fasta file! %d" % (len(list(records)))

        if b.email != None and  re.match(r"[^@]+@[^@]+\.[^@]+", b.email):
            try:
                message = "Dear NGPhylogeny user, \n\n"
                if b.status != b.FINISHED :
                    message= message + "Your NGPhylogeny BLAST job finished with errors.\n\n"
                else:
                    message=message + "Your NGPhylogeny BLAST job finished successfuly.\n"
                please = 'Please visit http://%s%s to check results\n\n' % ("ngphylogeny.fr", reverse('blast_view', kwargs={'pk':b.id}))
                message = message + please
                message = message + "Thank you for using ngphylogeny.fr\n\n"
                message = message + "NGPhylogeny.fr development team.\n"
                send_mail(
                    'NGPhylogeny.fr BLAST results',
                    message,
                    'ngphylogeny@pasteur.fr',
                    [b.email],
                    fail_silently=False,
                )
                print(message)
            except SMTPException as e:
                logging.warning("Problem with smtp server : %s" % (e))
            except Exception as e:
                logging.warning("Unknown Problem while sending e-mail: %s" % (e))
    except Exception as e:
        b.status=BlastRun.ERROR
        b.message=str(e)
    b.save()


# Every day at 2am, clears analyses older than 14 days
@periodic_task(run_every=(crontab(hour="02", minute="00", day_of_week="*")))
def deleteoldblastruns():
    logger.info("Start old blast deletion task")
    datecutoff = datetime.now() - timedelta(days=14)
    for e in BlastRun.objects.filter(deleted=False).filter(date__lte=datecutoff):
        e.soft_delete()
    logger.info("Old blast deletion task finished")
