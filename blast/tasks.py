from __future__ import absolute_import

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

from io import StringIO

import logging

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
def launchblast(blastrunid, sequence, prog, db, evalue):
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
            b.database=db
            b.blastprog = prog
            b.status=BlastRun.RUNNING
            b.save()

            result_handle = NCBIWWW.qblast(prog, db, sequence)
            blast_records = NCBIXML.parse(result_handle)

            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < evalue:
                            s = BlastSubject(subject_id=alignment.title.split(" ")[0],
                                             subject_seq=hsp.sbjct.replace("-",""),
                                             blastrun=b)
                            s.save()
            b.status=BlastRun.FINISHED
        else:
            b.status=BlastRun.ERROR
            b.message = "More than one record in the fasta file! %d" % (len(list(records)))
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
        e.deleted = True
        e.blastsubject_set.all().delete()
        e.save()
    logger.info("Old blast deletion task finished")

