from __future__ import absolute_import

from django.db.models import Q
from django.core.mail import send_mail
from django.core.urlresolvers import reverse
from django.core.cache import cache

from smtplib import SMTPException

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

from io import StringIO
import shutil
import logging
import re
import time
import tempfile

from celery import shared_task
from celery.decorators import periodic_task
from celery.utils.log import get_task_logger
from celery.schedules import crontab

from datetime import timedelta, datetime

from galaxy.decorator import galaxy_connection
from bioblend.galaxy.tools.inputs import inputs

from .models import BlastRun, BlastSubject
from .msa import PseudoMSA

from utils import biofile

logger = get_task_logger(__name__)

LOCK_EXPIRE = 60 * 5 # Lock expires in 5 minutes


## It should be alone on a celery queue with only 1 cpu
## Otherwise, may run too many jobs on ncbi server
@shared_task
def launch_ncbi_blast(blastrunid, sequence, prog, db, evalue, coverage, maxseqs):
    """
    Celery task that will launch a blast on the public blast server
    """
    logging.info("Blasting %s on %s" % (prog, db))
    b = BlastRun.objects.get(id=blastrunid)
    try:
        fasta_io = StringIO(sequence)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        if len(records) == 1:
            b.query_id = biofile.cleanseqname(records[0].id)
            b.query_seq = records[0].seq
            b.evalue = evalue
            b.coverage = coverage
            b.database = db
            b.blastprog = prog
            b.maxseqs = maxseqs
            b.status = BlastRun.RUNNING
            b.save()

            blast_inputtype = BlastRun.blast_inputtype(BlastRun.NCBI, prog)
            blast_type = BlastRun.blast_type(BlastRun.NCBI, prog)

            # We check alphabet of given sequence
            if ((blast_inputtype == "nt" and not biofile.check_nt(b.query_seq)) or
                (blast_inputtype == "aa" and not biofile.check_aa(b.query_seq))):
                b.status = BlastRun.ERROR
                b.message = "The given sequence has the wrong alphabet. Program %s expects %s sequence" % (
                    blast_type, blast_inputtype)
            else:
                rh = NCBIWWW.qblast(prog, db, sequence)
                tmp_file = tempfile.NamedTemporaryFile()
                shutil.copyfileobj(rh, tmp_file)
                tmp_file.flush()

                query_seq_bk = b.query_seq
                frame = 1
                if blast_type == 'blastx' or blast_type == 'tblastx' :
                    frame=majorityQueryFrame(tmp_file.name)
                    b.query_seq = biofile.translate(str(b.query_seq), frame)
                    b.save()
    
                result_handle = open(tmp_file.name, "r")
                blast_records = NCBIXML.parse(result_handle)
                ms = PseudoMSA(b.query_id, b.query_seq, query_seq_bk, frame, blast_type)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            e_val = hsp.expect
                            leng = float(hsp.align_length) / float(len(str(b.query_seq)))
                            if e_val < evalue and leng >= coverage:
                                ms.add_hsp(alignment.title.split(" ")[0], hsp)

                nseq=0
                
                if blast_type == 'blastx' or blast_type == 'tblastx' :
                    ms.crop_alignment(maxseqs)
                    b.query_seq = "".join(ms.query_seq)
                    b.save()
                
                for id, seq, fullseq in ms.first_n_max_score_sequences(maxseqs):
                    s = BlastSubject(subject_id=id,
                                     subject_seq=seq,
                                     subject_fullseq=fullseq,
                                     blastrun=b)
                    s.save()
                    nseq+=1

                if nseq>0:
                    b.tree = b.build_nj_tree()
                    b.status = BlastRun.FINISHED
                    b.save()
                else:
                    b.status = BlastRun.ERROR
                    b.message = "Blast Search returned no results"
                    b.save()
                    
        else:
            b.status = BlastRun.ERROR
            b.message = "More than one record in the fasta file! %d" % (
                len(list(records)))

        if b.email is not None and re.match(r"[^@]+@[^@]+\.[^@]+", b.email):
            try:
                message = "Dear NGPhylogeny user, \n\n"
                if b.status != b.FINISHED:
                    message = message + "Your NGPhylogeny BLAST job finished with errors.\n\n"
                else:
                    message = message + "Your NGPhylogeny BLAST job finished successfuly.\n"
                please = 'Please visit http://%s%s to check results\n\n' % (
                    "ngphylogeny.fr", reverse('blast_view', kwargs={'pk': b.id}))
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
            except SMTPException as e:
                logging.warning("Problem with smtp server : %s" % (e))
            except Exception as e:
                logging.warning(
                    "Unknown Problem while sending e-mail: %s" % (e))
    except Exception as e:
        logging.exception(str(e))
        b.status = BlastRun.ERROR
        b.message = str(e)
    b.save()
    time.sleep(30)

@shared_task
def launch_pasteur_blast(blastrunid, sequence, prog, db, evalue, coverage, maxseqs):
    """
    Celery task that will launch a blast on the pasteur Galaxy Server
    """
    logging.info("Blasting %s with %s on %s" % (sequence, prog, db))
    b = BlastRun.objects.get(id=blastrunid)
    try:
        fasta_io = StringIO(sequence)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        if len(records) == 1:
            galaxycon = galaxy_connection()
            galaxycon.nocache = True
            history = galaxycon.histories.create_history(name="BlastXplorer")
            
            b.history = history.get("id")
            b.query_id = biofile.cleanseqname(records[0].id)
            b.query_seq = records[0].seq
            b.evalue = evalue
            b.coverage = coverage
            b.database = db
            b.blastprog = prog
            b.maxseqs = maxseqs
            b.status = BlastRun.PENDING
            b.save()

            blast_type = BlastRun.blast_type(BlastRun.PASTEUR, prog)
            blast_inputtype = BlastRun.blast_inputtype(BlastRun.PASTEUR, prog)
            
            # We check alphabet of given sequence
            if ((blast_inputtype == "nt" and not biofile.check_nt(b.query_seq)) or
                (blast_inputtype == "aa" and not biofile.check_aa(b.query_seq))):
                b.status = BlastRun.ERROR
                b.message = "The given sequence has the wrong alphabet. Program %s expects %s sequence" % (
                    blast_type, blast_inputtype)
            elif blast_type is not None:
                tmp_file = tempfile.NamedTemporaryFile()
                tmp_file.write(sequence)
                tmp_file.flush()
                if biofile.is_fasta_one_seq(tmp_file.name):
                    ## Upload input query file to galaxy
	            outputs = galaxycon.tools.upload_file(path=tmp_file.name,file_name="blastinput.fasta",history_id=history.get("id"),file_type="fasta")
	            file_id = outputs.get('outputs')[0].get('id')
                    ## Configuring job
	            tool_inputs=inputs()
	            tool_inputs.set_dataset_param("query",file_id)
	            tool_inputs.set_param("db_opts|database", db)
	            tool_inputs.set_param("blast_type", blast_type)
	            tool_inputs.set_param("evalue_cutoff", evalue)
	            tool_inputs.set_param("output|out_format", "5")
                    ## Running blast job
	            outputs=galaxycon.tools.run_tool(history_id=history.get("id"),tool_id=prog,tool_inputs=tool_inputs)
                    b.history_fileid = outputs.get("outputs")[0].get("id")
                else:
                    b.status=BlastRun.ERROR
                    b.message="Bad input FASTA file format"
            else:
                b.status=BlastRun.ERROR
                b.message="Wrong blast program %s" % (prog)
            b.save()
        else:
            b.status = BlastRun.ERROR
            b.message = "More than one record in the fasta file! %d" % (
                len(list(records)))
    except Exception as e:
        logging.exception(str(e))
        b.status = BlastRun.ERROR
        b.message = str(e)
    b.save()
    time.sleep(30)
    

@shared_task
def build_tree(blastrunid):
    try:
        b = BlastRun.objects.get(id=blastrunid)
        b.status = BlastRun.RUNNING
        b.tree = ""
        b.save()
        b.tree = b.build_nj_tree()
        b.status = BlastRun.FINISHED
        b.save()
    except Exception as e:
        logging.exception(str(e))
        b.status = BlastRun.ERROR
        b.message = str(e)
        b.save()

@periodic_task(run_every=(crontab(hour="02", minute="00", day_of_week="*")))
def deleteoldblastruns():
    """
    Every day at 2am, clears analyses older than 14 days
    """
    logger.info("Start old blast deletion task")
    datecutoff = datetime.now() - timedelta(days=14)
    for e in BlastRun.objects.filter(deleted=False).filter(date__lte=datecutoff):
        if e.history != "":
            deletegalaxyhistory(e.history)
        e.soft_delete()
        e.save()
    logger.info("Old blast deletion task finished")


@periodic_task(run_every=(crontab(hour="*", minute="*", day_of_week="*")))
def checkblastruns():
    """
    Every minutes, check running pasteur blast runs
    """
    logger.info("Start pasteur blast task check")

    ## To be sure that the task is not reexecuted in parallel while
    ## the previous one is still running
    lock_id = "lock_ngphylo_blastmonitoring"
    acquire_lock = lambda: cache.add(lock_id, "true", LOCK_EXPIRE)
    release_lock = lambda: cache.delete(lock_id)

    if acquire_lock():
        pass
    else:
        return

    try:
        galaxycon = galaxy_connection()
        galaxycon.nocache = True
        
        for b in BlastRun.objects.filter(deleted=False, server=BlastRun.PASTEUR).filter(Q(status=BlastRun.PENDING) | Q(status=BlastRun.RUNNING)):
            # State of the output file we want (blast XML)
            dataset=galaxycon.histories.show_dataset(b.history,b.history_fileid)
            state=dataset.get('state')
            infos=dataset.get('misc_info')
            b.message=infos
    
            if state == 'ok':
                b.status=BlastRun.FINISHED
                blast_type = BlastRun.blast_type(BlastRun.PASTEUR, b.blastprog)
                ## Download the result file from galaxy first...
                tmp_file = tempfile.NamedTemporaryFile()
                galaxycon.datasets.download_dataset(b.history_fileid,tmp_file.name,False)
                query_seq_bk = b.query_seq
                frame = 1
                if blast_type == 'blastx' or blast_type == 'tblastx' :
                    frame=majorityQueryFrame(tmp_file.name)
                    b.query_seq = biofile.translate(str(b.query_seq), frame)
                    b.save()
                
                result_handle = open(tmp_file.name, "r")
                blast_records = NCBIXML.parse(result_handle)
                ms = PseudoMSA(b.query_id, b.query_seq, query_seq_bk, frame, blast_type)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            e_val = hsp.expect
                            leng = float(hsp.align_length) / float(len(b.query_seq))
                            if e_val < b.evalue and leng >= b.coverage:
                                ms.add_hsp(biofile.newick_clean(alignment.title), hsp)

                if blast_type == 'blastx' or blast_type == 'tblastx' :
                    ms.crop_alignment(b.maxseqs)
                    b.query_seq = "".join(ms.query_seq)
                    b.save()

                nseq=0
                for id, seq, fullseq in ms.first_n_max_score_sequences(b.maxseqs):
                    s = BlastSubject(subject_id=id,
                                     subject_seq=seq,
                                     subject_fullseq=fullseq,
                                     blastrun=b)
                    s.save()
                    nseq+=1

                if nseq>0:
                    b.tree = b.build_nj_tree()
                    b.status = BlastRun.FINISHED
                    b.save()
                else:
                    b.status = BlastRun.ERROR
                    b.message = "Blast Search returned no results"
                    b.save()
            elif state == 'queued' or state == 'new':
                b.status=BlastRun.PENDING
            elif state == 'running':
                b.status=BlastRun.RUNNING
            else:
                b.status=BlastRun.ERROR
            b.save()
    
            if b.email is not None and re.match(r"[^@]+@[^@]+\.[^@]+", b.email) and (b.status == BlastRun.ERROR or b.status == BlastRun.FINISHED):
                try:
                    message = "Dear NGPhylogeny user, \n\n"
                    if b.status != b.FINISHED:
                        message = message + "Your NGPhylogeny BLAST job finished with errors.\n\n"
                    else:
                        message = message + "Your NGPhylogeny BLAST job finished successfuly.\n"
                    please = 'Please visit http://%s%s to check results\n\n' % (
                        "ngphylogeny.fr", reverse('blast_view', kwargs={'pk': b.id}))
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
                except SMTPException as e:
                    logging.warning("Problem with smtp server : %s" % (e))
                except Exception as e:
                    logging.warning(
                        "Unknown Problem while sending e-mail: %s" % (e))
    except Exception as e:
        b.status=BlastRun.ERROR
        b.message=str(e)
        b.save()
        logger.info("Error while checking blast run: %s" % (e))
        logging.exception("message")

    release_lock()
    logger.info("Pasteur blast runs checked")

@shared_task
def deletegalaxyhistory(historyid):
    """
    Celery task that will delete an history on the galaxy server in background
    """
    logging.info("Deleting history %s" % (historyid))
    try:
        galaxycon = galaxy_connection()
        galaxycon.nocache = True
        galaxycon.histories.delete_history(historyid, purge=True)
    except Exception as e:
        logging.warning("Problem while deleting history: %s" % (e))


def majorityQueryFrame(blastfile):
    """
    It takes a blast result file and returns the query frame that
    is the most frequent in all HSP
    """
    frames = dict()
    result_handle = open(blastfile, "r")
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.frame[0] in frames:
                    frames[hsp.frame[0]]+=1
                else:
                    frames[hsp.frame[0]]=1

    max_frame = None
    nb_frames = 0
    for k,v in frames.items():
        if nb_frames == 0 or nb_frames<v:
            max_frame = k
            nb_frames = v
    return max_frame

