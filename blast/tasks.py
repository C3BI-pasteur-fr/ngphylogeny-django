from __future__ import absolute_import

from django.conf import settings
from django.db.models import Q
from django.core.cache import cache
from django.utils import timezone

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
from celery.exceptions import SoftTimeLimitExceeded
from celery.utils.log import get_task_logger

from datetime import timedelta

from galaxy.decorator import galaxy_connection
from bioblend.galaxy.tools.inputs import inputs

from .emails import send_blast_completion_email
from .models import BlastRun, BlastSubject
from .msa import PseudoMSA

from utils import biofile

logger = get_task_logger(__name__)

LOCK_EXPIRE = 60 * 5 # Lock expires in 5 minutes

# checkblastruns() gives up on a Pasteur run still pending/running past
# this age. There's no timeout at all on the actual Galaxy-side blast
# computation otherwise (unlike launch_ncbi_blast/launch_pasteur_blast's
# own soft_time_limit/time_limit, which only cover *submitting* the
# job) - a real search against a big database (e.g. blastn vs nt) can
# legitimately take a while, but a run that's still "running" after
# this long is more likely stuck on Galaxy's/the cluster's side than
# genuinely still computing. 3 hours is a guess at "generous enough for
# a real, slow-but-legitimate search, bounded enough to actually
# recover" - adjust based on real observed run times if this turns out
# to be too tight or too loose. Reads settings.
# NGPHYLO_PASTEUR_BLAST_STALE_HOURS (settings/base.py) - configurable
# via the PASTEUR_BLAST_STALE_HOURS GitLab CI/CD variable, 3 if unset
# (this task's original hardcoded value).
PASTEUR_RUN_STALE_AFTER = timedelta(
    hours=settings.NGPHYLO_PASTEUR_BLAST_STALE_HOURS)


## It should be alone on a celery queue with only 1 cpu
## Otherwise, may run too many jobs on ncbi server
#
# time_limit/soft_time_limit: NCBIWWW.qblast() below polls NCBI in a bare
# `while True:` loop with no timeout or retry cap of its own (checked
# directly against biopython 1.70's actual source, this project's pinned
# version - see CLAUDE.md's "Known dependency ceilings") - it relies
# entirely on NCBI eventually sending a recognizable "READY" (or
# no-Status) response. If NCBI is slow, the query gets stuck server-side,
# or a response comes back in a shape this old biopython version doesn't
# recognize as done, this call - and therefore this task - hangs
# indefinitely, with nothing in this codebase to ever notice or recover.
# Since this queue is deliberately meant to run one task at a time (see
# above), one stuck run blocks every subsequent NCBI BLAST submission
# behind it too, forever, until someone manually restarts the
# celery-worker process. soft_time_limit raises SoftTimeLimitExceeded
# inside the task (caught below, so the run gets marked ERROR with a
# clear message and the worker moves on to the next queued task);
# time_limit is Celery's own hard SIGKILL backstop shortly after, in case
# the soft one doesn't get a chance to run (e.g. blocked in a C
# extension). 10 minutes is a guess at "generous enough for a real,
# slow-but-legitimate NCBI search, bounded enough to actually recover" -
# adjust based on real observed run times if this turns out to be too
# tight or too loose.
@shared_task(soft_time_limit=600, time_limit=660)
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
            b.query_length = len(records[0].seq)
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
                # Same branded HTML template as the workflow job-completion
                # email (workspace.emails), not a hand-built plain-text
                # message - see CLAUDE.md's "BLAST completion email" note.
                send_blast_completion_email(b, b.email)
            except SMTPException as e:
                logging.warning("Problem with smtp server : %s" % (e))
            except Exception as e:
                logging.warning(
                    "Unknown Problem while sending e-mail: %s" % (e))
    except SoftTimeLimitExceeded:
        # See this task's own soft_time_limit comment above - a clear,
        # specific message here beats the generic except below's bare
        # str(e) (which for this exception is just an empty
        # "SoftTimeLimitExceeded()").
        logging.warning(
            "NCBI BLAST run %s exceeded the time limit and was aborted" %
            (blastrunid))
        b.status = BlastRun.ERROR
        b.message = ("NCBI took too long to respond and this search was "
                      "aborted. NCBI's public server can be slow or "
                      "congested - please try again, or try a smaller/"
                      "more specific query.")
    except Exception as e:
        logging.exception(str(e))
        b.status = BlastRun.ERROR
        b.message = str(e)
    b.save()
    time.sleep(30)

# soft_time_limit/time_limit: same reasoning as launch_ncbi_blast's own
# comment above, but the likely hang location differs - the blast
# computation itself runs asynchronously on Galaxy once submitted, and
# is separately monitored (with no timeout of its own - see CLAUDE.md's
# BLAST notes) by checkblastruns(). A hang here is most likely in the
# (network-bound) create_history/upload_file/run_tool calls that submit
# the job in the first place. 10 minutes, matching launch_ncbi_blast.
@shared_task(soft_time_limit=600, time_limit=660)
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
            b.query_length = len(records[0].seq)
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
                # mode='w': sequence is a plain str (the pasted/uploaded
                # FASTA text) - NamedTemporaryFile() defaults to binary
                # mode, which raises "a bytes-like object is required, not
                # 'str'" here. Same Python 2->3 bug class as CLAUDE.md's
                # "Code paths only a real Galaxy run exercises" section.
                tmp_file = tempfile.NamedTemporaryFile(mode='w')
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
    except SoftTimeLimitExceeded:
        logging.warning(
            "Pasteur BLAST run %s exceeded the time limit and was "
            "aborted" % (blastrunid))
        # b.history (an in-memory attribute set right after
        # create_history() returns, whether or not it's been saved yet -
        # see above) tells us whether a Galaxy history actually got
        # created before the timeout fired. Clean it up rather than
        # leaving an orphaned history nothing will ever reference again
        # - queued separately (not called directly) so this already
        # timed-out task doesn't also block on deleting it.
        if b.history:
            deletegalaxyhistory.delay(b.history)
        b.status = BlastRun.ERROR
        b.message = ("Submitting this search to the Pasteur Galaxy "
                      "server took too long and it was aborted. Please "
                      "try again.")
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

@shared_task
def deleteoldblastruns():
    """
    Every day at 2am, clears analyses older than BlastRun.RETENTION_DAYS
    """
    logger.info("Start old blast deletion task")
    # timezone.now, not datetime.now: USE_TZ=True is on - see BlastRun.date's
    # own comment in blast/models.py for the same fix/reasoning.
    datecutoff = timezone.now() - timedelta(days=BlastRun.RETENTION_DAYS)
    for e in BlastRun.objects.filter(deleted=False).filter(date__lte=datecutoff):
        if e.history != "":
            # Queued, not called directly - same reasoning as the
            # submission/staleness timeouts' own cleanup: don't let one
            # slow/unresponsive Galaxy history-delete hold up the rest of
            # this batch.
            deletegalaxyhistory.delay(e.history)
        e.soft_delete()
        # Re-derived here (not just trusted from submission time) so
        # rows predating BlastRun.query_length, or from a launch path
        # that somehow never set it, still keep this before it's lost
        # for good.
        if e.query_length is None:
            e.query_length = len(e.query_seq or "")
        # Frees space on rows old enough to be cleaned up anyway - safe
        # for the daily report (workspace/reports.py), which only ever
        # reads BlastRun's date/deleted/id, never query_seq/tree.
        e.query_seq = ""
        e.tree = ""
        e.save()
    logger.info("Old blast deletion task finished")


@shared_task
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
    except Exception as e:
        logger.info("Error while connecting to galaxy: %s" % (e))
        logging.exception("message")
        release_lock()
        return

    # Excludes history_fileid='': launch_pasteur_blast() saves the run as
    # PENDING right after creating its Galaxy history, but only sets
    # history_fileid afterwards, once the (network-bound) file upload +
    # tool run calls complete - a real race with this task's own 1-minute
    # schedule. Without this exclude, show_dataset(b.history, '') turns
    # into a GET on Galaxy's history *contents list* endpoint (trailing
    # empty dataset id) instead of a single dataset, which returns a
    # list, not a dict - crashing with "'list' object has no attribute
    # 'get'". That row is picked up again on the next pass once
    # history_fileid is set.
    #
    # Each run is also processed in its own try/except: previously the
    # entire loop shared one try/except, so a single run's failure (this
    # race included) silently aborted checking of every other
    # pending/running run in the same pass too.
    for b in BlastRun.objects.filter(
            deleted=False, server=BlastRun.PASTEUR
        ).exclude(history_fileid='').filter(
            Q(status=BlastRun.PENDING) | Q(status=BlastRun.RUNNING)):
        try:
            if timezone.now() - b.date > PASTEUR_RUN_STALE_AFTER:
                logging.warning(
                    "Pasteur BLAST run %s has been pending/running for "
                    "over %s - giving up on it" % (b.id, PASTEUR_RUN_STALE_AFTER))
                # Queued separately (not called directly), same reasoning
                # as launch_pasteur_blast's own SoftTimeLimitExceeded
                # cleanup: don't add another blocking Galaxy call to a
                # run we've already decided to abandon.
                if b.history:
                    deletegalaxyhistory.delay(b.history)
                b.status = BlastRun.ERROR
                b.message = (
                    "This search has been running on the Pasteur Galaxy "
                    "server for longer than expected and was aborted. "
                    "Please try again, possibly with a smaller/more "
                    "specific query.")
                b.save()
                continue

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
                    # Same branded HTML template as the workflow
                    # job-completion email (workspace.emails) - see
                    # CLAUDE.md's "BLAST completion email" note.
                    send_blast_completion_email(b, b.email)
                except SMTPException as e:
                    logging.warning("Problem with smtp server : %s" % (e))
                except Exception as e:
                    logging.warning(
                        "Unknown Problem while sending e-mail: %s" % (e))
        except Exception as e:
            logging.warning(
                "Problem while checking blast run %s: %s" % (b.id, e))
            logging.exception("message")
            b.status=BlastRun.ERROR
            b.message=str(e)
            b.save()

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

