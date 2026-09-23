# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os

from django.db import models
from django.conf import settings
from django.utils import timezone

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo

import uuid
import textwrap
import logging
from io import StringIO

import re

from .trees import prot_dist, nucl_dist


class BlastRun(models.Model):
    PENDING = 'P'
    RUNNING = 'R'
    FINISHED = 'F'
    ERROR = 'E'
    RUNSTATUS = (
        (PENDING, 'Pending'),
        (RUNNING, 'Running'),
        (FINISHED, 'Finished'),
        (ERROR, 'Error')
    )

    PASTEUR = 'pasteur'
    NCBI = 'ncbi'
    BLASTSERVERS =(
        (PASTEUR, 'Pasteur'),
        (NCBI, 'NCBI'),
    )
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    email = models.CharField(null=True, max_length=100)
    # timezone.now, not datetime.now: USE_TZ=True is on, and a naive
    # datetime here throws "RuntimeWarning: DateTimeField BlastRun.date
    # received a naive datetime while time zone support is active" on
    # every save - same bug class already fixed for Workflow.date (see
    # CLAUDE.md's "Workflow duplicates and the Celery cleanup jobs").
    # default= is Python-side only, not part of the DB schema, so this
    # needs no migration/manual ALTER TABLE on an already-deployed DB.
    date = models.DateTimeField(default=timezone.now, blank=True)
    query_id = models.CharField(null=True, max_length=1000)
    query_seq = models.TextField(null=True)
    # Set at submission time (blast/tasks.py's launch_ncbi_blast/
    # launch_pasteur_blast) and re-derived defensively at cleanup time
    # (deleteoldblastruns(), right before query_seq is cleared to free
    # space) so the sequence length survives even once the sequence
    # itself doesn't - null=True since existing rows predate this field
    # and are never backfilled (migrations aren't committed/data-
    # migrated in this project - see CLAUDE.md).
    query_length = models.PositiveIntegerField(null=True, blank=True)
    evalue = models.FloatField(default=0.00001)
    coverage = models.FloatField(default=0.8)
    maxseqs = models.PositiveIntegerField(default=10)
    database = models.CharField(max_length=100, default='swissprot')
    blastprog = models.CharField(max_length=100, default='blastp')
    # max_length=250, not 20: matches workflows.Workflow.id_galaxy's
    # existing convention for the same kind of value (an opaque
    # Galaxy-provided encoded id). 20 was too tight for this Galaxy
    # server's actual dataset ids - launch_pasteur_blast() would run the
    # blast job for real, then crash saving the result:
    # "django.db.utils.DataError: value too long for type character
    # varying(20)" on history_fileid, leaving the run stuck showing
    # PENDING in NGPhylogeny while it kept running/finished on Galaxy.
    history = models.CharField(max_length=250) # If pasteur blast: galaxy history id
    history_fileid= models.CharField(max_length=250) # If pasteur blast: output file galaxy id
    status = models.CharField(max_length=1, default=PENDING, choices=RUNSTATUS)
    server = models.CharField(max_length=50, default=NCBI, choices=BLASTSERVERS)
    message = models.TextField(null=True)
    deleted = models.BooleanField(default=False)
    tree = models.TextField(null=True)

    class Meta:
        indexes = [
            # Same reasoning as WorkspaceHistory's own
            # wsph_running_jobs_idx (see workspace/models.py) -
            # workspace.views.running_jobs_view's BlastRun query filters
            # on exactly these two columns, with no index. A partial
            # index over just the still-pending/running rows stays tiny
            # regardless of how large blast_blastrun grows overall.
            models.Index(
                fields=['date'],
                # Literal 'P'/'R', not the PENDING/RUNNING class
                # constants - a nested Meta class body doesn't have
                # access to BlastRun's own namespace, only the module's.
                condition=models.Q(status__in=['P', 'R'], deleted=False),
                name='blastrun_running_idx',
            ),
        ]

    def format_sequence(self):
        return re.sub("\*$","",('\n'.join(textwrap.wrap(self.query_seq, 60)))).rstrip()

    def sequence(self):
        return re.sub("\*$","",self.query_seq.rstrip())
        
    def to_fasta(self):
        """
        Returns all full sequences in Fasta format.
        Considers also insertions in query sequence
        """
        fasta = ">%s\n" % self.query_id
        fasta += "%s\n" % self.format_sequence()
        for s in self.blastsubject_set.all():
            fasta += ">%s\n" % s.subject_id
            fasta += "%s\n" % s.format_fullsequence()
        return fasta

    def status_str(self):
        for (code, desc) in self.RUNSTATUS:
            if self.status == code:
                return desc
        return 'Error'

    def server_str(self):
        # Was comparing against self.status (a RUNSTATUS code, e.g. 'P'/
        # 'R') instead of self.server (a BLASTSERVERS code, 'pasteur'/
        # 'ncbi') - two entirely different code spaces that can never
        # match, so this always fell through to 'Error' regardless of
        # the actual server. No current caller (found while fixing the
        # BLASTSERVERS mislabeling right above - both NCBI and Pasteur
        # runs displayed as "Pasteur" anywhere that used it), but worth
        # fixing alongside rather than leaving broken next to the fix.
        for (code, desc) in self.BLASTSERVERS:
            if self.server == code:
                return desc
        return 'Error'

    def finished(self):
        '''
        Nor running anymore (success or error)
        '''
        return self.status == self.FINISHED or self.status == self.ERROR

    def soft_delete(self):
        self.deleted = True
        self.blastsubject_set.all().delete()
        self.save()

    def is_prot(self):
        return self.blastprog in ['blastp', 'blastx', 'tblastn', 'tblastx']

    def build_nj_tree(self):
        dm = self.distance_matrix()
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        treeio = StringIO()
        Phylo.write(tree, treeio, 'newick')
        treestr = treeio.getvalue()
        treeio.close()
        return treestr

    def distance_matrix(self):
        names = []
        matrix = []
        seqs = []
        names.append(str(self.query_id))
        seqs.append("".join(self.query_seq))
        for s in self.blastsubject_set.all():
            id = s.subject_id
            seq = s.subject_seq
            names.append(str(id))
            seqs.append("".join(seq))
        for i in range(0, len(names)):
            matrix.append([])
            for j in range(0, i+1):
                d = 0.0
                if i != j:
                    if self.is_prot():
                        d = prot_dist(seqs[i], seqs[j])
                    else:
                        d = nucl_dist(seqs[i], seqs[j])
                matrix[i].append(d)
        return DistanceMatrix(names=names, matrix=matrix)

    @staticmethod
    def blast_servers():
        context = dict()
        for server in settings.BLASTS:
            if settings.BLASTS.get(server).get('activated'):
                name = settings.BLASTS.get(server).get('name')
                context.update({server : name})
        return context

    @staticmethod
    def blast_progs(server):
        context = dict()
        blast = settings.BLASTS.get(server)
        if blast is not None and blast.get('activated'):
            progs = blast.get('progs')
            for prog in progs:
                context.update({prog : progs.get(prog).get('name')})
        return context

    @staticmethod
    def blast_dbs(server, prog):
        context = dict()
        blast = settings.BLASTS.get(server)
        if blast is not None and blast.get('activated'):
            blastprog = blast.get('progs').get(prog)
            if blastprog is not None:
                blastdbs = blastprog.get('dbs')
                for db in blastdbs:
                    context.update({db : blastdbs.get(db)})
        return context

    @staticmethod
    def blast_type(server, prog):
        blast = settings.BLASTS.get(server)
        if blast is not None and blast.get('activated'):
            blastprog = blast.get('progs').get(prog)
            if blastprog is not None:
                type = blastprog.get('type')
                return type
        return None

    @staticmethod
    def blast_inputtype(server, prog):
        """
        returns the input type of the given program
        Ex: blast_inputtype(BlastRun.PASTEUR, "blastx") returns "nt"
        """
        blast = settings.BLASTS.get(server)
        if blast is not None and blast.get('activated'):
            blastprog = blast.get('progs').get(prog)
            if blastprog is not None:
                inputtype = blastprog.get('input')
                return inputtype
        return None
    
    @staticmethod
    def blast_example(server, prog):
        context = []
        blast = settings.BLASTS.get(server)
        if blast is not None and blast.get('activated'):
            blastprog = blast.get('progs').get(prog)
            if blastprog is not None:
                testdata = blastprog.get('test_data')
                if testdata is not None:
                    filestr = os.path.join(settings.TESTDATA_DIR,testdata)
                    if os.path.exists(filestr):
                        with open(filestr, "r") as f:
                            context.append(f.read())
        return context

class BlastSubject(models.Model):
    subject_id = models.CharField(max_length=1000)
    subject_seq = models.TextField()
    subject_fullseq = models.TextField()
    blastrun = models.ForeignKey(BlastRun, on_delete=models.CASCADE)

    def format_sequence(self):
        unalignseq = self.subject_seq.replace("-", "")
        return ('\n'.join(textwrap.wrap(unalignseq, 60))).rstrip()

    def format_fullsequence(self):
        unalignseq = self.subject_fullseq.replace("-", "")
        return ('\n'.join(textwrap.wrap(unalignseq, 60))).rstrip()

    def fullsequence(self):
        unalignseq = self.subject_fullseq.replace("-", "")
        return unalignseq.rstrip()
    
