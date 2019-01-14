# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os

from django.db import models
from django.conf import settings

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo

from datetime import datetime
import uuid
import textwrap
import logging
import StringIO

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
        (NCBI, 'Pasteur'),
    )
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    email = models.CharField(null=True, max_length=100)
    date = models.DateTimeField(default=datetime.now, blank=True)
    query_id = models.CharField(null=True, max_length=1000)
    query_seq = models.TextField(null=True)
    evalue = models.FloatField(default=0.00001)
    coverage = models.FloatField(default=0.8)
    maxseqs = models.PositiveIntegerField(default=10)
    database = models.CharField(max_length=100, default='swissprot')
    blastprog = models.CharField(max_length=100, default='blastp')
    history = models.CharField(max_length=20) # If pasteur blast: galaxy history id
    history_fileid= models.CharField(max_length=20) # If pasteur blast: output file galaxy id 
    status = models.CharField(max_length=1, default=PENDING, choices=RUNSTATUS)
    server = models.CharField(max_length=50, default=NCBI, choices=BLASTSERVERS)
    message = models.TextField(null=True)
    deleted = models.BooleanField(default=False)
    tree = models.TextField(null=True)

    def format_sequence(self):
        return ('\n'.join(textwrap.wrap(self.query_seq, 60))).rstrip()

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
        for (code, desc) in self.BLASTSERVERS:
            if self.status == code:
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
        treeio = StringIO.StringIO()
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
