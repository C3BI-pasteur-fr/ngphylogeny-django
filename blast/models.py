# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

from datetime import datetime
import uuid
import textwrap

# Create your models here.
class BlastRun(models.Model):
    PENDING   = 'P'
    RUNNING   = 'R'
    FINISHED  = 'F'
    ERROR = 'E'
    RUNSTATUS = (
        (PENDING, 'Pending'),
        (RUNNING, 'Running'),
        (FINISHED,'Finished'),
        (ERROR,   'Error')
    )
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    date = models.DateTimeField(default=datetime.now, blank=True)
    query_id = models.CharField(null=True, max_length=1000)
    query_seq = models.TextField(null=True)
    evalue = models.FloatField(default=0.00001)
    database = models.CharField(max_length=100,default='swissprot')
    blastprog = models.CharField(max_length=100, default='blastp')
    status = models.CharField(max_length=1,default=PENDING, choices=RUNSTATUS)
    message = models.TextField(null=True)
    deleted = models.BooleanField(default=False)

    def format_sequence(self):
        return ('\n'.join(textwrap.wrap(self.query_seq, 60))).rstrip()

    def to_fasta(self):
        fasta= ">%s\n" % self.query_id
        fasta+= "%s\n" % self.format_sequence()
        for s in self.blastsubject_set.all():
            fasta+= ">%s\n" % s.subject_id
            fasta+= "%s\n" % s.format_sequence()
        return fasta

    def status_str(self):
        for (code,desc) in self.RUNSTATUS:
            if self.status == code:
                return desc
        return 'Error'

    def finished(self):
        '''
        Nor running anymore (success or error)
        '''
        return self.status== self.FINISHED or self.status == self.ERROR

    def soft_delete(self):
        self.deleted = True
        self.blastsubject_set.all().delete()
        self.save()

    
class BlastSubject(models.Model):
    subject_id = models.CharField(max_length=1000)
    subject_seq = models.TextField()
    blastrun = models.ForeignKey(BlastRun, on_delete=models.CASCADE)

    def format_sequence(self):
        return ('\n'.join(textwrap.wrap(self.subject_seq, 60))).rstrip()
