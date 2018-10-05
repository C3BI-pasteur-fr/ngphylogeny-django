# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.views.generic.edit import FormView
from django.views.generic import DetailView, DeleteView, ListView
from django.views.generic.detail import BaseDetailView
from django.urls import reverse_lazy
from  django.http import HttpResponseRedirect

from .forms import BlastForm
from tasks import launchblast
from .models import BlastRun, BlastSubject
from galaxy.decorator import connection_galaxy
from workspace.views import create_history

def add_blast_id_to_session(request, blastrunid):
    if not request.session.get('blastruns'):
        request.session['blastruns']=[]
    rundict = request.session['blastruns']
    if blastrunid not in rundict:
        rundict.append(blastrunid)
        request.session["blastruns"] = rundict

    
class BlastView(FormView, ListView):
    template_name = "blast/blast.html"
    form_class = BlastForm
    success_url = 'blast_view'
    runid=""
    model=BlastRun

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        seq = form.cleaned_data.get('sequence')
        prog= form.cleaned_data.get('program')
        db= form.cleaned_data.get('database')
        evalue= form.cleaned_data.get('evalue')
        coverage= form.cleaned_data.get('coverage')
        email= form.cleaned_data.get('email')
        b = BlastRun()
        b.query_id=""
        b.query_sequence=""
        b.email = email
        b.save()
        self.runid=b.id
        launchblast.delay(b.id, seq, prog, db, evalue, coverage)
        add_blast_id_to_session(self.request, str(b.id))
        #form.send_email()
        return super(BlastView, self).form_valid(form)

    def get_success_url(self):
        return reverse_lazy( self.success_url , kwargs={'pk': str(self.runid)})

    # We take only blast runs whose id are in the session
    def get_queryset(self):
        sessionids= []
        if self.request.session.get('blastruns'):
            sessionids=self.request.session.get('blastruns')
        queryset = self.model.objects.filter(pk__in=sessionids, deleted=False).order_by('-date')
        return queryset


class BlastRunView(DetailView):
    model=BlastRun
    template_name = 'blast/result.html'

    def get_object(self):
        o = super(BlastRunView,self).get_object()
        add_blast_id_to_session(self.request,str(o.id))
        return o
        

class DeleteBlastSubjectView(DeleteView):
    model=BlastSubject
    success_url='blast_view'
    
    def get_success_url(self):
        """ Go to the BlastRun page"""
        blastrun = self.object.blastrun
        return reverse_lazy( self.success_url, kwargs={'pk': blastrun.id})
    
    def get(self, request, *args, **kwargs):
        """ No confirmation template """
        return self.post(request, *args, **kwargs)

class DeleteBlastRunView(DeleteView):
    model=BlastRun
    success_url = reverse_lazy('blast_form')

    def get(self, request, *args, **kwargs):
        """ No confirmation template """
        return self.post(request, *args, **kwargs)

    # Overrides the default delete method to prevent actual
    # deletion, but instead mark it as deleted and delete the
    # subject sequence results
    def delete(self, request, *args, **kwargs):
        """
        Calls the delete() method on the fetched object and then
        redirects to the success URL.
        """
        self.get_object().soft_delete()
        return HttpResponseRedirect(self.success_url)
    

class BlastRunFasta(DetailView):
    model=BlastRun
    template_name = 'blast/result_fasta.txt'
    content_type='text/plain'
    
    def get_object(self):
        o = super(BlastRunFasta,self).get_object()
        add_blast_id_to_session(self.request,str(o.id))
        return o

