# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.views.generic.edit import FormView
from django.views.generic import DetailView, DeleteView, TemplateView, View
from django.urls import reverse_lazy
from django.http import HttpResponseRedirect
from django.conf import settings
from django.http import HttpResponse

from .forms import BlastForm
from tasks import launch_ncbi_blast, launch_pasteur_blast, build_tree
from .models import BlastRun, BlastSubject

import json

def add_blast_id_to_session(request, blastrunid):
    if not request.session.get('blastruns'):
        request.session['blastruns'] = []
    rundict = request.session['blastruns']
    if blastrunid not in rundict:
        rundict.append(blastrunid)
        request.session["blastruns"] = rundict

def available_blasts_servers(request):
    """
    Ajax: return possible blasts servers : {id:name}
    """
    context = BlastRun.blast_servers()
    return HttpResponse(json.dumps(context), content_type='application/json')

def available_blasts_progs(request, server):
    """
    Ajax: return possible blasts progs : {id:name}
    """
    context = BlastRun.blast_progs(server)
    return HttpResponse(json.dumps(context), content_type='application/json')

def available_blasts_dbs(request, server, prog):
    """
    Ajax: return possible blasts databases : {id, name}
    """
    
    context = BlastRun.blast_dbs(server, prog)
    return HttpResponse(json.dumps(context), content_type='application/json')

def blast_example(request, server, prog):
    """
    Ajax: return fasta sequence (prot or nt depending on prog)
    """
    context = BlastRun.blast_example(server, prog)
    return HttpResponse(json.dumps(context), content_type='application/json')


class BlastView(FormView, TemplateView):
    template_name = "blast/blast.html"
    form_class = BlastForm
    success_url = 'blast_view'
    runid = ""
    model = BlastRun

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        seq = form.cleaned_data.get('sequence')
        prog = form.cleaned_data.get('program')
        server = form.cleaned_data.get('server')
        db = form.cleaned_data.get('database')
        maxseqs = form.cleaned_data.get('hitslimit')
        evalue = form.cleaned_data.get('evalue')
        coverage = form.cleaned_data.get('coverage')
        email = form.cleaned_data.get('email')
        b = BlastRun()
        b.query_id = ""
        b.query_seq = ""
        b.server = server
        b.email = email
        b.save()
        self.runid = b.id
        if server == 'pasteur':
            launch_pasteur_blast.delay(b.id, seq, prog, db, evalue, coverage, maxseqs)
        else:
            launch_ncbi_blast.delay(b.id, seq, prog, db, evalue, coverage, maxseqs)
        add_blast_id_to_session(self.request, str(b.id))
        return super(BlastView, self).form_valid(form)

    def get_success_url(self):
        return reverse_lazy(self.success_url, kwargs={'pk': str(self.runid)})

    # We take only blast runs whose id are in the session
    def get_context_data(self, **kwargs):
        context = super(BlastView, self).get_context_data(**kwargs)
        sessionids = []
        if self.request.session.get('blastruns'):
            sessionids = self.request.session.get('blastruns')
        queryset = self.model.objects.filter(pk__in=sessionids, deleted=False).order_by('-date')
        context["blast_runs"] = queryset
        return context


class BlastRunView(DetailView):
    model = BlastRun
    template_name = 'blast/result.html'

    def get_object(self):
        o = super(BlastRunView, self).get_object()
        add_blast_id_to_session(self.request, str(o.id))
        return o


class DeleteBlastSubjectView(DeleteView):
    model = BlastSubject
    success_url = 'blast_view'

    def get_success_url(self):
        """ Go to the BlastRun page"""
        blastrun = self.object.blastrun
        """ And relaunch the pseudo tree building """
        blastrun.tree = ""
        blastrun.status = blastrun.RUNNING
        blastrun.save()
        build_tree.delay(blastrun.id)
        return reverse_lazy(self.success_url, kwargs={'pk': blastrun.id})

    def get(self, request, *args, **kwargs):
        """ No confirmation template """
        return self.post(request, *args, **kwargs)


class DeleteBlastSequences(View):

    def post(self, request, *args, **kwargs):
        id = self.kwargs['pk']
        blastrun = BlastRun.objects.get(pk=id)

        for seqid in request.POST.getlist('todelete[]'):
            BlastSubject.objects.filter(blastrun__id=id, subject_id=seqid).delete()
            
        blastrun.tree = ""
        blastrun.status = blastrun.RUNNING
        blastrun.save()
        build_tree.delay(id)
        
        return HttpResponseRedirect(reverse_lazy('blast_view', kwargs=kwargs))

    def get(self, request, *args, **kwargs):
        return HttpResponseRedirect(reverse_lazy('blast_view', kwargs=kwargs))
    
class DeleteBlastRunView(DeleteView):
    model = BlastRun
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
    model = BlastRun
    template_name = 'blast/result_fasta.txt'
    content_type = 'text/plain'

    def get_object(self):
        o = super(BlastRunFasta, self).get_object()
        add_blast_id_to_session(self.request, str(o.id))
        return o

