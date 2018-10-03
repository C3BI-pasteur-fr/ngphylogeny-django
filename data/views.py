import urllib
import json

try:
    # Python 3:
    from urllib.parse import urlparse

except ImportError:
    # Python 2:
    import urlparse

import tempfile
import requests
from django.http import StreamingHttpResponse
from django.http import JsonResponse
from django.shortcuts import render, redirect
from django.urls import reverse_lazy
from django.utils.decorators import method_decorator
from django.views.generic.edit import FormView
from .models import ExampleFile
from .forms import UploadForm
from galaxy.decorator import connection_galaxy
from workspace.views import get_or_create_history
from blast.models import BlastRun

class UploadMixin(object):
    def upload_content(self, content, history_id=None, name="pasted_data"):
        """
            send content into galaxy history: return galaxy response
        """
        if history_id:
            self.history_id = history_id
        else:
            self.history_id = get_or_create_history(self.request)

        return self.request.galaxy.tools.paste_content(content=content, file_name=name,
                                                       history_id=self.history_id)

    def upload_file(self, file, history_id=None):
        """
            upload file into galaxy history: return galaxy response
        """

        tmpfile = tempfile.NamedTemporaryFile()
        for chunk in file.chunks():
            tmpfile.write(chunk)
        tmpfile.flush()

        if history_id:
            self.history_id = history_id
        else:
            self.history_id = get_or_create_history(self.request)

        return self.request.galaxy.tools.upload_file(path=tmpfile.name, file_name=file.name, history_id=self.history_id)


@method_decorator(connection_galaxy, name="dispatch")
class UploadView(UploadMixin, FormView):
    """
       Upload file into Galaxy Server
    """

    template_name = 'upload_form.html'
    form_class = UploadForm
    success_url = reverse_lazy("home")

    # Add blast runs id to kwargs
    # Used by UploadForm to get blastruns ids
    def get_form_kwargs(self):
        kwargs = super(UploadView, self).get_form_kwargs()
        if self.request.session.get('blastruns'):
            blastruns = []
            for r in self.request.session['blastruns']:
                b = BlastRun.objects.get(pk=r)
                if b != None and b.status==b.FINISHED :
                    blastruns.append(r)
            kwargs['blastruns'] = blastruns
        return kwargs
    
    def form_valid(self, form):
        if form.cleaned_data.get('input_file') is not None:
            myfile = form.cleaned_data.get('input_file')
            outputs = self.upload_file(myfile)
        elif form.cleaned_data.get('pasted_text') is not None:
            p_content = form.cleaned_data['pasted_text']
            outputs = self.upload_content(p_content)
        elif form.cleaned_data.get('blast_run') is not None:
            # We treat blastrun id as fasta content
            b = BlastRun.objects.get(pk=form.cleaned_data.get('blast_run'))
            p_content = b.to_fasta()
            outputs = self.upload_content(p_content)

        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': self.history_id}, )

        return super(UploadView, self).form_valid()
    

@connection_galaxy
def download_file(request, file_id):
    """permet a l'utilisateur de telecharger le fichier grace a l'api"""
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    name = "error"
    if isinstance(data, dict):
        dlurl = data.get('download_url')
        name = data.get('name')
        if not name:
            name = "download"
        if dlurl:
            url = urlparse.urljoin(gi.base_url, dlurl)
            response = urllib.urlopen(url)
            stream_response = StreamingHttpResponse(response.read())
            stream_response['Content-Disposition'] = 'attachment; filename=' + name
        else:
            stream_response = StreamingHttpResponse("No file download URL corresponds to the given dataset id " + file_id)

    else:
        stream_response = StreamingHttpResponse(data)
    return stream_response

@connection_galaxy
def display_raw(request, file_id):
    """Display file content to the web browser """
    stream_response = download_file(request, file_id)
    del stream_response['Content-Disposition']
    return stream_response


@connection_galaxy
def display_file(request, file_id):
    """Display file content to the web browser """
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    if isinstance(data, dict):
        historyid = data.get('history_id')
        if historyid:
            if request.is_ajax():
                return display_raw(request, file_id)
            else:
                return render(request, 'display.html', {'history_id': historyid})
    return render(request, 'error.html', {'errortitle': 'Error querying galaxy', 'errormessage': data})

@connection_galaxy
def display_params(request, file_id):
    """Display tool run parameters """
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    if isinstance(data, dict):
        job_id = data.get('creating_job')
        if job_id:
            if request.is_ajax():
                job = gi.jobs.show_job(job_id)
                return JsonResponse(job)
            else:
                return render(request, 'display_params.html', {'history_id': data.get('history_id')})
    return render(request, 'error.html', {'errortitle': 'Error querying galaxy', 'errormessage': data})

@connection_galaxy
def display_msa(request, file_id):
    """Display multiple alignment file content to the web browser """
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    if isinstance(data, dict):
        historyid = data.get('history_id')
        if historyid:
            if request.is_ajax():
                return display_raw(request, file_id)
            else:
                return render(request, 'msaviz/msa.html', {'history_id': historyid})
    return render(request, 'error.html', {'errortitle': 'Error querying galaxy', 'errormessage': data})


@connection_galaxy
def tree_visualization(request, file_id):

    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    if isinstance(data, dict):
        dlurl = data.get('download_url')
        historyid = data.get('history_id')
        if dlurl and historyid:
            url = urlparse.urljoin(gi.base_url, dlurl)
            response = urllib.urlopen(url)
            return render(request,
                          template_name='treeviz/tree.html',
                          context={'newick_tree': response.read(),
                                   'history_id': historyid})
    return render(request, 'error.html', {'errortitle': 'Error querying galaxy', 'errormessage': data})


@connection_galaxy
def export_to_itol(request, file_id):
    # retrieve newick from galaxy server
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)

    if isinstance(data, dict):
        dlurl = data.get('download_url')
        if dlurl:
            url = urlparse.urljoin(gi.base_url, dlurl)
            response = urllib.urlopen(url)
            tmpfile = tempfile.NamedTemporaryFile()
            tmpfile.write(response.read())
            tmpfile.flush()
            # send file to itol server
            url_itol = 'https://itol.embl.de/upload.cgi'
            payload = {'tname': "", 'tfile': open(tmpfile.name, 'rb'), }
            r = requests.post(url_itol, files=payload)
            return redirect(r.url)
    return render(request, 'error.html', {'errortitle': 'Error querying galaxy', 'errormessage': data})

def get_example(request):
    """return example file according to ext_file found """

    tool_id = ''
    if request.POST:
        # ext_file = ast.literal_eval((request.POST.get('ext_file', [])))
        tool_id = request.POST.get('tool_id', '')
        tool_input = request.POST.get('input_name', '')
    if tool_id and tool_input:
        example_file = ExampleFile.objects.filter(toolinputdata__tool_id=tool_id,
                                                  toolinputdata__name=tool_input).first()
    else:
        ext_file = request.GET.get('ext_file')
        example_file = ExampleFile.objects.filter(ext=ext_file).first()

    if example_file:
        stream_response = StreamingHttpResponse(open(example_file.upload.path).read())
        return stream_response

    return StreamingHttpResponse()

@connection_galaxy
def add_file_to_session(request, file_id):
    """
    Adds a Galaxy result file to the session,
    in order to use it in an other workflow
    then redirects to the history detail page
    """
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    if isinstance(data, dict):
        if not request.session.get('files'):
            request.session['files']={}
        fdict = request.session['files']
        if file_id not in fdict:
            print "data:"
            print json.dumps(data)
            fdict[file_id]={'id': file_id, 'ext' : data.get('file_ext'), 'history' : data.get('history_id'), 'name': data.get('name')}
        return redirect('history_detail', history_id=data.get('history_id'))
    return render(request, 'error.html', {'errortitle': 'Error while adding file to session', 'errormessage': 'File id does not exist'})

@connection_galaxy
def remove_file_from_session(request, file_id):
    """
    Remove a Galaxy result file from the session
    """
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    if isinstance(data, dict):
        if not request.session.get('files'):
            request.session.files={}
        fdict = request.session['files']
        if file_id in fdict:
            del fdict[file_id]
        return redirect('history_detail', history_id=data.get('history_id'))
    return render(request, 'error.html', {'errortitle': 'Error while remove file from session', 'errormessage': 'File id does not exist'})
