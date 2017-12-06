import urllib

try:
    # Python 3:
    from urllib.parse import urlparse

except ImportError:
    # Python 2:
    import urlparse

import tempfile
import requests
from django.http import StreamingHttpResponse
from django.shortcuts import render, redirect
from django.urls import reverse_lazy
from django.utils.decorators import method_decorator
from django.views.generic.edit import FormView
from .models import ExampleFile
from .forms import UploadForm, PastedContentForm
from galaxy.decorator import connection_galaxy
from workspace.views import get_or_create_history


class UploadMixin(object):
    def upload_content(self, content, history_id=None):
        """
            send content into galaxy history: return galaxy response
        """
        if history_id:
            self.history_id = history_id
        else:
            self.history_id = get_or_create_history(self.request)

        return self.request.galaxy.tools.paste_content(content=content, file_name="pasted_data",
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

    def form_valid(self, form):
        myfile = form.cleaned_data['file']
        outputs = self.upload_file(myfile)
        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': self.history_id}, )

        return super(UploadView, self).form_valid()


@method_decorator(connection_galaxy, name="dispatch")
class ImportPastedContentView(UploadMixin, FormView):
    """
       Import user pasted content into Galaxy Server
    """

    form_class = PastedContentForm
    success_url = reverse_lazy("home")

    def form_valid(self, form):
        p_content = form.cleaned_data['pasted_text']
        outputs = self.upload_content(p_content)
        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': self.history_id}, )

        return super(ImportPastedContentView, self).form_valid()


@connection_galaxy
def download_file(request, file_id):
    """permet a l'utilisateur de telecharger le fichier grace a l'api"""
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    url = urlparse.urljoin(gi.base_url, data['download_url'])

    response = urllib.urlopen(url)
    stream_response = StreamingHttpResponse(response.read())
    stream_response['Content-Disposition'] = 'attachment; filename=' + data["name"]

    # TODO test if bigDATA
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

    if request.is_ajax():
        return display_raw(request, file_id)
    else:

        return render(request, 'display.html', {'history_id': data.get('history_id')})


@connection_galaxy
def display_msa(request, file_id):
    """Display multiple alignment file content to the web browser """
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)

    if request.is_ajax():
        return display_raw(request, file_id)
    else:
        return render(request, 'msaviz/msa.html', {'history_id': data.get('history_id')})


@connection_galaxy
def tree_visualization(request, file_id):

    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)

    url = urlparse.urljoin(gi.base_url, data['download_url'])
    response = urllib.urlopen(url)

    return render(request, template_name='treeviz/tree.html', context={'newick_tree': response.read(),
                                                                       'history_id': data.get('history_id')})


@connection_galaxy
def export_to_itol(request, file_id):
    # retrieve newick from galaxy server
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)

    url = urlparse.urljoin(gi.base_url, data['download_url'])
    response = urllib.urlopen(url)

    tmpfile = tempfile.NamedTemporaryFile()
    tmpfile.write(response.read())
    tmpfile.flush()

    # send file to itol server
    url_itol = 'http://itol.embl.de/upload.cgi'
    payload = {'tname': "", 'tfile': open(tmpfile.name, 'rb'), }
    r = requests.post(url_itol, files=payload)

    return redirect(r.url)


@connection_galaxy
def export_file(request):
    pass


def get_example(request, ext_file=""):
    """return example file according to ext_file found """

    tool = ''
    if request.POST:
        # ext_file = ast.literal_eval((request.POST.get('ext_file', [])))
        tool = request.POST.get('tool_id', '')
    if tool:
        example_file = ExampleFile.objects.filter(toolinputdata__tool_id=tool).first()
    else:
        ext_file = [str(ext).strip() for ext in ext_file]
        example_file = ExampleFile.objects.filter(ext__in=ext_file).first()

    if example_file:
        stream_response = StreamingHttpResponse(open(example_file.upload.path).read())
        return stream_response

    return StreamingHttpResponse()
