import tempfile
import urllib
import urlparse
from django.shortcuts import render
from django.core.urlresolvers import reverse_lazy
from django.http import StreamingHttpResponse
from django.views.generic import FormView
from forms import UploadForm
from django.utils.decorators import method_decorator
from account.decorator import connection_galaxy


@method_decorator(connection_galaxy, name="dispatch")
class UploadView(FormView):
    """
        Upload file into Galaxy Server
    """

    template_name = 'upload_form.html'
    form_class = UploadForm
    success_url = reverse_lazy("home")

    # TODO Page d'erreur

    def form_valid(self, form):

        myfile = form.cleaned_data['file']
        tmpfile = tempfile.NamedTemporaryFile()
        for chunk in myfile.chunks():
            tmpfile.write(chunk)
        tmpfile.flush()

        gi = self.request.galaxy

        try:
            history_id = gi.histories.get_most_recently_used_history().get('id')

        except:

            # Create a new galaxy history and delete older if the user is not authenticated
            history_id = gi.histories.create_history().get('id')

        hist = gi.tools.upload_file(tmpfile.name, history_id, file_name=myfile.name)
        self.success_url = reverse_lazy("history_detail", kwargs={'history_id': hist.get('id', history_id)}, )

        # TODO
        # supprime l'history en cas d'erreur
        # gi.histories.delete_history(history_id)

        return super(UploadView, self).form_valid(form)


@connection_galaxy
def download_file(request, file_id):

    """permet a l'utilisateur de telecharger le fichier grace a l'api"""
    gi = request.galaxy
    data = gi.datasets.show_dataset(dataset_id=file_id)
    url = urlparse.urljoin(gi.base_url, data['download_url'])

    response = urllib.urlopen(url)
    stream_response = StreamingHttpResponse(response.read())
    stream_response['Content-Disposition'] = 'attachment; filename='+data["name"]

    # TODO test if bigDATA
    return stream_response

@connection_galaxy
def display_file(request, file_id):
    """Affiche le fichier dans le navigateur"""
    if request.is_ajax():
        stream_response = download_file(request, file_id)
        del stream_response['Content-Disposition']
        return stream_response
    else:
        return render(request, 'display.html')

@connection_galaxy
def export_file(request):
    pass
