from django.shortcuts import render
from django.views.generic import TemplateView


# Create your views here.

def galaxy_connection_error_view(request):

    error_msg = "Connection attempts with the Galaxy server at "+request.galaxy_server+" failed. " \
                "Please check Galaxy server is properly configured and online"

    error = {"message": error_msg}

    return render(request, 'error.html', {'error': error})


class HttpResponseError(TemplateView):

    template_name = "error.html"