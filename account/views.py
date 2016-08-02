from django.shortcuts import render
from django.conf import settings

# Create your views here.

def galaxy_connection_error_view(request):

    error_msg = "Connection attempts with the Galaxy server at "+settings.GALAXY_SERVER_URL+" failed. "\
                "Please check Galaxy server is properly configured and online" \

    error = {"message": error_msg}

    return render(request, 'error.html', {'error': error})

