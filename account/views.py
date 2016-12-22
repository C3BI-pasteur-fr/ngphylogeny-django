from django.conf import settings
from django.views.generic import UpdateView
from django.contrib.messages.views import SuccessMessageMixin
from django.shortcuts import get_object_or_404
from django.shortcuts import render
from models import GalaxyUser

# Create your views here.

def galaxy_connection_error_view(request):

    error_msg = "Connection attempts with the Galaxy server at "+settings.GALAXY_SERVER_URL+" failed. " \
                "Please check Galaxy server is properly configured and online"

    error = {"message": error_msg}

    return render(request, 'error.html', {'error': error})


class GalaxyUserUpdateApiKey(SuccessMessageMixin, UpdateView):


    model = GalaxyUser
    template_name = "account/user_info.html"
    fields = ['api_key', ]
    success_url = '/account'
    success_message = "Your Api Key was updated successfully"

    def get_object(self, queryset=None):
        return get_object_or_404(GalaxyUser, user=self.request.user)