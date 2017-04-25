from django.contrib.messages.views import SuccessMessageMixin
from django.shortcuts import render
from django.views.generic import UpdateView, TemplateView

from models import GalaxyUser, Server


# Create your views here.

def galaxy_connection_error_view(request):

    error_msg = "Connection attempts with the Galaxy server at "+request.galaxy_server+" failed. " \
                "Please check Galaxy server is properly configured and online"

    error = {"message": error_msg}

    return render(request, 'error.html', {'error': error})


class HttpResponseError(TemplateView):

    template_name = "error.html"



class UpdateApiKey(SuccessMessageMixin, UpdateView):
    """

    """
    queryset = GalaxyUser.objects.filter(anonymous=True)
    template_name = "account/user_info.html"
    fields = ['api_key', ]
    success_url = '/account'
    success_message = "Your Api Key was updated successfully"

    def get_object(self, queryset=queryset):

        galaxy_server = Server.objects.get(current=True)
        return GalaxyUser.objects.get_or_create(user=self.request.user, galaxy_server=galaxy_server)[0]