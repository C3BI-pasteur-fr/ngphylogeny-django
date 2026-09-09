from django.conf.urls import url
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import LoginView, LogoutView

from galaxy.views import UpdateApiKey

urlpatterns = [
    url(r'^login$', LoginView.as_view(template_name="account/login.html"), name='login'),
    url(r'logout$', LogoutView.as_view(next_page="/"), name='logout'),
    url(r'^$', login_required(UpdateApiKey.as_view()), name='account'),
    #url(r'^create_account/$', create_account),
    #url(r'^success/$', TemplateView.as_view(template_name='success.html'))
]