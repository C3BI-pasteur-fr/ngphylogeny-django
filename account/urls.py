from django.urls import re_path
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import LoginView, LogoutView

from galaxy.views import UpdateApiKey

urlpatterns = [
    re_path(r'^login$', LoginView.as_view(template_name="account/login.html"), name='login'),
    re_path(r'logout$', LogoutView.as_view(next_page="/"), name='logout'),
    re_path(r'^$', login_required(UpdateApiKey.as_view()), name='account'),
    #url(r'^create_account/$', create_account),
    #url(r'^success/$', TemplateView.as_view(template_name='success.html'))
]