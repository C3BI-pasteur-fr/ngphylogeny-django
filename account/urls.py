
from django.conf.urls import url
from django.contrib.auth.views import login, logout
from django.contrib.auth.decorators import login_required
from django.views.generic import TemplateView


from views import *

urlpatterns = [
    url(r'^login$', login, {'template_name': "account/login.html" }, name='login'),
    url(r'logout$', logout, {'next_page': "/"}, name='logout'),
    url(r'^$',login_required( GalaxyUserUpdateApiKey.as_view() ), name='account'),
    #url(r'^create_account/$', create_account),
    #url(r'^success/$', TemplateView.as_view(template_name='success.html'))
]