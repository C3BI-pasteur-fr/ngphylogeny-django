"""NGPhylogeny_fr URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import include
from django.conf.urls import url
from django.contrib import admin
from django.views.generic import TemplateView

from account import urls as account_urls
from data import urls as data_urls
from galaxy import urls as galaxy_urls
from surveys import urls as surveys_urls
from tools import urls as tool_urls
from workflows import urls as workflows_urls
from workspace import urls as workspace_urls
from blast import urls as blast_urls

urlpatterns = [
    url(r'^admin/', admin.site.urls, name='admin'),  # Django ADMIN URLS
    url(r'^about$', TemplateView.as_view(template_name="about.html"), name="about"),
    url(r'^about/', include(surveys_urls)),
    url(r'^documentation$', TemplateView.as_view(template_name="documentation.html"), name="documentation"),
    url(r'^analysis$', TemplateView.as_view(template_name="phylogeny_analysis_choices.html"), name="analysis_list"),
    url(r'^status$', TemplateView.as_view(template_name="status.json",content_type='application/json'), name="status"),
    url(r'^galaxy/', include(galaxy_urls)),
    url(r'^account/', include(account_urls)),
    url(r'^tools/', include(tool_urls)),
    url(r'^data/', include(data_urls)),
    url(r'^workflows/', include(workflows_urls)),
    url(r'^workspace/', include(workspace_urls)),
    url(r'^blast/',include(blast_urls)),
    url(r'^$', TemplateView.as_view(template_name="home.html"), name="home"),
    # url(r'.*', TemplateView.as_view(template_name="maintenance.html"), name="home"),
]
