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
from django.urls import include, re_path
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
    re_path(r'^admin/', admin.site.urls, name='admin'),  # Django ADMIN URLS
    re_path(r'^about$', TemplateView.as_view(template_name="about.html"), name="about"),
    re_path(r'^about/', include(surveys_urls)),
    re_path(r'^documentation$', TemplateView.as_view(template_name="documentation.html"), name="documentation"),
    re_path(r'^analysis$', TemplateView.as_view(template_name="phylogeny_analysis_choices.html"), name="analysis_list"),
    re_path(r'^status$', TemplateView.as_view(template_name="status.json",content_type='application/json'), name="status"),
    re_path(r'^galaxy/', include(galaxy_urls)),
    re_path(r'^account/', include(account_urls)),
    re_path(r'^tools/', include(tool_urls)),
    re_path(r'^data/', include(data_urls)),
    re_path(r'^workflows/', include(workflows_urls)),
    re_path(r'^workspace/', include(workspace_urls)),
    re_path(r'^blast/',include(blast_urls)),
    re_path(r'^$', TemplateView.as_view(template_name="home.html"), name="home"),
    # url(r'.*', TemplateView.as_view(template_name="maintenance.html"), name="home"),
]
