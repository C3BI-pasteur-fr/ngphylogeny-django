"""NGPhylogeny_fr URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
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
from django.conf.urls import url
from django.contrib import admin
from django.views.generic import TemplateView, RedirectView

from django.conf.urls import include
from data import urls as data_urls
from workspace import urls as workspace_urls
from tools import urls as tool_urls
from workflows import urls as workflows_urls
from django.conf import settings


urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'galaxy/$', RedirectView.as_view(url=settings.GALAXY_SERVER_URL), name="galaxy_url"),
    url(r'galaxy/(?P<id>[\w-]+)$', RedirectView.as_view(url=settings.GALAXY_SERVER_URL+'dataset/errors?id=%(id)s'), name="galaxy_error_url"),
    url(r'^documentation$', TemplateView.as_view(template_name="documentation.html"), name="documentation"),
    url(r'^tools/', include(tool_urls)),
    url(r'^data/', include(data_urls)),
    url(r'^workflows/', include(workflows_urls)),
    url(r'^workspace/', include(workspace_urls)),
    url(r'^$', TemplateView.as_view(template_name="home.html"), name="home"),
]
