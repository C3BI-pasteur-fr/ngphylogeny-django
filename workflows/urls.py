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
from views import *
from django.views.generic import TemplateView

urlpatterns = [

    url(r'^oneclick_mode$', WorkflowOneClickListView.as_view(), name="workflow_oneclick_list"),
    url(r'^oneclick_mode/(?P<slug>[\w-]+)/$', WorkflowOneClickView.as_view(), name="workflow_oneclick_form"),
    #url(r'^oneclick_mode$', workflows_oneclick_mode_build, name="workflows_oneclick"),

    url(r'advanced_mode$', workflows_advanced_mode_build, name="workflows_advanced"),
    url(r'^alacarte_mode$', workflows_alacarte_mode_build, name="workflows_alacarte"),
    #url(r'^(?P<slug>[\w-]+)/$', LaunchGalaxyWorkflowView.as_view(), name="workflows"),

    #url(r'^$', workflow_form, name="workflow_form"),
]


