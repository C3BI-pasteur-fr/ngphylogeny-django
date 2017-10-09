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

from views.oneclick import WorkflowStartedView, WorkflowOneClickListView, WorkflowOneClickFormView
from views.wkadvanced import WorkflowAdvancedListView, WorkflowAdvancedRedirectView
from views.wkmaker import WorkflowsMarkerRedirectView, workflows_alacarte_build

urlpatterns = [

    url(r'^quickstart$', WorkflowStartedView.as_view(), name="get_started_workflow"),
    url(r'^oneclick/$', WorkflowOneClickListView.as_view(), name="workflow_oneclick_list"),
    url(r'^oneclick/(?P<slug>[\w-]+)$', WorkflowOneClickFormView.as_view(), name="workflow_oneclick_form"),
    url(r'^advanced/$', WorkflowAdvancedListView.as_view(), name="workflows_advanced"),
    url(r'^advanced/(?P<slug>[\w-]+)$', WorkflowAdvancedRedirectView.as_view(), name="workflows_advanced_step"),

    url(r'^alacarte$', workflows_alacarte_build, name="workflows_alacarte"),
    url(r'^wkmake/(?P<id>[\w-]+)$', WorkflowsMarkerRedirectView.as_view(), name="workflow_maker_form"),

]
