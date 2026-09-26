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

from django.urls import re_path

from .views.generic import RerunWorkflow
from .views.wkoneclick import WorkflowOneClickListView
from .views.wkoneclick import WorkflowOneClickFormView
from .views.wkadvanced import WorkflowAdvancedListView
from .views.wkadvanced import WorkflowAdvancedFormView
from .views.wkmaker import WorkflowMakerView
from .views.wkmaker import workflows_alacarte_build

urlpatterns = [
    re_path(r'^oneclick/$', WorkflowOneClickListView.as_view(),
        name="workflow_oneclick_list"),
    re_path(r'^oneclick/(?P<slug>[\w-]+)$', WorkflowOneClickFormView.as_view(),
        name="workflow_oneclick_form"),
    re_path(r'^advanced/$', WorkflowAdvancedListView.as_view(),
        name="workflows_advanced"),
    re_path(r'^advanced/(?P<slug>[\w-]+)$',
        WorkflowAdvancedFormView.as_view(),
        name="workflows_advanced_fullsteps"),
    re_path(r'^alacarte$', workflows_alacarte_build,
        name="workflows_alacarte"),
    re_path(r'^rerun/(?P<id>[\w-]+)$', RerunWorkflow.as_view(),
        name="workflow_rerun"),
    re_path(r'^wkmake/(?P<id>[\w-]+)$', WorkflowMakerView.as_view(),
        name="workflow_maker_form"),
]
