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

from .views import HistoryDetailView, PreviousHistoryListView, \
    get_dataset_toolprovenance, GalaxyErrorView, \
    get_dataset_citations, \
    WorkspaceDeleteView, WorkspaceRenameView, \
    WorkspaceChangeEmailView, DeleteAllHistories

urlpatterns = [
    url(r'^history$',
        HistoryDetailView.as_view(), name="history_current_detail"),
    url(r'^histories$',
        PreviousHistoryListView.as_view(), name="previous_analyses"),
    url(r'^history/(?P<history_id>[\w-]+)$',
        HistoryDetailView.as_view(), name="history_detail"),
    url(r'^history/(?P<history_id>[\w-]+)/rename$',
        WorkspaceRenameView.as_view(), name="history_rename"),
    url(r'^history/(?P<history_id>[\w-]+)/email$',
        WorkspaceChangeEmailView.as_view(), name="change_email"),
    url(r'^history/(?P<history_id>[\w-]+)/delete$',
        WorkspaceDeleteView.as_view(), name="history_delete"),
    url(r'^histories/deleteall$',
        DeleteAllHistories.as_view(), name="history_delete_all"),
    url(r'^history/provenance/(?P<history_id>[\w-]+)$',
        get_dataset_toolprovenance, name="get_dataset_tool"),
    url(r'^history/citations/(?P<history_id>[\w-]+)$',
        get_dataset_citations, name="get_dataset_citations"),
    url(r'^history/galaxyerror/(?P<id>[\w-]+)$',
        GalaxyErrorView.as_view(),
        name="galaxy_error_url"),
]
