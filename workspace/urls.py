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

from .views import HistoryDetailView, HistoryContentRefreshView, PreviousHistoryListView, \
    get_dataset_toolprovenance, GalaxyErrorView, \
    get_dataset_citations, get_dataset_citations_bibtex, get_dataset_citations_txt, \
    WorkspaceDeleteView, WorkspaceRenameView, \
    WorkspaceChangeEmailView, DeleteAllHistories, daily_report_view, \
    running_jobs_view, WorkspacePermalinkView, export_rocrate

urlpatterns = [
    re_path(r'^report$',
        daily_report_view, name="daily_report"),
    re_path(r'^running$',
        running_jobs_view, name="running_jobs"),
    re_path(r'^history$',
        HistoryDetailView.as_view(), name="history_current_detail"),
    re_path(r'^histories$',
        PreviousHistoryListView.as_view(), name="previous_analyses"),
    re_path(r'^permalink/(?P<token>[\w:=-]+)$',
        WorkspacePermalinkView.as_view(), name="workspace_permalink"),
    re_path(r'^history/(?P<history_id>[\w-]+)$',
        HistoryDetailView.as_view(), name="history_detail"),
    re_path(r'^history/(?P<history_id>[\w-]+)/refresh$',
        HistoryContentRefreshView.as_view(), name="history_content_refresh"),
    re_path(r'^history/(?P<history_id>[\w-]+)/rename$',
        WorkspaceRenameView.as_view(), name="history_rename"),
    re_path(r'^history/(?P<history_id>[\w-]+)/email$',
        WorkspaceChangeEmailView.as_view(), name="change_email"),
    re_path(r'^history/(?P<history_id>[\w-]+)/delete$',
        WorkspaceDeleteView.as_view(), name="history_delete"),
    re_path(r'^histories/deleteall$',
        DeleteAllHistories.as_view(), name="history_delete_all"),
    re_path(r'^history/provenance/(?P<history_id>[\w-]+)$',
        get_dataset_toolprovenance, name="get_dataset_tool"),
    re_path(r'^history/citations/(?P<history_id>[\w-]+)$',
        get_dataset_citations, name="get_dataset_citations"),
    re_path(r'^history/citations/bibtex/(?P<history_id>[\w-]+)$',
        get_dataset_citations_bibtex, name="get_dataset_citations_bibtex"),
    re_path(r'^history/citations/text/(?P<history_id>[\w-]+)$',
        get_dataset_citations_txt, name="get_dataset_citations_txt"),
    re_path(r'^history/(?P<history_id>[\w-]+)/rocrate$',
        export_rocrate, name="history_rocrate"),
    re_path(r'^history/galaxyerror/(?P<id>[\w-]+)$',
        GalaxyErrorView.as_view(),
        name="galaxy_error_url"),
]
