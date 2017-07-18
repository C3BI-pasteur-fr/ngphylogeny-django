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

from views import HistoryDetailView, PreviousHistoryListView, get_dataset_toolprovenance, GalaxyErrorView

urlpatterns = [
    url(r'^history$', HistoryDetailView.as_view(), name="history_current_detail"),
    url(r'^histories$', PreviousHistoryListView.as_view(), name="previous_analyses"),

    url(r'^history/(?P<history_id>[\w-]+)$', HistoryDetailView.as_view(), name="history_detail"),
    url(r'^history/provenance/(?P<history_id>[\w-]+)$', get_dataset_toolprovenance, name="get_dataset_tool"),
    url(r'^history/galaxyerror/(?P<id>[\w-]+)$', GalaxyErrorView.as_view(url='/dataset/errors?id=%(id)s'),
        name="galaxy_error_url"),

]
