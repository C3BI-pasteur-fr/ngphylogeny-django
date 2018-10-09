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
from .views import BlastView, BlastRunView, DeleteBlastRunView
from .views import DeleteBlastSubjectView, BlastRunFasta
from .views import DeleteBlastSequences


urlpatterns = [
    url(r'^$', BlastView.as_view(),
        name="blast_form"),
    url(r'^(?P<pk>[\w-]+)$', BlastRunView.as_view(),
        name="blast_view"),
    url(r'^(?P<pk>[\w-]+)/deletemulti$', DeleteBlastSequences.as_view(),
        name="blast_delete_seqs"),
    url(r'^subject/(?P<pk>[\w-]+)/delete$', DeleteBlastSubjectView.as_view(),
        name="blast_subject_delete"),
    url(r'^(?P<pk>[\w-]+)/fasta$', BlastRunFasta.as_view(),
        name="blast_fasta"),
    url(r'^(?P<pk>[\w-]+)/delete$', DeleteBlastRunView.as_view(),
        name="blast_delete"),

]
