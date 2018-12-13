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
from .views import available_blasts_progs, available_blasts_dbs
from .views import available_blasts_servers, blast_example


urlpatterns = [
    url(r'^dbs/(?P<server>\w+)/(?P<prog>[\w/\.]+)$',
        available_blasts_dbs, name="available_blasts_dbs"),
    url(r'^progs/(?P<server>\w+)$',
        available_blasts_progs, name="available_blasts_progs"),
    url(r'^servers$',
        available_blasts_servers, name="available_blasts_servers"),
    url(r'^example/(?P<server>\w+)/(?P<prog>[\w/\.]+)$',
        blast_example, name="blast_example"),
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

def available_blasts_progs(request, server):
    """
    Ajax: return possible blasts progs : {id:name}
    """
    context = BlastRun.blast_progs(server)
    return HttpResponse(json.dumps(context), content_type='application/json')
