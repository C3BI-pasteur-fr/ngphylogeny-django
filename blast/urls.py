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
from .views import BlastView, BlastRunView, DeleteBlastRunView
from .views import DeleteBlastSubjectView, BlastRunFasta
from .views import DeleteBlastSequences
from .views import available_blasts_progs, available_blasts_dbs
from .views import available_blasts_servers, blast_example


urlpatterns = [
    # prog is a Galaxy toolshed tool id (e.g.
    # "toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/
    # ncbi_blastn_wrapper/2.14.1+galaxy2" for Pasteur's current wrappers -
    # see settings.BLASTS) inserted as-is into the URL by
    # templates/blast/blast.html's JS (url.replace('wildcard2', progid),
    # not encodeURIComponent-escaped). [\w/\.]+ (\w = letters/digits/_)
    # didn't allow the literal "+" a Galaxy version suffix like
    # "+galaxy2" contains, or "-" (common in toolshed owner/repo names,
    # not currently used here but a very plausible future one) - either
    # 404s this route entirely, breaking both the database dropdown and
    # the example-sequence button for any prog whose id has one.
    re_path(r'^dbs/(?P<server>\w+)/(?P<prog>[\w/.+-]+)$',
        available_blasts_dbs, name="available_blasts_dbs"),
    re_path(r'^progs/(?P<server>\w+)$',
        available_blasts_progs, name="available_blasts_progs"),
    re_path(r'^servers$',
        available_blasts_servers, name="available_blasts_servers"),
    re_path(r'^example/(?P<server>\w+)/(?P<prog>[\w/.+-]+)$',
        blast_example, name="blast_example"),
    re_path(r'^$', BlastView.as_view(),
        name="blast_form"),
    re_path(r'^(?P<pk>[\w-]+)$', BlastRunView.as_view(),
        name="blast_view"),
    re_path(r'^(?P<pk>[\w-]+)/deletemulti$', DeleteBlastSequences.as_view(),
        name="blast_delete_seqs"),
    re_path(r'^subject/(?P<pk>[\w-]+)/delete$', DeleteBlastSubjectView.as_view(),
        name="blast_subject_delete"),
    re_path(r'^(?P<pk>[\w-]+)/fasta$', BlastRunFasta.as_view(),
        name="blast_fasta"),
    re_path(r'^(?P<pk>[\w-]+)/delete$', DeleteBlastRunView.as_view(),
        name="blast_delete"),
]
