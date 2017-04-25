from django.conf.urls import url

from views import *

urlpatterns = [
    url(r'^upload/$', UploadView.as_view(), name='upload'),

    url(r'^display/(?P<file_id>[\w-]+)$', display_file, name="display_file"),
    url(r'^download/(?P<file_id>[\w-]+)$', download_file, name="download_file"),
    url(r'^export/(?P<file_id>[\w-]+)$', export_file, name="export_file"),
]
