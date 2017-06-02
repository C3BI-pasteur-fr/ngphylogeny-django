from django.conf.urls import url

from views import *

urlpatterns = [
    url(r'^upload/$', UploadView.as_view(), name='upload'),

    url(r'^display/(?P<file_id>[\w-]+)$', display_file, name="display_file"),
    url(r'^displaytree/(?P<file_id>[\w-]+)$', tree_visualization, name="display_tree"),
    url(r'^download/(?P<file_id>[\w-]+)$', download_file, name="download_file"),
    url(r'^export_to_itol/(?P<file_id>[\w-]+)$', export_to_itol, name="export_to_itol"),
    url(r'^export/(?P<file_id>[\w-]+)$', export_file, name="export_file"),
]
