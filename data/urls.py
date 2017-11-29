from django.conf.urls import url

from .views import ImportPastedContentView, display_file, download_file, tree_visualization, export_file, \
    export_to_itol, display_msa, display_raw, get_example

urlpatterns = [
    url(r'^upload/$', ImportPastedContentView.as_view(), name='upload'),

    url(r'^display/(?P<file_id>[\w-]+)$', display_file, name="display_file"),
    url(r'^displayraw/(?P<file_id>[\w-]+)$', display_raw, name="display_raw"),
    url(r'^displaytree/(?P<file_id>[\w-]+)$', tree_visualization, name="display_tree"),
    url(r'^displaymsa/(?P<file_id>[\w-]+)$', display_msa, name="display_msa"),
    url(r'^example/$', get_example, name="get_example"),
    # url(r'^example/(?P<ext_file>[\w-]+)', get_example , name="get_example"),

    url(r'^download/(?P<file_id>[\w-]+)$', download_file, name="download_file"),
    url(r'^export_to_itol/(?P<file_id>[\w-]+)$', export_to_itol, name="export_to_itol"),
    url(r'^export/(?P<file_id>[\w-]+)$', export_file, name="export_file"),
]
