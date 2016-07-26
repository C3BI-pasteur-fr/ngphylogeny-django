from django.conf.urls import url
from views import *

urlpatterns = [
    url(r'^$', ToolListView.as_view(), name='tools'),

    url(r'^id/(?P<pk>[\w-]+)$', ToolDetailView.as_view(), name="tool_detail"),
    url(r'^id/(?P<pk>[\w-]+)/json$', ToolJSONView.as_view(), name="tool_detail_json"),
]
