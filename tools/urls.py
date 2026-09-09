from django.urls import re_path

from .views import ToolListView, tool_exec_view, ToolDetailView, get_tool_name

urlpatterns = [
    re_path(r'^$', ToolListView.as_view(), name='tools'),
    re_path(r'^tool/(?P<pk>[\w-]+)/form$', tool_exec_view, name="tool_form"),
    re_path(r'^tool/(?P<pk>[\w-]+)$',
        ToolDetailView.as_view(), name="tool_detail"),
    re_path(r'^tool/galaxy_id/$', get_tool_name, name="get_tool_name")
]
