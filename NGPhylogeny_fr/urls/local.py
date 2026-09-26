# import debug_toolbar
from django.urls import include, re_path
from django.views.generic import TemplateView

from . import urlpatterns as base_url

urlpatterns = [
                  # url(r'^__debug__/', include(debug_toolbar.urls)),
                  re_path(r'^error/404$', TemplateView.as_view(template_name="404.html"), name="error404"),
                  re_path(r'^error/500$', TemplateView.as_view(template_name="500.html"), name="error500"),
              ] + base_url
