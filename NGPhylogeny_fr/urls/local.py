import debug_toolbar
from django.conf.urls import include
from django.conf.urls import url
from django.views.generic import TemplateView

from . import urlpatterns as base_url

urlpatterns = [
                  url(r'^__debug__/', include(debug_toolbar.urls)),
                  url(r'^error/404$', TemplateView.as_view(template_name="404.html"), name="error404"),
                  url(r'^error/500$', TemplateView.as_view(template_name="500.html"), name="error500"),
              ] + base_url
