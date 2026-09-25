from django.urls import re_path
from django.contrib.auth.views import LogoutView

from .views import AccountCreateView, AccountDeleteView, AccountDetailView, AccountLoginView

urlpatterns = [
    re_path(r'^login$', AccountLoginView.as_view(), name='login'),
    re_path(r'logout$', LogoutView.as_view(next_page="/"), name='logout'),
    re_path(r'^$', AccountDetailView.as_view(), name='account'),
    re_path(r'^create$', AccountCreateView.as_view(), name='create_account'),
    re_path(r'^delete$', AccountDeleteView.as_view(), name='delete_account'),
]