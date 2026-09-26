from django.urls import re_path
from django.contrib.auth.views import LogoutView

from .views import (
    AccountCreateView, AccountDeleteView, AccountDetailView, AccountLoginView,
    AccountPasswordResetView, AccountPasswordResetDoneView,
    AccountPasswordResetConfirmView, AccountPasswordResetCompleteView)

urlpatterns = [
    re_path(r'^login$', AccountLoginView.as_view(), name='login'),
    re_path(r'logout$', LogoutView.as_view(next_page="/"), name='logout'),
    re_path(r'^$', AccountDetailView.as_view(), name='account'),
    re_path(r'^create$', AccountCreateView.as_view(), name='create_account'),
    re_path(r'^delete$', AccountDeleteView.as_view(), name='delete_account'),
    re_path(r'^password-reset$',
        AccountPasswordResetView.as_view(), name='password_reset'),
    re_path(r'^password-reset/done$',
        AccountPasswordResetDoneView.as_view(), name='password_reset_done'),
    # [^/]+ each, matching django.contrib.auth.urls's own reference
    # urlconf exactly (path("reset/<uidb64>/<token>/", ...) - a plain
    # str converter, not a stricter hand-picked regex).
    re_path(r'^reset/(?P<uidb64>[^/]+)/(?P<token>[^/]+)/$',
        AccountPasswordResetConfirmView.as_view(), name='password_reset_confirm'),
    re_path(r'^reset/done$',
        AccountPasswordResetCompleteView.as_view(), name='password_reset_complete'),
]