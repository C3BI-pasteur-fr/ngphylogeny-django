import os

from .base import *

DEBUG = False

# Set when deployed behind a TLS-terminating reverse proxy (e.g. Caddy) -
# see SECURE_PROXY_SSL_HEADER in base.py, which makes Django trust that
# proxy's X-Forwarded-Proto header to know a request is actually HTTPS.
# Without CSRF_TRUSTED_ORIGINS matching, Django's CSRF check computes the
# expected origin as http://... (unaware of the proxy) while the browser
# sends Origin: https://..., failing every POST with "Referer checking
# failed". SESSION/CSRF_COOKIE_SECURE are opt-in here (not in base.py)
# because turning them on without HTTPS actually being served would break
# cookies entirely - fine behind Caddy, wrong for a plain-HTTP deployment.
NGPHYLO_HTTPS_HOST = os.environ.get('NGPHYLO_HTTPS_HOST')
if NGPHYLO_HTTPS_HOST:
    CSRF_TRUSTED_ORIGINS = ['https://' + NGPHYLO_HTTPS_HOST]
    SESSION_COOKIE_SECURE = True
    CSRF_COOKIE_SECURE = True
