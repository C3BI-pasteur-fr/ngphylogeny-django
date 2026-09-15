from django.conf import settings
from django.template.response import TemplateResponse


class MaintenanceModeMiddleware:
    """
    Serves templates/maintenance.html for every request, instead of
    routing normally, whenever settings.NGPHYLO_MAINTENANCE_MODE is True -
    toggled via the NGPHYLO_MAINTENANCE_MODE env var (see CLAUDE.md's
    "Kubernetes deployment" section for the MAINTENANCE GitLab CI/CD
    variable that drives it there). A 503 status, not 200 - the correct
    code for "intentionally unavailable", distinguishable by monitoring/
    uptime checks from either a real error or the site actually being up.

    Must stay after WhiteNoiseMiddleware in MIDDLEWARE (settings/base.py):
    whitenoise fully handles and returns static asset requests before
    they'd ever reach this middleware, which is exactly what lets
    maintenance.html's own static assets (CSS, the maintenance image)
    keep loading while everything else is blocked.
    """
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        if settings.NGPHYLO_MAINTENANCE_MODE:
            return TemplateResponse(
                request, 'maintenance.html', status=503)
        return self.get_response(request)
