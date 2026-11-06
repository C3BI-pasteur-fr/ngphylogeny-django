from django.conf import settings
from django.template.response import TemplateResponse

from galaxy.outage import is_galaxy_outage
from workflows.views.generic import galaxy_unavailable_response


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

    /status is deliberately exempt: manifest.yaml's web Deployment points
    both its readinessProbe and livenessProbe at it. Without this
    exemption, enabling maintenance mode made the probes see the same 503
    the middleware sends everywhere else, read that as "the container is
    broken" (not "intentionally in maintenance"), and endlessly restart
    it via failed liveness checks - the Service ends up with zero ready
    endpoints and the site goes fully down (HAProxy's "no available
    server") *because* of turning maintenance mode on, the opposite of
    the intended effect. Only caught by actually enabling
    MAINTENANCE=true against a real deployment, not by any test written
    before that (which exercised the middleware/template, not the
    Kubernetes probe interaction).
    """
    EXEMPT_PATHS = {'/status'}

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        if (settings.NGPHYLO_MAINTENANCE_MODE
                and request.path not in self.EXEMPT_PATHS):
            # Rendered explicitly, not left to Django's usual
            # auto-render-a-TemplateResponse step - that step only
            # applies to the innermost view call, not to whatever an
            # outer middleware like this one returns without ever
            # calling self.get_response(). Without this, .content stays
            # unrendered (SimpleTemplateResponse raises
            # ContentNotRenderedError on access rather than serving it
            # empty) and .templates stays empty - caught by
            # assertTemplateUsed in tests before ever reaching a real
            # deployment with MAINTENANCE=true.
            return TemplateResponse(
                request, 'maintenance.html', status=503).render()
        return self.get_response(request)


class GalaxyUnavailableMiddleware:
    """
    Turns a Galaxy outage raised by any view into the friendly 503 page
    instead of a 500 - one place rather than a try/except around each of
    the many views that call Galaxy. Only outages, see is_galaxy_outage:
    anything else still propagates to Django's normal 500 handling.

    Must stay last in MIDDLEWARE (settings/base.py): process_exception
    runs innermost-first, and the page needs the session/auth/messages
    middleware to have already run for base.html to render.
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        return self.get_response(request)

    def process_exception(self, request, exception):
        if is_galaxy_outage(exception):
            return galaxy_unavailable_response(request)
        return None
