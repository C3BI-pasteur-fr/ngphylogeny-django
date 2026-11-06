"""
Telling a Galaxy outage apart from a real answer, and backing off while
one lasts - see CLAUDE.md's "Cutting Galaxy load" section.
"""
import requests
from bioblend.galaxy.client import ConnectionError
from django.core.cache import cache

# Raised by bioblend itself (ConnectionError, once its own retries are
# exhausted) or by the underlying requests call (e.g. a ReadTimeout - see
# GalaxyUser.GALAXY_REQUEST_TIMEOUT in galaxy/models.py) whenever the
# configured Galaxy server can't actually be reached (e.g. Galaxy's own
# maintenance mode) - see CLAUDE.md's "pod restart" section for why
# ConnectionError alone isn't enough.
GALAXY_UNREACHABLE_EXCEPTIONS = (
    ConnectionError, requests.exceptions.RequestException)

GALAXY_POLLING_PAUSE_SECONDS = 60
GALAXY_POLLING_PAUSE_KEY = 'galaxy_polling_paused'


def is_galaxy_outage(exception):
    """
    True if this exception means Galaxy can't be reached right now, as
    opposed to a real answer from Galaxy about this specific request.
    bioblend raises the same ConnectionError for both: a maintenance page
    arrives as a 200 with a non-JSON body, or a 302 to that page on a
    POST (both seen in production), while a 4xx (a stale id, a dataset
    we can't access) is a genuine reply that "try again in a few minutes"
    would only misdescribe, and that shouldn't be hidden from error logs.
    """
    if not isinstance(exception, GALAXY_UNREACHABLE_EXCEPTIONS):
        return False
    status_code = getattr(exception, 'status_code', None)
    return status_code is None or not 400 <= status_code < 500


def pause_galaxy_polling():
    """
    Background polling only: anything a visitor explicitly asks for
    (submitting a run, opening a file) still goes through.
    """
    cache.set(GALAXY_POLLING_PAUSE_KEY, True, GALAXY_POLLING_PAUSE_SECONDS)


def galaxy_polling_paused():
    return cache.get(GALAXY_POLLING_PAUSE_KEY) is True
