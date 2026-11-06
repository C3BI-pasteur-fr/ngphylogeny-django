import hashlib
import json
import re
from functools import wraps
from urllib.parse import urlsplit

import requests
from bioblend.galaxy import GalaxyInstance as BioblendGalaxyInstance
from bioblend.galaxyclient import ConnectionError
from django.core.cache import cache
from requests import Request
from requests_toolbelt import MultipartEncoder


def is_cacheable_url(url):
    """
    Default-deny: only tool and workflow definitions, which nothing here
    edits in place. Anything carrying live state (dataset/job/history
    status) must stay fresh, and /build depends on the history it's asked
    about.
    """
    path = urlsplit(url).path
    return (bool(re.search(r'/api/workflows/[^/]+$', path))
            or ('/api/tools/' in path and not path.endswith('/build')))


def is_cacheable_response(response):
    """
    Galaxy's maintenance page arrives as a 200 with an HTML body - caching
    that would keep showing an outage for a full timeout after it ended.
    """
    if response.status_code != 200:
        return False
    try:
        response.json()
    except ValueError:
        return False
    return True


def get_request_cache(func):
    @wraps(func)
    def wrapper(self, url, **kwargs):

        if hasattr(self, "nocache") or not is_cacheable_url(url):
            return func(self, url, **kwargs)

        # Built from the request itself, never from response.url: bioblend
        # sends the API key as a header, so it isn't in the response's
        # URL, and this cache could never hit when it was stored that
        # way. Hashed so the key isn't readable in the cache table.
        params = dict(kwargs.get('params') or {})
        params['key'] = self.key
        request_url = Request('GET', url, params=params).prepare().url
        cache_key = 'galaxy_get:%s' % hashlib.sha256(
            request_url.encode()).hexdigest()

        cached = cache.get(cache_key)
        if cached is not None:
            return cached
        response = func(self, url, **kwargs)
        if is_cacheable_response(response):
            cache.set(cache_key, response)
        return response
    return wrapper


class GalaxyInstance(BioblendGalaxyInstance):
    """
    Override GalaxyInstance make_get_request method to use django cache
    """

    @get_request_cache
    def make_get_request(self, url, **kwargs):
        return super(GalaxyInstance, self).make_get_request(url, **kwargs)


class GalaxyInstanceAnonymous(GalaxyInstance):
    """
    Overload bioblend GalaxyClient methods to accept cookies on http requests
    """

    def __init__(self, url, galaxysession):
        super(GalaxyInstanceAnonymous, self).__init__(url, key=None, email=None, password=None)
        self.galaxysession = galaxysession
        # bioblend >=1.0 dropped GalaxyClient.default_params along with
        # query-string ?key= auth (replaced by the x-api-key header).
        # This class never had an API key to begin with - it
        # authenticates via the galaxysession cookie instead - so this
        # preserves the previous (effectively no-op) fallback exactly.
        self.default_params = {'key': self.key}

    def make_post_request(self, url, payload, params=None, files_attached=False, ):
        if params is not None and params.get('key', False) is False:
            params['key'] = self.key
        else:
            params = self.default_params

        # Compute data, headers, params arguments for request.post,
        # leveraging the requests-toolbelt library if any files have
        # been attached.
        if files_attached:
            payload.update(params)
            payload = MultipartEncoder(fields=payload)
            headers = self.json_headers.copy()
            headers['Content-Type'] = payload.content_type
            post_params = {}
        else:
            payload = json.dumps(payload)
            headers = self.json_headers
            post_params = params

        r = requests.post(url, data=payload, headers=headers,
                          verify=self.verify, params=post_params, cookies=dict(galaxysession=self.galaxysession))

        if r.status_code == 200:
            return r.json()
        # @see self.body for HTTP response body
        raise ConnectionError("Unexpected response from galaxy: %s" % r.status_code, body=r.text)

    def make_get_request(self, url, **kwargs):
        params = kwargs.get('params')
        if params is not None and params.get('key', False) is False:
            params['key'] = self.key
        else:
            params = self.default_params

        kwargs['cookies'] = dict(galaxysession=self.galaxysession)
        kwargs['params'] = params
        kwargs.setdefault('verify', self.verify)
        r = requests.get(url, **kwargs)
        return r
