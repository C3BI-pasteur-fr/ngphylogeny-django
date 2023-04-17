import json
from functools import wraps

import requests
from bioblend.galaxy import GalaxyInstance as BioblendGalaxyInstance
from bioblend.galaxyclient import ConnectionError
from django.core.cache import cache
from requests import Request
from requests_toolbelt import MultipartEncoder


def get_request_cache(func):
    @wraps(func)
    def wrapper(self, url, **kwargs):

        if hasattr(self, "nocache"):
            return func(self, url, **kwargs)

        params = dict(kwargs.get('params') or {})
        params['key'] = self.key
        pre_req = Request('GET', url, params=params).prepare()

        c = cache.get(pre_req.url)
        if c:
            return c
        else:
            r = func(self, url, **kwargs)
            cache.set(r.url, r)
        return r
    return wrapper


class GalaxyInstance(BioblendGalaxyInstance):
    """
    Override GalaxyInstance make_get_request method to use django cache
    """
    def __init__(self, url, galaxysession):
        super(BioblendGalaxyInstance, self).__init__(url, key=None, email=None, password=None)
        self.json_header['x-api-key'] = self.key
        
    
    @get_request_cache
    def make_get_request(self, url, **kwargs):
        """
        Make a GET request using the provided ``url``.
        Keyword arguments are the same as in requests.request.
        If ``verify`` is not provided, ``self.verify`` will be used.
        If the ``params`` are not provided, use ``default_params`` class field.
        If params are provided and the provided dict does not have ``key`` key,
        the default ``self.key`` value will be included in what's passed to
        the server via the request.
        :rtype: requests.Response
        :return: the response object.
        """
        params = kwargs.get('params')
        if params is not None and params.get('key', False) is False:
            params['key'] = self.key
        else:
            params = self.default_params
        headers=self.json_headers
        kwargs['params'] = params
        kwargs.setdefault('verify', self.verify)
        kwargs.setdefault('timeout', self.timeout)
        r = requests.get(url, headers=headers, **kwargs)
        return r


class GalaxyInstanceAnonymous(GalaxyInstance):
    """
    Overload bioblend GalaxyClient methods to accept cookies on http requests
    """

    def __init__(self, url, galaxysession):
        super(GalaxyInstanceAnonymous, self).__init__(url, key=None, email=None, password=None)
        self.galaxysession = galaxysession

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
        
        headers=self.json_headers
        kwargs['cookies'] = dict(galaxysession=self.galaxysession)
        kwargs['params'] = params
        kwargs.setdefault('verify', self.verify)
        r = requests.get(url, headers=headers, **kwargs)
        return r
