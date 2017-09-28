import json

import requests
from bioblend.galaxy import GalaxyInstance as GalaxyI
from bioblend.galaxyclient import ConnectionError
from django.core.cache import cache
from django.utils.decorators import method_decorator
from requests import Request
from requests_toolbelt import MultipartEncoder


def get_request_cache(func):
    def wrapper(url, **kwargs):
        pre_req = Request('GET', url, params=kwargs.get('params', '')).prepare()
        print pre_req.url
        c = cache.get(pre_req.url)
        if c:
            return c
        else:
            r = func(url, **kwargs)
            cache.set(pre_req.url, r, 60 * 15)
        return r

    return wrapper


@method_decorator(get_request_cache, name="make_get_request")
class GalaxyInstance(GalaxyI):
    """
    Overlaod GalaxyInstance make_get_request method to use django cache
    """
    pass


class GalaxyInstanceAnonymous(GalaxyInstance):
    """
    Overriding of bioblend GalaxyClient methods to accept cookies on http requests
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

        kwargs['cookies'] = dict(galaxysession=self.galaxysession)
        kwargs['params'] = params
        kwargs.setdefault('verify', self.verify)
        r = requests.get(url, **kwargs)
        return r
