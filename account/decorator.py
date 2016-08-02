from __future__ import unicode_literals

from django.conf import settings
from galaxylib import GalaxyInstanceAnonymous, GalaxyInstance
from django.core.exceptions import ObjectDoesNotExist

from views import galaxy_connection_error_view

import requests


def get_galaxy_session_id():
    """
    Try to connect with Galaxy server and return galaxy session information
    """
    # TODO Galaxy no response page
    id_galaxysession = ""
    connection = requests.get(settings.GALAXY_SERVER_URL)
    # connect to galaxy server and retrieve galaxysession id stored in a cookie
    if connection.status_code == 200:
        for cookie in connection.cookies:
            id_galaxysession = str(cookie.value)
    else:
        raise Exception("Galaxy server response:" + connection.status.code)

    return id_galaxysession


def connection_galaxy(view_function):
    """Initiating Galaxy connection"""

    def wrapper(request, *args, **kwargs):

        try:
            if request.user.is_anonymous():
                """if user is anonymous galaxyuser does not exit"""
                raise ObjectDoesNotExist

            request.galaxy = GalaxyInstance(url=settings.GALAXY_SERVER_URL,
                                            key=request.user.galaxyuser.api_key)
        except ObjectDoesNotExist:
            """if user is not galaxyuser"""
            try:
                id_galaxysession = request.session.setdefault(settings.GALAXY_SESSION_ID,
                                                              get_galaxy_session_id())

                request.galaxy = GalaxyInstanceAnonymous(url=settings.GALAXY_SERVER_URL,
                                                         galaxysession=id_galaxysession)
            except requests.ConnectionError:

                return galaxy_connection_error_view(request)

        return view_function(request, *args, **kwargs)

    return wrapper
