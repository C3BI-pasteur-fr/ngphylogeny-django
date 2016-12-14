from __future__ import unicode_literals

from django.conf import settings
from galaxylib import GalaxyInstanceAnonymous, GalaxyInstance
from django.core.exceptions import ObjectDoesNotExist

from views import galaxy_connection_error_view, galaxy_user_apikey_settings

import requests


def get_galaxy_session_id():
    """
    Open an anonymous connection with Galaxy server and save Galaxy session information
    """
    connection = requests.get(settings.GALAXY_SERVER_URL)
    # connect to galaxy server and retrieve galaxysession ID stored in a cookie
    if connection.status_code == 200:
        return connection.cookies.get('galaxysession')

    else:
        # TODO Galaxy no response page
        raise Exception("Galaxy server response:" + connection.status.code)


def connection_galaxy(view_function):
    """Initiating Galaxy connection"""

    def wrapper(request, *args, **kwargs):

        if request.user.is_authenticated():
            """Try to use related Galaxy user information"""
            if request.user.galaxyuser.api_key:
                request.galaxy = GalaxyInstance(url=settings.GALAXY_SERVER_URL, key=request.user.galaxyuser.api_key)
            else:

                return galaxy_user_apikey_settings(request)


        elif request.user.is_anonymous():
            """If user is not an authenticated galaxyuser"""
            try:
                id_galaxysession = request.session.setdefault(settings.GALAXY_SESSION_ID, get_galaxy_session_id())
                request.galaxy = GalaxyInstanceAnonymous(url=settings.GALAXY_SERVER_URL, galaxysession=id_galaxysession)

            except requests.ConnectionError:

                return galaxy_connection_error_view(request)

        return view_function(request, *args, **kwargs)

    return wrapper
