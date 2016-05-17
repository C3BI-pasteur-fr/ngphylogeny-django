from __future__ import unicode_literals

from django.conf import settings
from galaxylib import GalaxyInstanceAnonymous, GalaxyInstance

import requests


class GalaxySessionMiddleware(object):
    """
    Middleware for utilizing Galaxy Anonymous authentication.
    """
    def process_request(self, request):

        if request.user.is_anonymous():
            galaxysession = request.session.setdefault(settings.GALAXY_SESSION_ID, {})

            if not galaxysession:
                galaxysession = self.new_connexion(request)

            request.galaxy = GalaxyInstanceAnonymous(url=settings.GALAXY_SERVER_URL, galaxysession=galaxysession)
        else:
            request.galaxy = GalaxyInstance(url=settings.GALAXY_SERVER_URL,
                                            key=request.user.galaxyuser.api_key)

    def new_connexion(self, request):
        """
        Try to connect with Galaxy server and save session information
        """
        connection = requests.get(settings.GALAXY_SERVER_URL)
        if connection.status_code == 200:

            for cookie in connection.cookies:
                galaxysession = str(cookie.value)
        else:
            raise Exception("Galaxy server response:" + connection.status.code)

        # update the galaxy session
        request.session[settings.GALAXY_SESSION_ID] = galaxysession
        # mark the session as "modified" to make sure it is saved
        request.session.modified = True
        return galaxysession

