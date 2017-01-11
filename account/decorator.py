from __future__ import unicode_literals

from galaxylib import  GalaxyInstance
from django.shortcuts import redirect, get_object_or_404
from account.models import GalaxyConf, GalaxyUser

# from views import galaxy_connection_error_view
# import requests
# from galaxylib import GalaxyInstanceAnonymous
# from django.conf import settings
#
# Allow anonymous access deprecated
# def get_galaxy_session_id():
#     """
#     Open an anonymous connection with Galaxy server and save Galaxy session information
#     """
#
#     connection = requests.get(settings.GALAXY_SERVER_URL)
#     # connect to galaxy server and retrieve galaxysession ID stored in a cookie
#     if connection.status_code == 200:
#         return connection.cookies.get('galaxysession')
#
#     else:
#         # TODO Galaxy no response page
#         raise Exception("Galaxy server response:" + connection.status.code)


def connection_galaxy(view_function):
    """Initiating Galaxy connection"""

    def wrapper(request, *args, **kwargs):

        try:
            galaxy_conf = get_object_or_404(GalaxyConf, active=True)
        except:

            raise Exception("NGPhylogeny server is not properly configured,"
                   "please ensure that the Galaxy server is correctly set up")

        if request.user.is_authenticated():
            """Try to use related Galaxy user information"""

            try:
                gu = GalaxyUser.objects.get(user=request.user, galaxy_server__galaxyconf__active=True )

                if gu.api_key:
                    request.galaxy = GalaxyInstance(url=galaxy_conf.galaxy_server.url, key=gu.api_key)
                else:
                    return redirect('account')

            except Exception, e:
                print e
                gu = GalaxyUser(user=request.user,
                                galaxy_server=galaxy_conf.galaxy_server)
                gu.save()
                return redirect('account')

        elif request.user.is_anonymous():

            """If user is not an authenticated user default galaxy user"""
            try:

                request.galaxy = GalaxyInstance(url=galaxy_conf.galaxy_server.url, key=galaxy_conf.global_api_key)

            except:
                return redirect('account')

                # Allow anonymous access deprecated
                # """If user is not an authenticated galaxyuser use cookies"""
                # try:
                #    id_galaxysession = request.session.setdefault(settings.GALAXY_SESSION_ID, get_galaxy_session_id())
                #    request.galaxy = GalaxyInstanceAnonymous(url=settings.GALAXY_SERVER_URL, galaxysession=id_galaxysession)
                #
                # except requests.ConnectionError:
                #
                #    return galaxy_connection_error_view(request)

        return view_function(request, *args, **kwargs)

    return wrapper
