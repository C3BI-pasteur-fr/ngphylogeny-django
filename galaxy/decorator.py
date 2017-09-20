from __future__ import unicode_literals

import logging

from django.http import Http404, HttpResponseGone
from django.shortcuts import redirect

from galaxy.models import Server, GalaxyUser

logger = logging.getLogger(__name__)


def connection_galaxy(view_function):
    """Initiating Galaxy connection"""

    def wrapper(request, *args, **kwargs):

        try:
            if hasattr(request, 'galaxy_server'):
                galaxy_server = request.galaxy_server

            elif request.session.get('galaxy_server'):
                galaxy_server = Server.objects.get(id=request.session.get('galaxy_server'))

            else:
                #by default use the current Galaxy server
                galaxy_server = Server.objects.get(current=True)

        except Server.DoesNotExist :
            msg = "NGPhylogeny server is not properly configured, " \
                  "please ensure that the Galaxy server is correctly set up"
            logger.exception(msg)
            raise Http404(msg)

        except Exception as e:
            logger.exception(e)
            raise HttpResponseGone()

        request.galaxy_server = galaxy_server
        request.session['galaxy_server'] = galaxy_server.id

        if request.user.is_authenticated():
            """Try to use related Galaxy user information"""

            try:
                """get or create Galaxy user onfly"""
                gu, created = GalaxyUser.objects.get_or_create(user=request.user, galaxy_server=galaxy_server)

                """If the key api is not defined, prompts the user to define it"""
                if gu.api_key:
                    request.galaxy = gu.get_galaxy_instance()
                else:
                    return redirect('galaxy_account')

            except GalaxyUser.DoesNotExist :
                msg = "NGPhylogeny server is not properly configured, " \
                      "please ensure that the Galaxy server is correctly set up"
                logger.exception("Galaxy User is not set")

                raise Http404(msg)

            except Exception as e:
                logger.exception("Galaxy account Error")
                return HttpResponseGone()

        elif request.user.is_anonymous():
            """If user is not an authenticated, use the anonymous Galaxy user set"""
            try:
                gu = GalaxyUser.objects.get(anonymous=True, galaxy_server=galaxy_server)
                request.galaxy = gu.get_galaxy_instance()

            except GalaxyUser.DoesNotExist :
                msg = "NGPhylogeny server is not properly configured, " \
                      "please ensure that the Galaxy server is correctly set"

                logger.exception("Anonymous user not set")
                raise Http404(msg)

            except Exception as e:
                logger.exception("Galaxy anonymous account Error: "+e)
                return HttpResponseGone()

        return view_function(request, *args, **kwargs)

    return wrapper
