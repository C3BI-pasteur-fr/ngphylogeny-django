from __future__ import unicode_literals

import logging

from django.http import Http404, HttpResponseGone

from galaxy.models import Server, GalaxyUser
from galaxy.galaxylib import GalaxyInstance

logger = logging.getLogger(__name__)


def connection_galaxy(view_function):
    """Initiating Galaxy connection"""

    def wrapper(request, *args, **kwargs):

        try:
            if request.session.get('server_change') is True:
                galaxy_server = Server.objects.get(current=True)
                request.session['server_change'] = False

            elif hasattr(request, 'galaxy_server'):
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
            return HttpResponseGone()

        request.galaxy_server = galaxy_server
        request.session['galaxy_server'] = galaxy_server.id

        # Every visitor - logged into an NGPhylogeny account or not -
        # authenticates to Galaxy through the same single, shared API
        # key (the "anonymous" GalaxyUser row): NGPhylogeny holds one
        # Galaxy identity, independent of any individual NGPhylogeny
        # account. This used to differ for authenticated users (each
        # expected their own personal GalaxyUser/api_key, redirecting
        # to set one up otherwise) - that redirect target
        # ('galaxy_account') was itself a dead URL name that's never
        # existed in this project's URLconf (the real name is
        # 'account'), so an authenticated user with no personal key
        # actually got an uncaught NoReverseMatch 500, not a clean
        # redirect - never caught because nothing exercised that path.
        # WorkspaceHistory.user (workspace/views.py's create_history())
        # still correctly tracks the real logged-in NGPhylogeny account
        # separately from this - this only concerns which Galaxy
        # credential is used to actually talk to Galaxy.
        try:
            gu = GalaxyUser.objects.get(anonymous=True, galaxy_server=galaxy_server)
            request.galaxy = gu.get_galaxy_instance

        except GalaxyUser.DoesNotExist:
            msg = "NGPhylogeny server is not properly configured, " \
                  "please ensure that the Galaxy server is correctly set up"
            logger.exception("Anonymous user not set")
            raise Http404(msg)

        except Exception as e:
            logger.exception("Galaxy account error: %s" % e)
            return HttpResponseGone()

        return view_function(request, *args, **kwargs)

    return wrapper

def galaxy_connection():
    """Initiating Galaxy connection"""
    galaxy_instance = None
    try:
        galaxy_server = Server.objects.get(current=True)
        gu = GalaxyUser.objects.get(anonymous=True, galaxy_server=galaxy_server)
        galaxy_instance = gu.get_galaxy_instance
    except Server.DoesNotExist as s:
        msg = "NGPhylogeny server is not properly configured, " \
              "please ensure that the Galaxy server is correctly set up"
        logger.exception(msg)
        raise s
    except Exception as e:
        logger.exception(e)
        raise e
    return galaxy_instance

def galaxy_connection_simple(server, key):
    """Initiating Galaxy connection"""
    return GalaxyInstance(url=server, key=key)
