from __future__ import unicode_literals

from django.shortcuts import redirect, get_object_or_404
from account.models import GalaxyConf, GalaxyUser


def connection_galaxy(view_function):
    """Initiating Galaxy connection"""

    def wrapper(request, *args, **kwargs):

        try:

            if hasattr(request, 'galaxyconf'):
                galaxy_conf = request.galaxyconf

            else:
                galaxy_conf = get_object_or_404(GalaxyConf, active=True)
                request.galaxyconf = galaxy_conf

        except Exception, e:
            print e
            raise Exception("NGPhylogeny server is not properly configured,"
                            "please ensure that the Galaxy server is correctly set up")

        if request.user.is_authenticated():
            """Try to use related Galaxy user information"""

            try:
                """get or create galaxy user onfly"""
                gu, created = GalaxyUser.objects.get_or_create(user=request.user, galaxy_server=galaxy_conf.galaxy_server)

                """If the key api is not defined, prompts the user to define it"""
                if gu.api_key:
                    request.galaxy = gu.get_galaxy_instance()
                else:
                    return redirect('account')

            except Exception, e:
                print e
                return redirect('account')

        elif request.user.is_anonymous():
            """If user is not an authenticated, use anonymous galaxyuser set"""
            try:
                request.galaxy = galaxy_conf.anonymous_user.get_galaxy_instance()

            except:
                return redirect('account')

        return view_function(request, *args, **kwargs)

    return wrapper
