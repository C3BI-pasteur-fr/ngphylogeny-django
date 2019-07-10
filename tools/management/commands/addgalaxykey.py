import requests
import re
from django.core.management.base import BaseCommand, CommandError
from django.contrib.auth.models import User;
from slugify import slugify

from galaxy.models import Server
from galaxy.models import GalaxyUser
from workflows.models import Workflow
from tools.models import Tool
from tools.models import ToolFlag


class Command(BaseCommand):
    help = 'Adds Galaxy Key to given user and given galaxy server'
    requires_system_checks = True

    def add_arguments(self, parser):
        parser.add_argument('--galaxyurl')
        parser.add_argument('--user')
        parser.add_argument('--galaxykey')

    def handle(self, *args, **options):
        galaxy_url = options.get('galaxyurl')
        galaxy_key = options.get('galaxykey')
        user = options.get('user')

        galaxy_server = None
        if galaxy_key is None:
            raise CommandError(
                    'Please give a galaxy key with --galaxykey')
        
        if galaxy_url:
            galaxy_server, created = Server.objects.get_or_create(
                url=galaxy_url)
        else:
            try:
                galaxy_server = Server.objects.get(current=True)
            except Server.DoesNotExist:
                raise CommandError(
                    'Server Galaxy does not exist, please use --galaxyurl')
        
        u = User.objects.filter(username=user)
        if u.count() == 0:
            raise CommandError('User {} does not exist, please use --user'.format(user))
        
        g = Server.objects.filter(url=galaxy_url)
        if g.count() == 0:
            raise CommandError('Galaxy server {} does not exist, please use --galaxyurl'.format(galaxy_url))

        try:
            gu=GalaxyUser.objects.get(user=u.first())
            gu.galaxy_server=g.first()
            gu.api_key=galaxy_key
            gu.anonymous=True
            gu.save()
            self.stdout.write("User %s already present" % user)
        except GalaxyUser.DoesNotExist:
            gu=GalaxyUser(user=u.first(),galaxy_server=g.first(),api_key=galaxy_key,anonymous=True)
            gu.save()
