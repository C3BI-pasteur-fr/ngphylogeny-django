from django.conf import settings
from django.core.management.base import BaseCommand

from account.models import GalaxyServer, GalaxyConf
from tools.models import Tool


class Command(BaseCommand):
    help = 'Import Galaxy tools to NGPhylogeny'

    def handle(self, *args, **options):

        self.handle_noargs()

    def handle_noargs(self):

        galaxy_server = ""
        try:
            galaxy_active_conf = GalaxyConf.objects.get(active=True)
            galaxy_server = galaxy_active_conf.galaxy_server

        except:
            galaxy_server, created = GalaxyServer.objects.get_or_create(url=settings.GALAXY_SERVER_URL)

        Tool.import_tools(galaxy_server)