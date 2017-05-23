from django.conf import settings
from django.core.management.base import BaseCommand

from galaxy.models import Server
from tools.models import Tool


class Command(BaseCommand):
    help = 'Import Galaxy tools to NGPhylogeny'

    def handle(self, *args, **options):

        self.handle_noargs()

    def handle_noargs(self):

        galaxy_server = ""
        try:

            galaxy_server = Server.objects.get(current=True)

        except:
            galaxy_server, created = Server.objects.get_or_create(url=settings.GALAXY_SERVER_URL)

        Tool.import_tools(galaxy_server)
