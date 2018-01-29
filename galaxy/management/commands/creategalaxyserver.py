import requests
from django.core.management.base import BaseCommand

from galaxy.models import Server


class Command(BaseCommand):
    help = 'Create or activate Galaxy Server to NGPhylogeny'
    requires_system_checks = True

    def add_arguments(self, parser):
        # Named (optional) arguments
        parser.add_argument('--url')
        parser.add_argument('--name', help="Server Name")
        parser.add_argument('--activate', action='store_true', help="force re-import tools")
        parser.add_argument('--interactive', action='store_true')

    def handle(self, *args, **options):

        galaxy_url = options.get('galaxyurl')
        galaxy_name = options.get('name')
        current = bool(options.get('activate'))
        interactive = bool(options.get('interactive'))

        if interactive:
            new_url = raw_input('Enter the URL of the Galaxy server you want to add: ')
            galaxy_url = requests.get(new_url).url

        if galaxy_url:
            galaxy_server, created = Server.objects.get_or_create(url=galaxy_url)

            if created:
                if not galaxy_name:
                    galaxy_name = "Galaxy Server " + str(len(Server.objects.all()))

                galaxy_server.name = galaxy_name
                self.stdout.write(
                    self.style.SUCCESS("Galaxy Server {} created as {}".format(galaxy_server.url, galaxy_name)))

            if current:
                galaxy_server.current = current
                self.stdout.write(self.style.SUCCESS("Galaxy Server {} activated".format(galaxy_server.name)))

            galaxy_server.save()
