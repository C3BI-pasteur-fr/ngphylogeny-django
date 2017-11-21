import requests
from django.core.management.base import BaseCommand, CommandError

from galaxy.models import Server
from tools.models import Tool


class Command(BaseCommand):
    help = 'Import Galaxy tools to NGPhylogeny'
    requires_system_checks = True

    def add_arguments(self, parser):
        # Named (optional) arguments
        parser.add_argument('--galaxyurl')
        parser.add_argument('-q', '--query', help="Find tools according to query by default 'phylogeny' ")
        parser.add_argument('--id', help="Find tool according to galaxy tool id")

    def handle(self, *args, **options):

        galaxy_url = options.get('galaxyurl')
        query = options.get('query')
        tool_id = options.get('id')
        tools_list = []

        if galaxy_url:
            galaxy_server, created = Server.objects.get_or_create(url=galaxy_url)

        else:
            try:
                galaxy_server = Server.objects.get(current=True)
            except Server.DoesNotExist:
                raise CommandError('Server Galaxy does not exist, please use --galaxyurl')

        tools_url = '%s/%s/%s/' % (galaxy_server.url, 'api', 'tools')

        if query:
            connection = requests.get(tools_url, params={'q': query})
        else:
            connection = requests.get(tools_url, params={'q': "phylogeny"})

        if tool_id:
            tools_url_id = '%s/%s' % (tools_url, tool_id)
            connection = requests.get(tools_url_id)

        if connection.status_code == 200:
            tools_ids = connection.json()

            if tools_ids:

                tools_list = []

                if type(tools_ids) is list:
                    tools_list = tools_ids
                    self.stdout.write("%s" % ('\n'.join(tools_ids)))

                if type(tools_ids) is dict:
                    tools_list.append(tools_ids.get("id"))
                    self.stdout.write("%s" % (tools_ids.get("id")))

                response = raw_input('Do you want (re)import this tool(s)? [y/N]: ')

                if response.lower() == 'y':
                    tools_imported = Tool.import_tools(galaxy_server, tools=tools_list)
                    self.stdout.write(self.style.SUCCESS("%s successfully imported new tools." % (len(tools_imported))))

            else:
                self.stdout.write("No result was found for query  %s" % (options.get('query')))
