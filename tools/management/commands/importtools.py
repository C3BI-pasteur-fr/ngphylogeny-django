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
        parser.add_argument('--force', action='store_true', help="force re-import tools")

    def handle(self, *args, **options):

        galaxy_url = options.get('galaxyurl')
        query = options.get('query')
        tool_id = options.get('id')
        tools_found = []

        if galaxy_url:
            galaxy_server, created = Server.objects.get_or_create(url=galaxy_url)

        else:
            try:
                galaxy_server = Server.objects.get(current=True)
            except Server.DoesNotExist:
                raise CommandError('Server Galaxy does not exist, please use --galaxyurl')

        tools_url = '%s/%s/%s/' % (galaxy_server.url, 'api', 'tools')

        if tool_id:
            # search specific tool
            tools_url_id = '%s/%s/' % (tools_url, tool_id)
            connection = requests.get(tools_url_id).json()
            if 'traceback' in connection:
                raise CommandError('Tool {} was not found in {}'.format(tool_id, galaxy_server.url))

            tools_found.append(connection.get('id'))

        else:
            # fetch list of tools
            connection = requests.get(tools_url, params={'in_panel': "false"})
            if connection.status_code == 200:
                tools_list = connection.json() or []

            for tool in tools_list:
                _m = [
                    tool.get('id').lower(),
                    tool.get('name').lower(),
                    str(tool.get('panel_section_name')).lower()
                ]
                if query in _m:
                    tools_found.append(tool.get('id'))

        if tools_found:

            self.stdout.write("%s" % ('\n'.join(tools_found)))
            response = raw_input('Do you want (re)import this tool(s)? [y/N] from {}:'.format(galaxy_server))

            if response.lower() == 'y':

                import_tools_report = Tool.import_tools(galaxy_server, tools=tools_found, force=options.get('force'))

                if import_tools_report['new']:
                    self.stdout.write(
                        self.style.SUCCESS("%s successfully imported new tools." % (len(import_tools_report['new']))))

                    for ti in import_tools_report['new']:
                        self.stdout.write(
                            self.style.SUCCESS("%s %s successfully imported" % (ti.name, ti.version))
                        )

                else:
                    self.stdout.write(self.style.WARNING("No new tool has been imported"))

                for tool_id in import_tools_report['error']:
                    self.stdout.write(
                        self.style.ERROR("import of %s has failed " % (tool_id))
                    )

        else:
            self.stdout.write("No result was found for query  %s" % (options.get('query')))
