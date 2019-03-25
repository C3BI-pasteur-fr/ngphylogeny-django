import requests
import sys

from django.core.management.base import BaseCommand, CommandError

from galaxy.models import Server
from tools.models import Tool
from tools.models import ToolFlag
from tools.models import ToolInputData

class Command(BaseCommand):
    help = 'Import Galaxy tools to NGPhylogeny'
    requires_system_checks = True
    input_fields = []
    flags = []

    def add_arguments(self, parser):
        # Named (optional) arguments
        parser.add_argument('--galaxyurl')
        parser.add_argument(
            '-q', '--query',
            help="Find tools according to query by default 'phylogeny' ")
        parser.add_argument(
            '--id', help="Find tool according to galaxy tool id")
        parser.add_argument('--force', action='store_true',
                            help="force re-import tools")
        parser.add_argument(
            '--flags',
            help="Imports mapping between tools and flags using given file")
        parser.add_argument(
            '--inputfields',
            help="Imports input fields from an input file")

    def read_flag_file(self, flagfile):
        with open(flagfile) as f:
            for line in f:
                link = line.strip().split("\t")
                self.flags.append(link)

    def read_input_fields_file(self, fieldfile):
        with open(fieldfile) as f:
            for line in f:
                link = line.strip().split("\t")
                self.input_fields.append(link)

    def handle(self, *args, **options):
        galaxy_url = options.get('galaxyurl')
        flagfile = options.get('flags')
        fieldfile = options.get('inputfields')
        query = options.get('query')
        force = options.get('force')
        tool_id = options.get('id')

        self.read_flag_file(flagfile)
        self.read_input_fields_file(fieldfile)

        self.import_tools(galaxy_url, query, tool_id, force)
        self.associate_flags()
        self.set_input_fields()

    def import_tools(self, galaxy_url, query, tool_id, force):
        tools_found = []

        if galaxy_url:
            galaxy_server, created = Server.objects.get_or_create(
                url=galaxy_url)
        else:
            try:
                galaxy_server = Server.objects.get(current=True)
            except Server.DoesNotExist:
                raise CommandError(
                    'Server Galaxy does not exist, please use --galaxyurl')

        tools_url = '%s/%s/%s/' % (galaxy_server.url, 'api', 'tools')
        if tool_id:
            # search specific tool
            tools_url_id = '%s/%s/' % (tools_url, tool_id)
            connection = requests.get(tools_url_id).json()
            if 'traceback' in connection:
                raise CommandError('Tool {} was not found in {}'.format(
                    tool_id, galaxy_server.url))
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
            else:
                sys.exit("Cannot connect to Galaxy server")
        if tools_found:
            self.stdout.write("%s" % ('\n'.join(tools_found)))
            response = 'y' if force else raw_input(
                'Do you want (re)import this tool(s)? [y/N] from {}:'.format(galaxy_server))
            if response.lower() == 'y':
                import_tools_report = Tool.import_tools(
                    galaxy_server, tools=tools_found,
                    force=force)
                if import_tools_report['new']:
                    self.stdout.write(
                        self.style.SUCCESS(
                            "%s successfully imported new tools." %
                            (len(import_tools_report['new']))))
                    for ti in import_tools_report['new']:
                        self.stdout.write(
                            self.style.SUCCESS(
                                "%s %s successfully imported" %
                                (ti.name, ti.version))
                        )
                else:
                    self.stdout.write(self.style.WARNING(
                        "No new tool has been imported"))

                for tool_id in import_tools_report['error']:
                    self.stdout.write(
                        self.style.ERROR(
                            "import of %s has failed " % (tool_id))
                    )
        else:
            self.stdout.write("No result was found for query  %s" % (query))

    def associate_flags(self):
        for flaglink in self.flags:
            tool = Tool.objects.filter(
                name=flaglink[0]).first()
            for f in flaglink[1].split(","):
                self.stdout.write(self.style.SUCCESS(
                    "%s : %s" % (flaglink[0], f)))
                flag = ToolFlag.objects.filter(
                    verbose_name=f).first()
                flag.tool.add(tool)

    def set_input_fields(self):
        for input_field in self.input_fields:
            field = ToolInputData.objects.filter(
                tool__name=input_field[0], name=input_field[1]).first()
            if field:
                self.stdout.write(self.style.SUCCESS(
                    "%s : %s" % (input_field[0], input_field[1])))
                field.galaxy_input_data = True
                field.save()
