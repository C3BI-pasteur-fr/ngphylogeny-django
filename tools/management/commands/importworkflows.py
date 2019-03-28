import requests
import re
from django.core.management.base import BaseCommand, CommandError
from slugify import slugify

from galaxy.models import Server
from galaxy.models import GalaxyUser
from workflows.models import Workflow
from tools.models import Tool
from tools.models import ToolFlag


class Command(BaseCommand):
    help = 'Import Galaxy workflows into NGPhylogeny'
    requires_system_checks = True
    flags = []
    wfnames = []

    def add_arguments(self, parser):
        # Named (optional) arguments
        parser.add_argument('--galaxyurl')
        parser.add_argument('--wfnamefile')

    def handle(self, *args, **options):
        wffile = options.get('wfnamefile')
        galaxy_url = options.get('galaxyurl')
        self.read_wfname_file(wffile)
        self.import_workflows(galaxy_url)

    def read_wfname_file(self, wffile):
        with open(wffile) as f:
            for line in f:
                wfname = line.strip()
                self.wfnames.append(wfname)
        
    def import_workflows(self, galaxy_url):
        if galaxy_url:
            galaxy_server, created = Server.objects.get_or_create(
                url=galaxy_url)
        else:
            try:
                galaxy_server = Server.objects.get(current=True)
            except Server.DoesNotExist:
                raise CommandError(
                    'Server Galaxy does not exist, please use --galaxyurl')

        galaxy_key = GalaxyUser.objects.filter(
            galaxy_server__url=galaxy_url,
            anonymous=True)
        print(galaxy_key)
        workflows_url = '%s/%s/%s/?key=%s' % (
            galaxy_server.url,
            'api',
            'workflows',
            galaxy_key.first().api_key)

        # fetch list of tools
        connection = requests.get(workflows_url)
        print(workflows_url)
        print(connection.status_code)
        if connection.status_code == 200:
            wf_list = connection.json() or []
            for wf in wf_list:
                wfname = wf.get('name')
                wfid = wf.get('id')
                self.stdout.write(
                    self.style.SUCCESS(
                        "importing workflow %s" % (wfname)
                    )
                )
                if(re.search('oneclick', wfname, re.IGNORECASE) or
                   wfname in self.wfnames):
                    w = Workflow(
                        galaxy_server=galaxy_server,
                        id_galaxy=wfid,
                        name=wfname,
                        category='base',
                        description=wfname,
                        slug=slugify(wfname))
                    w.save()
        else:
            self.stdout.write("Problem while querying galaxy server")

    def associate_flags(self, flagfile):
        for flaglink in self.flags:
            tool = Tool.objects.filter(
                name=flaglink[0]).first()
            for f in flaglink[1].split(","):
                self.stdout.write(self.style.SUCCESS(
                    "%s : %s" % (flaglink[0], f)))
                flag = ToolFlag.objects.filter(
                    verbose_name=f).first()
                flag.tool.add(tool)
