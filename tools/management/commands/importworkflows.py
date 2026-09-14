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
    requires_system_checks = '__all__'
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
        api_key = galaxy_key.first().api_key
        workflows_url = '%s/%s/%s/?key=%s' % (
            galaxy_server.url,
            'api',
            'workflows',
            api_key)

        # This is a raw (non-bioblend) Galaxy call, like the ones in
        # data/views.py - bioblend>=1.4.0 moved auth from ?key= query
        # params to an x-api-key header automatically, but this command
        # predates that and never picked it up. Against a permissive
        # Galaxy (e.g. a local dev instance) the query param alone is
        # enough and this silently worked; against galaxy.pasteur.fr it
        # doesn't - this request just gets a non-200 back, which this
        # command already handles by logging "Problem while querying
        # galaxy server" and returning without importing anything or
        # raising, so importworkflows "succeeds" with zero workflows
        # imported and no error in the init Job's own exit code. Send
        # both, same as data/views.py's fix for the identical issue.
        connection = requests.get(workflows_url, headers={'x-api-key': api_key})
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
                    # Galaxy's workflow list also contains every per-run
                    # copy Workflow.duplicate() creates each time a user
                    # actually launches this workflow (same name, fresh
                    # Galaxy-side id, tracked locally as its own
                    # category='duplicated' row) - those must not be
                    # confused with the single canonical 'base' copy this
                    # command is meant to (re)import. Skip any Galaxy id
                    # already tracked locally under any category: this
                    # both avoids re-processing already-known duplicates
                    # (whose id_galaxy legitimately belongs to a different
                    # row and would otherwise collide on the id_galaxy
                    # unique constraint) and, for the base row itself,
                    # avoids pointlessly overwriting it with whichever
                    # same-named entry happens to be listed last.
                    if Workflow.objects.filter(
                            galaxy_server=galaxy_server,
                            id_galaxy=wfid).exists():
                        continue
                    # update_or_create, not a plain insert: Galaxy re-imports
                    # its bundled workflows (fresh id_galaxy each time) on
                    # every restart, so re-running this against the same
                    # server used to hit the unique constraint on slug with
                    # "duplicate key value violates unique constraint
                    # workflows_workflow_slug_key" instead of just
                    # refreshing the existing row's id_galaxy.
                    w, created = Workflow.objects.update_or_create(
                        galaxy_server=galaxy_server,
                        slug=slugify(wfname),
                        defaults={
                            'id_galaxy': wfid,
                            'name': wfname,
                            'category': 'base',
                            'description': wfname,
                        })
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
