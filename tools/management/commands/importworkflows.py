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
        parser.add_argument(
            '--wfids',
            help='Comma-separated Galaxy workflow ids to import directly '
                 '(one GET per id, via /api/workflows/{id}) instead of '
                 'listing and filtering the whole server workflow '
                 'collection by name. Precise and fast once you already '
                 'know a target Galaxy\'s stable base-workflow ids (e.g. '
                 'from a database dump) - bypasses the oldest-N heuristic '
                 'import_workflows() otherwise relies on.')

    def handle(self, *args, **options):
        galaxy_url = options.get('galaxyurl')
        wfids = options.get('wfids')
        if wfids:
            ids = [wfid.strip() for wfid in wfids.split(',') if wfid.strip()]
            self.import_workflows_by_id(galaxy_url, ids)
        else:
            self.read_wfname_file(options.get('wfnamefile'))
            self.import_workflows(galaxy_url)

    def read_wfname_file(self, wffile):
        with open(wffile) as f:
            for line in f:
                wfname = line.strip()
                self.wfnames.append(wfname)

    def _resolve_server_and_key(self, galaxy_url):
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
            galaxy_server__url=galaxy_server.url, anonymous=True)
        api_key = galaxy_key.first().api_key
        return galaxy_server, api_key

    def _save_base_workflow(self, galaxy_server, wfid, wfname):
        # Galaxy's workflow list also contains every per-run copy
        # Workflow.duplicate() creates each time a user actually launches
        # this workflow (same name, fresh Galaxy-side id, tracked locally
        # as its own category='duplicated' row) - those must not be
        # confused with the single canonical 'base' copy this command is
        # meant to (re)import. Skip any Galaxy id already tracked locally
        # under any category: this both avoids re-processing already-known
        # duplicates (whose id_galaxy legitimately belongs to a different
        # row and would otherwise collide on the id_galaxy unique
        # constraint) and, for the base row itself, avoids pointlessly
        # overwriting it with whichever same-named entry happens to be
        # processed last.
        if Workflow.objects.filter(
                galaxy_server=galaxy_server, id_galaxy=wfid).exists():
            return
        # update_or_create, not a plain insert: Galaxy re-imports its
        # bundled workflows (fresh id_galaxy each time) on every restart,
        # so re-running this against the same server used to hit the
        # unique constraint on slug with "duplicate key value violates
        # unique constraint workflows_workflow_slug_key" instead of just
        # refreshing the existing row's id_galaxy.
        Workflow.objects.update_or_create(
            galaxy_server=galaxy_server,
            slug=slugify(wfname),
            defaults={
                'id_galaxy': wfid,
                'name': wfname,
                'category': 'base',
                'description': wfname,
            })

    def import_workflows_by_id(self, galaxy_url, wfids):
        """
        Fetch and import specific workflows by their known Galaxy id, one
        GET per id (/api/workflows/{id}) - skips listing/filtering the
        server's whole workflow collection entirely, unlike
        import_workflows() below.
        """
        galaxy_server, api_key = self._resolve_server_and_key(galaxy_url)
        for wfid in wfids:
            wf_url = '%s/api/workflows/%s' % (galaxy_server.url, wfid)
            connection = requests.get(
                wf_url, headers={'x-api-key': api_key})
            if connection.status_code == 200:
                wf = connection.json()
                wfname = wf.get('name')
                self.stdout.write(
                    self.style.SUCCESS(
                        "importing workflow %s (%s)" % (wfname, wfid)))
                self._save_base_workflow(galaxy_server, wfid, wfname)
            else:
                self.stdout.write(
                    self.style.ERROR(
                        "Problem while fetching workflow %s (HTTP %s)" %
                        (wfid, connection.status_code)))

    def import_workflows(self, galaxy_url):
        galaxy_server, api_key = self._resolve_server_and_key(galaxy_url)

        # Real usage over years leaves Galaxy's full /api/workflows/ list
        # dominated by per-run duplicates (Workflow.duplicate() - see
        # workflows/models.py) - 800,000+ of them on galaxy.pasteur.fr,
        # none of them relevant here, since this command only ever
        # imports the canonical "<Tool> OneClick" base workflows. Fetching
        # the whole list on every single redeploy was a heavy, slow
        # request for nothing this command actually needs. Those base
        # workflows are long-lived (created once, e.g. 2019-01-16 on
        # galaxy.pasteur.fr, and never touched again) - they're reliably
        # among Galaxy's OLDEST workflows, not its newest, so
        # sort_by=create_time, sort_desc=false + a small limit finds them
        # without paging through the duplicate backlog at all. See
        # scripts/cleanup_old_galaxy_workflows.sh for the matching
        # oldest-first reasoning on the deletion side, and --wfids above
        # for a more precise/faster alternative once you already know a
        # target Galaxy's stable base-workflow ids.
        workflows_list_limit = 100
        workflows_url = '%s/%s/%s/?key=%s&limit=%d&offset=0&sort_by=create_time&sort_desc=false' % (
            galaxy_server.url,
            'api',
            'workflows',
            api_key,
            workflows_list_limit)

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
                    self._save_base_workflow(galaxy_server, wfid, wfname)
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
