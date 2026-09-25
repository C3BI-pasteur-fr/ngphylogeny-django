# NGPhylogeny.fr

A Django web application for phylogenetic analysis. NGPhylogeny.fr itself
does no bioinformatics computation — it's an orchestration/UI layer that
drives a separate [Galaxy](https://galaxyproject.org/) server (via the
`bioblend` API client) to run the actual tools (alignment, tree building,
bootstrapping, BLAST, etc.) and tracks job status/results locally.

# Installation

NGPhylogeny.fr requires python 3.8 (see the Dockerfile).

To install it:

* You may first install and activate a new conda environment:

```
conda create --name ngphylo python=3.8
source activate ngphylo
```

* Install all NGPhylogeny.fr dependencies:

```
pip install -r requirement.txt
apt-get install redis-server # Ubuntu / Debian
```

* Create Django databases

```
python manage.py makemigrations
python manage.py migrate
python manage.py createcachetable
python manage.py collectstatic
```

* Create admin user

```
python manage.py createsuperuser
```

* Load tool categories (flags)

```
python manage.py loaddata tool_flags
```

* Link NGPhylogeny.fr to a running Galaxy server

```
python manage.py creategalaxyserver --url=http://url_galaxy:port --activate
python manage.py addgalaxykey --user <superuser name> --galaxyurl http://url_galaxy:port --galaxykey <galaxy key>
```

* Import existing tools from Galaxy to NGPhylogeny.fr

```
python manage.py importtools --galaxyurl=http://url_galaxy:port \
                             --query="phylogeny" \
                             --flags=toolflags.txt \
                             --inputfields=toolfields.txt
```

* Two possibilities to link tools inputs and outputs

  1. Already computed links (preferable):

  ```
  python manage.py import_links --linkfile=toollinks.txt
  ```

  2. Links computed on the fly (based on extension compatiblity, deprecated)

  ```
  python manage.py compute_tools_links --ignore=txt
  ```
* Import Phylogeny workflows from the Galaxy Instance

It will import all workflows with name containing "oneclick" (case insensitive):
```
python manage.py importworkflows --galaxyurl=http://url_galaxy:port --wfnamefile=wfnames.txt
```

* Run Celery task queue

	1. Add following lines to `NGPhylogeny_fr/settings/local.py`:
	```
	EMAIL_HOST = 'smtp.server.url'
	EMAIL_PORT = 587
	EMAIL_HOST_USER = 'smtp.user'
	EMAIL_HOST_PASSWORD = 'smtp.pass'
	EMAIL_USE_TLS=<True|False>
	```

	2. Start celery
	```
	export PYTHONPATH=$PWD:$PYTHONPATH
	celery multi start 3 -l INFO -c:2 1 -c:3 1 -Q:1 default -Q:2 ncbi_blast -Q:3 monitor --app=NGPhylogeny_fr.celery:app
	celery beat --app=NGPhylogeny_fr.celery:app --loglevel=DEBUG --detach
	```

* Run the django server

```
python manage.py runserver
```

# Tests and linting

GitLab CI (`.gitlab-ci.yml`) runs both of these on every push/merge request:

```
python manage.py makemigrations
python manage.py migrate
python manage.py check
python manage.py test

flake8 --select=E9,F63,F7,F82 --exclude=migrations .
```

# Docker

```
docker compose up -d
```

This runs five services (see `docker-compose.yml`): Postgres, Redis, the
Django app (`web`, on http://localhost:8000), and a Celery worker + beat.
A one-shot `init` service runs migrations, loads the tool-flag fixtures, and
creates an admin user (`admin` / `password` by default - override via the
`NGPHYLO_ADMIN_USER`/`NGPHYLO_ADMIN_EMAIL`/`NGPHYLO_ADMIN_PASSWORD` env vars)
before the other services start.

The app works without one, but to actually run analyses it needs a Galaxy
server linked (with NGPhylogeny's tools/workflows installed - see the
[NGPhylogeny_fr_galaxytools](https://github.com/C3BI-pasteur-fr/ngphylogeny-galaxy)
repo, which has its own `docker compose up` deployment). Point this stack at
one by setting `NGPHYLO_GALAXY_URL` and `NGPHYLO_GALAXY_KEY` (a Galaxy admin
API key) before starting it:

```
NGPHYLO_GALAXY_URL=http://host.docker.internal:8080 \
NGPHYLO_GALAXY_KEY=<galaxy admin api key> \
docker compose up -d
```

`init` only links Galaxy and imports its tools/workflows on that first run;
to do it later against an already-running stack, `docker compose run init`
(or the underlying `manage.py creategalaxyserver`/`addgalaxykey`/`importtools`/
`importworkflows` commands - see CLAUDE.md) works too.

By default the app runs with `DEBUG=True` (`NGPhylogeny_fr.settings.local`)
so Django's dev server serves static files itself with no extra setup;
override `NGPHYLO_SETTINGS_MODULE=NGPhylogeny_fr.settings.prod` for
production-like behavior, but then something needs to serve `STATIC_ROOT`
(`/static/`) separately - this compose file doesn't set that up.

## Maintenance mode and BLAST server selection

`NGPHYLO_MAINTENANCE_MODE=True` serves `templates/maintenance.html` (503)
for every request instead of routing normally (see
`NGPhylogeny_fr/middleware.py`). Off by default; safe to flip on/off
without a rebuild since it's just an env var.

`NGPHYLO_PASTEUR_BLAST_ENABLED=True` activates the Institut Pasteur Galaxy
BLAST server option and deactivates NCBI's public one at the same time -
the two aren't independently toggleable (see `settings/base.py`'s `BLASTS`
dict). Off by default, meaning BLAST runs go to NCBI's public server.

`NGPHYLO_ACCOUNT_CREATION_ENABLED=True` re-enables public account sign-up
(`/account/create`) - off by default pending a real RGPD/privacy notice
(see `account/views.py`'s `AccountCreateView`); left off, that page serves
a 503 "temporarily unavailable" page instead.

`NGPHYLO_WORKSPACE_RETENTION_DAYS` sets how many days a finished analysis
is kept before `workspace.tasks.deleteoldgalaxyhistory`'s daily cleanup
removes it (`workspace/models.py`'s `WorkspaceHistory.RETENTION_DAYS`) -
also what the "days left" estimate on the Workspace/account pages is
computed from. Defaults to `14` if unset.

A few related, independently configurable cleanup/staleness cutoffs:
- `NGPHYLO_BLAST_RETENTION_DAYS` - days a BLAST run is kept before
  `blast.tasks.deleteoldblastruns`'s daily cleanup removes it
  (`BlastRun.RETENTION_DAYS`). Defaults to `7`.
- `NGPHYLO_WORKFLOW_RETENTION_DAYS` - days a non-base (per-run
  duplicated) Galaxy workflow *definition* is kept before
  `workflows.tasks.deleteoldgalaxyworkflows`'s daily cleanup removes it
  (`Workflow.RETENTION_DAYS`) - doesn't affect the history's own actual
  data. Defaults to `14`.
- `NGPHYLO_WORKFLOW_RUN_STALE_HOURS` - hours a still-running/queued
  analysis is left alone before it's force-cancelled
  (`workspace/tasks.py`'s `WORKFLOW_RUN_STALE_AFTER`). Defaults to `24`.
- `NGPHYLO_PASTEUR_BLAST_STALE_HOURS` - hours a still-running Pasteur
  BLAST search is polled before it's given up on and marked `ERROR`
  (`blast/tasks.py`'s `PASTEUR_RUN_STALE_AFTER`). Defaults to `3`.

## Email: job-completion notices and the daily report

Both job-completion emails and the daily workflow-usage report (HTML, with
charts - `workspace/reports.py`, sent by `workspace.tasks.send_daily_report`
at 8am UTC every day, see `CELERY_BEAT_SCHEDULE` in `settings/base.py`) need
working SMTP config: `NGPHYLO_EMAIL_HOST`/`NGPHYLO_EMAIL_PORT`/
`NGPHYLO_EMAIL_HOST_USER`/`NGPHYLO_EMAIL_HOST_PASSWORD`/`NGPHYLO_EMAIL_USE_TLS`.
Left unset (the default), sending just fails silently - the app itself is
unaffected either way.

The daily report additionally needs `NGPHYLO_REPORT_RECIPIENTS` (comma-
separated email addresses) set - without it, `send_daily_report` logs and
does nothing, so it's safe to leave the periodic task enabled on
deployments that don't want the report. `NGPHYLO_REPORT_FROM_EMAIL`
overrides the sender address (defaults to `ngphylogeny@pasteur.fr`).

The contact form (`surveys` app) similarly needs
`NGPHYLO_CONTACT_FORM_RECIPIENTS` (comma-separated) set to also email a
submitted message to staff, in addition to saving it as a `Feedback` row -
left unset, only the DB row is saved, same no-op-by-default convention as
the daily report above.

## Standalone (Django + Galaxy, all in one)

`docker-compose.yml` above needs a Galaxy server to already be running
somewhere. To test NGPhylogeny.fr completely from scratch - Django, Postgres,
Redis, *and* a Galaxy server with NGPhylogeny's own tools/workflows, all
brought up and wired together by one command - use
`docker-compose.standalone.yml` instead. This is a separate, self-contained
stack, not an overlay for the file above: run one or the other, don't combine
them with `-f`.

It needs the [NGPhylogeny_fr_galaxytools](https://github.com/C3BI-pasteur-fr/ngphylogeny-galaxy)
repo cloned as a sibling directory (`../NGPhylogeny_fr_galaxytools` by
default - override with `GALAXYTOOLS_DIR` if you keep it elsewhere), since
that's where the tool wrappers, the PhyML-SMS/Noisy combined-image
Dockerfiles, and the base workflow `.ga` files it needs actually live:

```
git clone https://github.com/C3BI-pasteur-fr/ngphylogeny-galaxy.git ../NGPhylogeny_fr_galaxytools
docker compose -f docker-compose.standalone.yml up -d
```

That one command chains everything through `depends_on` conditions: Postgres
and Galaxy come up, two one-shot services then build the PhyML-SMS/Noisy
combined images inside Galaxy's own internal Docker daemon and import the 4
base "\<Tool\> OneClick" workflows, and only then does NGPhylogeny's own
`init` (migrations, Galaxy linking, tool/workflow import) run - `web`/
`celery-*` wait for that. First boot takes a few minutes (pulling the Galaxy
image, conda/container tool dependency resolution on first use, etc.).
NGPhylogeny comes up on http://localhost:8000, Galaxy itself on
http://localhost:8080 (not needed for normal use, but there if you want to
poke at it directly - `admin` / `password` by default, same as NGPhylogeny's
own admin login, overridable via the same `NGPHYLO_ADMIN_*`/`GALAXY_ADMIN_*`
env vars documented at the top of the file).

Re-running `docker compose -f docker-compose.standalone.yml up -d` (without
`down -v` first) is safe: the workflow import step skips workflows that
already exist by name, and NGPhylogeny's own `init` is idempotent for the
same reasons the regular `docker-compose.yml` deployment is (see CLAUDE.md).

## Kubernetes (GitLab CI)

`.gitlab-ci.yml`'s `build`/`deploy-dev`/`deploy-prod` stages build this
same `Dockerfile` image and deploy it to a Kubernetes cluster via
`kubectl apply` of `manifest_datastores.yaml` (Postgres, Redis) and
`manifest.yaml` (the app itself - a one-shot init `Job`, `web`/
`celery-worker`/`celery-beat` Deployments, a Service, and an Ingress),
env-var-templated via `envsubst`. Galaxy itself is external - point
`NGPHYLO_GALAXY_URL`/the Secret's `galaxy-key` at an already-running
Galaxy server (the Pasteur Galaxy server, for this deployment), same as
every other deployment path in this repo.

`deploy-dev` is live and confirmed working end to end against
`ngphylogenyfr-dev` (a real OneClick workflow has been submitted and
completed there); `deploy-prod`'s namespace/domain are still a
placeholder pending real provisioning. See CLAUDE.md's "Kubernetes
deployment" section for the full detail, including what each manifest
does, the required GitLab CI/CD variables, and a list of real-cluster-
only issues (RBAC, runner selection, the migrations-vs-a-persistent-
database trap) hit getting `deploy-dev` working - worth reading before
setting up `deploy-prod`, since the same class of issues will likely
recur there.

`scripts/cleanup_old_galaxy_workflows.sh` is a standalone maintenance
script (not run by any pipeline) for bulk-deleting old, non-base
workflows directly from a Galaxy server's own `/api/workflows` - see its
own header comment and CLAUDE.md's "Workflow duplicates and the Celery
cleanup jobs" section for why this was needed (832,000+ accumulated rows
on `galaxy.pasteur.fr`) and how it's meant to be used.
