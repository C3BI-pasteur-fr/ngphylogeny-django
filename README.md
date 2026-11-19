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
`NGPHYLO_MAINTENANCE_ALTERNATE_URL` (optional) adds a box on that page
pointing visitors at a temporary instance (on Kubernetes, the
`MAINTENANCE_ALTERNATE_URL` GitLab CI/CD variable). Left empty, the page
shows only its generic message.

`NGPHYLO_PASTEUR_BLAST_ENABLED=True` activates the Institut Pasteur Galaxy
BLAST server option and deactivates NCBI's public one at the same time -
the two aren't independently toggleable (see `settings/base.py`'s `BLASTS`
dict). Off by default, meaning BLAST runs go to NCBI's public server.

`NGPHYLO_ACCOUNT_CREATION_ENABLED=True` re-enables public account sign-up
(`/account/create`) *and* the "forgot your password?" flow
(`/account/password-reset`) - both off by default for now (see
`account/views.py`'s `AccountCreationGateMixin`); left off, those pages
serve a 503 "temporarily unavailable" page instead, and the login page
hides both links.

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

`NGPHYLO_POW_CAPTCHA_DIFFICULTY` tunes the self-hosted proof-of-work
CAPTCHA on the contact form and account sign-up (`utils/powcaptcha.py` -
replaces django-simple-captcha's distorted-letter image, with no
external service/API key). It's the number of leading hex-zero digits a
solution's SHA-256 hash must have; each extra digit multiplies the
browser's expected search by 16x. Defaults to `5` (~1e6 expected
attempts, measured at well under a second up to a few seconds in a real
browser - see that module's docstring for how this was calibrated).

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

`NGPHYLO_DEFAULT_FROM_EMAIL` sets the sender address for anything sent
without its own explicit from address - currently just the password-reset
email (`django.contrib.auth`'s own `PasswordResetForm.save()`, no other
hook for this). Defaults to `ngphylogeny@pasteur.fr`, same as
`NGPHYLO_REPORT_FROM_EMAIL`'s own default; left unset, Django's own default
(`webmaster@localhost`) would otherwise leak into the From header.

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

### Standalone as a temporary production instance

The standalone stack can serve real traffic while the production Galaxy is
down (done on 2026-11-17, see CLAUDE.md). It needs an overlay next to
`docker-compose.standalone.yml` on the server, plus a `.env` and a
`Caddyfile`. Only ports 443 (and 80, if open) are public; Django and Galaxy
listen on localhost only.

`docker-compose.prod.yml`:

```yaml
x-prod-env: &prod-env
  DJANGO_SETTINGS_MODULE: NGPhylogeny_fr.settings.prod
  NGPHYLO_SECRET_KEY: "${NGPHYLO_SECRET_KEY}"
  NGPHYLO_HTTPS_HOST: "${NGPHYLO_HTTPS_HOST}"

services:
  db:            {restart: unless-stopped}
  redis:         {restart: unless-stopped}
  galaxy:
    restart: unless-stopped
    ports: !override
      - "127.0.0.1:8080:80"
  init:          {environment: *prod-env}
  celery-worker: {environment: *prod-env, restart: unless-stopped}
  celery-beat:   {environment: *prod-env, restart: unless-stopped}
  web:
    environment: *prod-env
    restart: unless-stopped
    command: ["uwsgi", "--http", ":8000", "--module", "NGPhylogeny_fr.wsgi", "--master",
              "--processes", "8", "--threads", "2", "--harakiri", "120", "--die-on-term"]
    ports: !override
      - "127.0.0.1:8000:8000"
  caddy:
    image: caddy:2
    restart: unless-stopped
    ports: ["80:80", "443:443"]
    volumes:
      - ./Caddyfile:/etc/caddy/Caddyfile:ro
      - caddy-data:/data
      - caddy-config:/config
    depends_on: [web]

volumes:
  caddy-data:
  caddy-config:
```

`Caddyfile` (the host name must resolve to the server, for Let's Encrypt;
port 443 alone is enough to get the certificate):

```
<hostname> {
    request_body {
        max_size 60MB
    }
    reverse_proxy web:8000
}
```

`.env` holds random values for `NGPHYLO_ADMIN_PASSWORD`,
`GALAXY_ADMIN_PASSWORD`, `GALAXY_ADMIN_KEY`, `GALAXY_BOOTSTRAP_KEY` and
`NGPHYLO_SECRET_KEY`, plus `NGPHYLO_HTTPS_HOST=<hostname>`. The defaults in
`docker-compose.standalone.yml` (`admin`/`password`, fixed keys) are for
local testing only. Copy production's email, retention and account
settings into it as well, and set `NGPHYLO_PASTEUR_BLAST_ENABLED=False`
if the Pasteur Galaxy is the one that's down, so BLAST goes to NCBI.

```
docker compose -f docker-compose.standalone.yml -f docker-compose.prod.yml up -d --build
```

Then, before announcing it:

- **Copy the tool size limits** (`max_nbseq`, `max_boot`,
  `max_length_x_nbseq`) from production: a fresh tool import sets them all
  to `-1`, unlimited.
- **Warm up every tool, one at a time**: run each OneClick workflow and each
  single tool once, sequentially. Galaxy installs each tool's dependencies on
  first use, which takes several minutes, and simultaneous first runs can
  leave a broken conda environment behind (see CLAUDE.md).
- **Add a host cron to clean up Galaxy's internal Docker**: this stack has
  no `docker-prune` sidecar, and every image build leaves its containers
  behind. Every 30 minutes, e.g. in `/etc/cron.d/galaxy-container-prune`:
  ```
  */30 * * * * root docker exec ngphylo-standalone-galaxy sh -c "docker container prune -f --filter until=10m && docker image prune -f" >> /var/log/galaxy-container-prune.log 2>&1
  ```
  Never `image prune -a`: it would delete `local/phyml-sms-full` and
  `local/noisy-full`, which are built on the server and can't be pulled.

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
completed there). `deploy-prod` serves `ngphylogeny.fr` from
`ngphylogenyfr-prod`, runs from `master`, and only starts when triggered
by hand from the pipeline page. See CLAUDE.md's "Kubernetes deployment"
section for the full detail, including what each manifest does, the
required GitLab CI/CD variables, and a list of real-cluster-only issues
(RBAC, runner selection, the migrations-vs-a-persistent-database trap)
worth reading before setting up another deployment.

`/status` shares its worker slots with every other request, so under
heavy load a healthy `web` pod can be slow to answer it. Its probes are
patient on purpose: readiness allows about 60s of failures, liveness
about 120s, matching uwsgi's `--harakiri`. There's only one `web` replica,
so tightening either one turns a slow site into an unavailable one. uwsgi
runs with `--die-on-term`, so a rollout stops the old pod straight away
rather than waiting out the 30s grace period.

Uploads are capped at 50 MiB (`NGPHYLO_MAX_UPLOAD_MB` in
`settings/base.py`). Django enforces this, and `proxy-body-size` in
`manifest.yaml` should be kept in step with it as a backstop.

If Galaxy is unreachable, pages that need it return a 503 "Galaxy
unavailable" page, and background polling pauses for 90s at a time. No
action is needed once Galaxy is back. A run's status on its history page
can lag by up to about 50s, because each run is checked at most every 30s
to limit the load on Galaxy (every 9s until its steps first appear).

The footer on every page shows the running image's git commit, linked
to its page on GitHub - `Dockerfile` bakes it in at build time
(`--build-arg GIT_COMMIT`, which `.gitlab-ci.yml`'s `build` job passes
as `$CI_COMMIT_SHA`) into a plain `GIT_COMMIT` file, since `.dockerignore`
excludes `.git` itself from the build context. A local `docker compose
build` (no CI context) leaves this unset, and the footer just shows no
version line rather than a broken link.

`scripts/cleanup_old_galaxy_workflows.sh` is a standalone maintenance
script (not run by any pipeline) for bulk-deleting old, non-base
workflows directly from a Galaxy server's own `/api/workflows` - see its
own header comment and CLAUDE.md's "Workflow duplicates and the Celery
cleanup jobs" section for why this was needed (832,000+ accumulated rows
on `galaxy.pasteur.fr`) and how it's meant to be used.
