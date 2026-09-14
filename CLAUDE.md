# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

NGPhylogeny.fr is a Django front-end for phylogenetics analysis. **It does no
bioinformatics computation itself** — it's an orchestration/UI layer that
drives a separate **Galaxy server** (via the `bioblend` API client) to run
the actual tools (alignment, tree building, bootstrapping, BLAST, etc.). Any
change to a workflow, tool run, or job-status feature usually means tracing
through both this Django app *and* what it's asking the Galaxy server to do
— there's no local execution to fall back on for understanding behavior.

## Commands

Stack: Python 3.8, Django 4.2, Celery 5 + Redis.

```bash
# Install deps (numpy must land before biopython, which needs it at build time
# and doesn't declare it — a plain `pip install -r requirement.txt` can fail
# on a clean env depending on pip's resolution order)
pip install numpy==1.24.4
pip install -r requirement.txt

# Local dev DB: sqlite unless NGPHYLO_DATABASE_HOST is set (see settings/local.py)
python manage.py makemigrations   # migrations are NOT committed to git - see below
python manage.py migrate
python manage.py createcachetable
python manage.py runserver

# Tests (currently: static-page + crispy/captcha rendering smoke tests in
# data/tests.py, Tool.can_run_on_data() boundary tests in tools/tests.py)
python manage.py test
python manage.py test data.tests.StaticPagesSmokeTest.test_pages_return_200   # single test

# Lint (matches CI exactly - narrow on purpose, this is a legacy codebase)
flake8 --select=E9,F63,F7,F82 --exclude=migrations .

# DJANGO_SETTINGS_MODULE defaults to NGPhylogeny_fr.settings.local (see manage.py);
# production (wsgi.py) uses NGPhylogeny_fr.settings.prod instead.
```

GitLab CI (`.gitlab-ci.yml`) runs exactly `flake8` then
`makemigrations && migrate && check && test` on `python:3.8-buster`. There's
no Docker build/push/deploy stage yet.

### One-time / operational management commands

These aren't used in normal dev loops but are how the app gets linked to a
Galaxy server and populated with tools/workflows (see `docker/init.sh` for
the full bootstrap sequence used in Docker):

```bash
python manage.py creategalaxyserver --url=http://url_galaxy:port --activate
python manage.py addgalaxykey --user <username> --galaxyurl <url> --galaxykey <key>
python manage.py importtools --galaxyurl=<url> --query="phylogeny" --flags=toolflags.txt --inputfields=toolfields.txt
python manage.py import_links --linkfile=toollinks.txt
python manage.py importworkflows --galaxyurl=<url> --wfnamefile=wfnames.txt
```

## Architecture

### Migrations are not committed

`.gitignore` excludes `*/migrations/0*.py` (only `migrations/__init__.py` is
tracked). `docker/init.sh` runs `makemigrations` fresh on every container
start (see "Docker" below). Don't expect `git log` on a migrations
directory to tell you anything, and always run `makemigrations` before
`migrate` in a fresh checkout.

### App responsibilities

- **`galaxy`** — `Server` (a configured Galaxy instance; `current=True` marks
  the one actually used — enforced as a singleton in `Server.save()`) and
  `GalaxyUser` (per-Django-user Galaxy API key, with an `anonymous=True` flag
  marking the shared key used for anonymous visitors). `galaxylib.py` wraps
  `bioblend.galaxy.GalaxyInstance`, adding Django-cache-backed GET caching
  and a `GalaxyInstanceAnonymous` variant that authenticates via a Galaxy
  session cookie instead of an API key (for anonymous users).
- **`tools`** — mirrors Galaxy's tool definitions locally (`Tool`,
  `ToolInputData`/`ToolOutputData` with EDAM format info, `Citation`,
  `ToolFlag`). `Tool.import_tool_io()` fetches live JSON from the Galaxy API
  to populate these on import. `Tool.can_run_on_data()` enforces
  per-tool size limits (max sequences, bootstrap replicates, etc.) before a
  job is allowed to run.
- **`workflows`** — `Workflow` mirrors a Galaxy workflow. The interesting
  piece is `WorkflowGalaxyFactory` (in `models.py`): given an ordered list of
  `Tool`s, it auto-chains them into Galaxy workflow JSON by matching
  EDAM input/output formats step-to-step, inserting a Galaxy
  `ChangeDatatypeAction` when formats don't quite match. Three submission
  UIs share this machinery: `views/wkoneclick.py` (preconfigured),
  `views/wkadvanced.py` (parametrized), `views/wkmaker.py` (build-your-own,
  drives `WorkflowGalaxyFactory` directly).
- **`workspace`** — `WorkspaceHistory` represents one Galaxy "history" (a
  job run): status, ownership, monitoring state. Celery tasks
  (`workspace/tasks.py`) poll Galaxy for job status and email the user on
  completion.
- **`blast`** — a parallel, mostly self-contained BLAST subsystem (own
  models, Celery tasks, pseudo-MSA construction in `msa.py`) supporting
  either the public NCBI API or an internal Pasteur Galaxy BLAST server
  (`settings.BLASTS` in `settings/base.py` configures both).
- **`data`** — `ExampleFile` plus generic file display/download/upload
  views shared across the other apps.
- **`account`**, **`surveys`** — thin (`account` has no models; it's
  login/API-key management. `surveys` is the feedback form).

### Async: Celery

Queues: `default`, `ncbi_blast`, `monitor` (routed via `CELERY_TASK_ROUTES`
in `settings/base.py`). Periodic jobs (old-history cleanup, BLAST run
monitoring, etc.) are declared in `CELERY_BEAT_SCHEDULE` in that same file —
**not** as `@periodic_task` decorators on the task functions (that API was
removed in Celery 5; the task functions are plain `@shared_task`). All
`CELERY_*` settings are read via
`app.config_from_object('django.conf:settings', namespace='CELERY')` in
`NGPhylogeny_fr/celery.py`, so a setting only takes effect if it has the
`CELERY_` prefix and the correct post-namespace name (e.g.
`CELERY_TASK_DEFAULT_QUEUE`, not `CELERY_DEFAULT_QUEUE`).

### Settings split

`NGPhylogeny_fr/settings/{base,local,prod}.py`. `DATABASES` lives in
`base.py` (env-var-based: `NGPHYLO_DATABASE_*` if `NGPHYLO_DATABASE_HOST` is
set, sqlite fallback otherwise) so both `local.py` and `prod.py` inherit the
same logic via `from .base import *`. `local.py` (the default — see
`manage.py`) additionally sets `DEBUG=True`; `prod.py` (used by `wsgi.py`)
sets `DEBUG=False`. Custom error templates (`templates/500.html` etc.) only
ever render when `DEBUG=False` — Django always shows the interactive
traceback otherwise, regardless of what templates exist, so reproducing a
production-looking error page locally means running with
`DJANGO_SETTINGS_MODULE=NGPhylogeny_fr.settings.prod`, not `.local`.

### `dictsort` doesn't call methods

Django 3.1 hardened `dictsort`/`dictsortreversed` to stop auto-calling
methods in a lookup chain (`"foo.bar.baz"` where an intermediate step is a
method, e.g. a related-manager's `.first()`) — unlike normal template
variable resolution, which still auto-calls. A filter argument like
`dictsort:"toolflag_set.first.verbose_name"` will now silently resolve to
`""` instead of raising, so a `{% regroup %}`/`{% for %}` over it just
renders as empty — no error, no page-load failure, easy to miss entirely
unless the list actually has data (`tools/views.py`'s `ToolListView` hit
this rendering as a permanently-empty tools page; see git history for the
fix — precompute the value as a plain attribute in the view instead of
relying on dictsort to call anything).

### Code paths only a real Galaxy run exercises

A full local Galaxy setup (any current version — verified against Galaxy
25.1, tools/workflows from `NGPhylogeny_fr_galaxytools` installed via its
`tool_conf.xml` + conda auto-install) surfaced three more Python 2→3 bugs
that `manage.py check`/`test` and synthetic bioblend calls never touch,
because they only trigger when a real workflow actually runs end to end:
- `requires_system_checks = True` (a bare bool) in `creategalaxyserver`,
  `addgalaxykey`, `importtools`, `importworkflows` — Django 4.x requires a
  list/tuple or the `'__all__'` sentinel here, not a bool. These commands
  are how a Galaxy server gets linked to the app at all (see the commands
  above), so this broke silently through every earlier phase.
- `Tool.import_tools()`'s citation text used to round-trip through
  `.encode('iso-8859-1').decode('utf8')` — a Python 2-era fix for a
  `requests`/`json` mojibake quirk that doesn't exist under Python 3
  (`requests.json()` already returns correctly-decoded text). Re-applying
  it crashed on any citation with a character outside Latin-1 (en dashes,
  curly quotes — both common in real citations), blocking most tool
  imports.
- `utils/biofile.py`'s `valid_fasta()` passed an uploaded file's binary
  handle straight to `Bio.SeqIO.parse()`. Biopython's `SimpleFastaParser`
  detects EOF by comparing a line to `""` (str), which never matches
  `b""` (bytes) under Python 3 — so parsing any real uploaded fasta file
  crashed with `IndexError` at the true end of the file, every time. Only
  worked in Python 2, where `bytes == str`. Fixed by reading the upload
  fully upfront and parsing through a text `io.StringIO` regardless of
  whether the source yields bytes or str.

With all three fixed, a real "PhyML" oneclick workflow submitted through
the actual UI ran for real against the local Galaxy: upload succeeds,
MAFFT/seqtype-detection/BMGE/PhyML execute as real (conda-installed)
Galaxy jobs. If you're touching upload/import/tool-linking code, don't
trust `manage.py test` alone — these paths need an actual Galaxy to run
against.

### Known dependency ceilings (don't casually bump these)

- **`biopython==1.70`** is pinned and can't be bumped past Python 3.9: its
  `Bio.trie` C extension assigns through the `Py_TYPE()` macro, which CPython
  3.10 made illegal. Bumping biopython also means removing all
  `Bio.Alphabet` usage (`blast/msa.py`, `utils/biofile.py`,
  `blast/tasks.py`) since that module was removed in biopython ≥1.78. This
  is why the project is pinned to Python 3.8 rather than something newer.
- **`celery[redis]==5.4.0`**, **`bioblend==1.4.0`** were bumped from very old
  pins (4.4.7, 0.10) as part of this project's Python 2→3 / Django
  1.11→4.2 migration. `bioblend`'s runtime behavior has since been verified
  live both against the real production Galaxy server (`galaxy.pasteur.fr`
  — auth moved from `?key=` query params to an `x-api-key` header,
  `get_tools`/`get_histories`, `Workflow.fetch_details()`/`show_workflow()`,
  `Workflow.duplicate()`/`import_workflow_dict`) and end to end against a
  local Galaxy 25.1 running a real oneclick workflow (see "Code paths only
  a real Galaxy run exercises" above) — both confirmed working.

### Docker

`docker compose up -d` runs five services from one lean `Dockerfile` (`db`
postgres, `redis`, `web`, `celery-worker`, `celery-beat` — no
nginx/uwsgi/compiled-redis bundled into the image, unlike the pre-2026
deployment). A one-shot `init` service (`docker/init.sh`) runs migrations,
seeds the admin user, and — only if `NGPHYLO_GALAXY_URL`/`NGPHYLO_GALAXY_KEY`
are set — links a Galaxy server and imports its tools/workflows, before
`web`/`celery-*` start (`depends_on: init: condition:
service_completed_successfully`); see that script for the exact command
sequence if reproducing it outside Docker. Verified end to end against a
from-scratch `docker compose up` — see README.md for deploy steps and the
Galaxy-linking env vars.

Two dockerignore gotchas worth knowing if this ever breaks again: unlike
`.gitignore`, a bare pattern in `.dockerignore` (`*.pyc`) only matches at the
build context *root*, not at any depth — nested files need the `**/` prefix
(`**/*.pyc`). And migrations/__init__.py must stay in the build context (only
the generated `0*.py` files are excluded, matching `.gitignore`) — excluding
the whole `*/migrations` directory silently breaks `makemigrations`'
auto-detection for every app that hasn't got a migration yet.

`docker-compose.standalone.yml` is a separate, self-contained alternative
(not an overlay — don't combine the two with `-f`) that also brings up a
Galaxy server itself, for testing entirely from scratch with no pre-existing
Galaxy needed. It requires `NGPhylogeny_fr_galaxytools` cloned as a sibling
directory (`../NGPhylogeny_fr_galaxytools`, override via `GALAXYTOOLS_DIR`)
and chains, via `depends_on` conditions, straight through: Postgres/Galaxy
come up → `galaxy-build-images` (a `docker:27-cli` sidecar with the host
Docker socket mounted, reaching into Galaxy's own container to build the
PhyML-SMS/Noisy combined images by reusing that repo's own
`docker/build-combined-images.sh`) and `galaxy-import-workflows`
(`docker/import_base_workflows.py` in this repo, using bioblend — already a
dependency here — to import the 4 base `.ga` workflows) → this repo's own
`init`. See README.md's "Standalone" section for usage.

### `NGPHYLO_SETTINGS_MODULE` vs `DJANGO_SETTINGS_MODULE` in a prod override

A production override (`docker-compose.prod.yml`, server-side only — see
IFB_CLOUD.md) needs to set `DJANGO_SETTINGS_MODULE` directly, **not**
`NGPHYLO_SETTINGS_MODULE`. The base `docker-compose.yml`'s
`DJANGO_SETTINGS_MODULE: "${NGPHYLO_SETTINGS_MODULE:-NGPhylogeny_fr.settings.local}"`
looks like it's meant to be overridden by setting `NGPHYLO_SETTINGS_MODULE`
in a second `-f` file's `environment:` block — it isn't. That `${...}`
substitution is resolved from the **host's** shell/`.env` at
`docker compose` config-parse time, before any file's `environment:` map
is merged; an override file setting `NGPHYLO_SETTINGS_MODULE` there just
adds an unused env var under that name to the container, while
`DJANGO_SETTINGS_MODULE` itself silently stays whatever the base file's
substitution resolved to (`local` — dev settings, `DEBUG=True` — unless
the host/`.env` happens to also define `NGPHYLO_SETTINGS_MODULE`, which
it won't by default). This ran the real IFB Cloud deployment in `DEBUG=True`
for its entire lifetime, undetected, until it broke something that only
manifests under `settings.prod` (`NGPHYLO_HTTPS_HOST` doesn't even exist
as a setting under `settings.local`, so anything reading it — e.g.
`workspace/emails.py`'s job-completion email link — silently fell back to
a wrong default host instead of raising). Fix in the override: set
`DJANGO_SETTINGS_MODULE: NGPhylogeny_fr.settings.prod` directly.

### `collectstatic` needs a volume shared between `init` and `web`

`docker/init.sh` runs `collectstatic` inside the one-shot `init`
container — but `init` and `web` are separate containers from the same
image, each with their own independent, ephemeral filesystem. Without a
shared volume for `STATIC_ROOT` (`/home/ngphylo/static` — see
`ngphylo-static`/`ngphylo-standalone-static` in `docker-compose.yml`/
`docker-compose.standalone.yml`), the files `init` collects vanish when
it exits, and `web` starts with nothing there. Under `DEBUG=True`
(`settings.local`, the default), this is invisible - Django's dev server
serves static files directly from each app's own `static/` source
directory and never touches `STATIC_ROOT` at all. It only 404s everything
under `STATIC_URL` once `DEBUG=False` (`settings.prod`) actually takes
effect and whitenoise starts serving from `STATIC_ROOT` instead - which is
exactly how this went unnoticed for so long, compounding with the
`DJANGO_SETTINGS_MODULE` bug above (a deployment that never actually ran
`settings.prod` never hit this either).

### Workflow duplicates and the Celery cleanup jobs

Every real OneClick/A La Carte run gives its own Galaxy-side copy of the
workflow it launches (`Workflow.duplicate()` in `workflows/models.py`) —
same name as the base workflow (e.g. `"FastME OneClick"`), fresh Galaxy id,
tracked locally as its own `category='duplicated'` row. Over time, real
usage means Galaxy's `/api/workflows/` list accumulates many entries
sharing the exact same name as the one true base workflow.
`importworkflows` matches by name (`re.search('oneclick', wfname, ...)`),
so it has to explicitly skip any Galaxy id it already knows about under any
category — otherwise it tries to collapse every same-named entry into the
single `category='base'` row via `slug=slugify(wfname)`, and crashes with
`duplicate key value violates unique constraint
"workflows_workflow_id_galaxy_key"` the moment one of those ids already
belongs to an existing `'duplicated'` row. Only surfaced by redeploying
against a real Galaxy that had accumulated enough real workflow runs — see
`ImportWorkflowsCommandTest.test_ignores_per_run_duplicates_sharing_the_same_name`
in `workflows/tests.py`.

The two Celery-beat cleanup jobs that delete old data
(`workspace.tasks.deleteoldgalaxyhistory`, `workflows.tasks.deleteoldgalaxyworkflows`
— `CELERY_BEAT_SCHEDULE` in `settings/base.py`, daily at 2am) call
`deletegalaxyworkflow()`/`deletegalaxyhistory()`, which each swallow their
own Galaxy API exceptions and just log a warning. Both cleanup tasks now
check the actual return value and only mark their Django rows deleted (or
hard-delete them) when the Galaxy-side delete really succeeded — a row that
fails is left `deleted=False` and retried on the next run. This used to be
unconditional: a transient Galaxy failure at exactly 2am made a row look
"cleaned up" in Django forever, since every future run only looks at
`deleted=False` rows, while the actual data could still be sitting on
Galaxy with nothing left to ever notice or retry it. See
`workspace/tests.py`/`workflows/tests.py` for the regression tests. Also
worth remembering: `Workflow.date`'s field default and both tasks' cutoffs
use `django.utils.timezone.now()`, not `datetime.now()` — `USE_TZ=True` is
on, and the naive version is exactly what throws `RuntimeWarning:
DateTimeField Workflow.date received a naive datetime while time zone
support is active` (harmless only by coincidence here, since
`TIME_ZONE='UTC'` matches the container's system clock).

`deleteoldgalaxyworkflows()` used to only delete a duplicated workflow if
it had **zero** associated `WorkspaceHistory` (i.e. created, e.g. the user
opened a submission form, but never actually run), on a 1-day cutoff -
anything that was actually run stayed in Galaxy forever unless its own
specific linked history happened to independently satisfy
`deleteoldgalaxyhistory`'s much narrower conditions (`finished=True`,
14-day cutoff, `workflow` FK actually set) at cleanup time. A real
production dump (`galaxy.pasteur.fr`, years of usage) turned up 832,000+
`category='duplicated'` rows never touched by either task - see
`scripts/cleanup_old_galaxy_workflows.sh` for the one-off bash/curl-based
cleanup this backlog needed (batches through Galaxy's own
`/api/workflows`, oldest-first, since the Django-tracked rows alone don't
reflect what's actually accumulated in Galaxy over that many years).
`deleteoldgalaxyworkflows()` now deletes any non-base workflow past a
7-day cutoff regardless of whether it was ever run - deleting a
workflow's Galaxy *definition* doesn't touch its history's actual data
(datasets/job outputs live in the History, not the Workflow), so there's
no need to wait for that history's own 14-day retention window. This
safely overlaps with `deleteoldgalaxyhistory`'s own per-history workflow
cleanup regardless of which of the two Celery tasks happens to run first
in a given 2am cycle: `w.delete()` here `SET_NULL`s any
`WorkspaceHistory.workflow` FK pointing to it, so `deleteoldgalaxyhistory`
correctly sees `workflow=None` and skips re-deleting anything already
gone if it runs second; if it runs first instead, it already marks the
row `deleted=True` itself, so this task's own `deleted=False` filter
skips it in turn.

### Daily workflow-usage report

`workspace.tasks.send_daily_report` (`CELERY_BEAT_SCHEDULE`, 8am UTC) emails
an HTML report — 7-day daily breakdown plus all-time totals, both by
category and, for OneClick, by which of the 4 tools — with matplotlib
charts embedded as base64 PNG `<img>` data URIs (not JS/SVG-based: email
clients generally don't execute JS, so this is the reliable approach). All
the data gathering/chart/HTML rendering lives in `workspace/reports.py`;
the task itself just calls `render_report_html()` and emails it — see that
module's docstrings for the full breakdown logic. No-ops (just logs) if
`NGPHYLO_REPORT_RECIPIENTS` isn't set, so it's safe to leave enabled
everywhere; see README.md's "Email" section for the env vars, and note
that plain job-completion emails need the same SMTP config and were never
actually wired into either `docker-compose.yml`'s or
`docker-compose.standalone.yml`'s env passthrough before this — check
`NGPHYLO_EMAIL_HOST` etc. are actually set on any deployment where email
(this report or otherwise) is expected to work.

One thing worth knowing if you touch the category breakdown:
`WorkspaceHistory.workflow_category` is **not** a clean one-to-one mapping
to "how the user submitted this" — see "Workflow duplicates and the Celery
cleanup jobs" above for why `'duplicated'` covers both the ordinary
Advanced-form path and an actual rerun. `reports.py`'s `CATEGORY_LABELS`
relabels it for display (`'duplicated'` → `"Advanced"`) but doesn't change
what's actually being counted.

### Kubernetes deployment (GitLab CI)

A third deployment path, alongside `docker-compose.yml` (local/dev) and
the IFB Cloud VM setup (`IFB_CLOUD.md`, untracked — see below): `.gitlab-ci.yml`'s
`build`/`deploy-dev`/`deploy-prod` stages build the same `Dockerfile` image,
push it to this project's own GitLab container registry, and `kubectl
apply` `manifest_datastores.yaml` (Postgres + Redis, each their own
Deployment+PVC/none) then `manifest.yaml` (a one-shot `Job` for
`docker/init.sh`, then `web`/`celery-worker`/`celery-beat` Deployments +
Service + Ingress) to a Kubernetes cluster. Structure mirrors
[drmab-web](https://github.com/evolbioinfo/drmab-web)'s own
`.gitlab-ci.yml`/`manifest.yaml`/`manifest_mysql.yaml` closely (same
`docker:dind` build job, same env-var-templated-via-`envsubst` deploy
jobs, same `kubectl patch ... labels: {date: ...}` forced-rollout trick),
adapted for this app's shape.

**Galaxy is external, not deployed by this** — per an explicit decision:
`NGPHYLO_GALAXY_URL`/`NGPHYLO_GALAXY_KEY` point at the Pasteur Galaxy
server, the same way drmab-web's own `GALAXYURL`/`GALAXYKEY` point at a
Galaxy instance it doesn't run either. Deploying Galaxy itself into
Kubernetes (its current job-execution model is privileged Docker-in-
Docker — see `NGPhylogeny_fr_galaxytools`'s CLAUDE.md — which doesn't
translate directly to Kubernetes without re-architecting around Galaxy's
own Kubernetes job runner or the official Galaxy Helm chart) was
deliberately out of scope here.

**`deploy-dev`'s cluster access is real** (since 2026-09): namespace
`ngphylogenyfr-dev`, GitLab Environment `k8sdev-ngphylogenyfr-dev`, domain
`ngphylogenyfr.dev.pasteur.cloud`, Deploy Token + kubectl context already
wired into a runner. **`deploy-prod`'s are still a placeholder guess**
(`ngphylogeny-prod`/`k8sprod-ngphylogeny`/`ngphylogeny.pasteur.cloud`,
following drmab-web's own naming convention) — nothing prod-side is
provisioned yet; confirm the real values once it is, and update both
`.gitlab-ci.yml` and (for the public hostname) manifest.yaml's
`NGPHYLO_HTTPS_HOST`/Ingress `host` to match, same as was just done for
dev. `deploy-dev` triggers on every push to `upgrade`; `deploy-prod`
requires a manual trigger from the pipeline page even then, on purpose —
nothing rolls out to production automatically.

**Static files are baked into the image at build time**
(`Dockerfile`'s `RUN python manage.py collectstatic --noinput`), not
shared via a volume between the init step and `web` the way
`docker-compose.yml`'s `ngphylo-static` volume does — Kubernetes' one-shot
init `Job` and the `web` Deployment are separate pods with no shared
filesystem by default, and provisioning a PVC just to share static assets
between them would be needless complexity when baking them into the image
works everywhere (docker-compose included — `docker/init.sh` still runs
its own `collectstatic` too, redundant but harmless there).

**Both Secrets (`ngphylogeny-credentials` in `manifest.yaml`,
`postgres-credentials` in `manifest_datastores.yaml`) are `stringData`
templated via `envsubst` from masked/protected GitLab CI/CD variables**
(`POSTGRES_PASSWORD`, `DJANGO_SECRET_KEY`, `ADMIN_PASSWORD`, `GALAXY_KEY`,
`EMAIL_HOST_PASSWORD` — see `.gitlab-ci.yml`'s `.deploy` comment for the
full list), not committed placeholder base64 values — unlike drmab-web's
`mysql-credentials`, which does check in a placeholder to decode/replace
later. Went straight to the more secure form here since this pipeline
already had a real cluster to target the moment it was written, rather
than leaving a live secret-rotation step for later. Same reasoning as
`NGPHYLO_GALAXY_KEY` going through a Secret at all in the first place
(unlike drmab-web's plain-`envsubst`'d `GALAXYKEY`) — no reason to leave
an API key or the DB/Django secrets less protected than they need to be.

**Every container (`init` Job, `web`, `celery-worker`, `celery-beat`) sets
`DJANGO_SETTINGS_MODULE` explicitly** to `NGPhylogeny_fr.settings.prod` in
its own `env:` list (no anchor/YAML-reuse across them — anchors don't
survive across `---`-separated documents in the same file, only within
one; `kubectl apply --dry-run` catching an unresolvable anchor reference is
exactly how this was found while writing these manifests). This is the
same variable-name trap documented in "`NGPHYLO_SETTINGS_MODULE` vs
`DJANGO_SETTINGS_MODULE`" above for the IFB Cloud deployment — getting it
wrong doesn't error, it just silently serves `DEBUG=True`.

### `upgrade` vs the old `master` branch

`upgrade` (this branch) is a from-scratch Python 2→3 / Django 1.11→4.2 /
Docker-Compose rewrite done independently of `master`, which stalled on the
old Python 2 stack. They share history up to a merge-base years back, so
`master` isn't something to merge in wholesale — but it still holds a
handful of genuine bugfixes made against the real production instance that
predate the rewrite and were worth cross-checking for. That comparison
turned up:

- **Raw (non-bioblend) Galaxy HTTP calls need the `x-api-key` header added
  by hand.** `bioblend>=1.4.0` moved auth from `?key=` query params to an
  `x-api-key` header automatically for anything routed through it — but
  `data/views.py`'s `download_file`/`tree_visualization`/`export_to_itol`
  build and fetch URLs directly via `urllib.request.urlopen()`, bypassing
  bioblend entirely, so they never picked that up. Whether this is a hard
  requirement depends on the target Galaxy's `require_login`/dataset
  permission settings: tested against the local dev Galaxy image
  (`require_login` unset → permissive), these endpoints return `200` even
  fully anonymous, so the header isn't *load-bearing* there — but sending
  an invalid key does get a hard `401` (Galaxy checks it if present, it
  just doesn't require one), and `gi.key` here is always a real stored key
  (never `None`/garbage — `GalaxyUser.get_galaxy_instance()` raises rather
  than hand back an instance with no key), so there's no regression risk in
  adding it. Fixed to keep this consistent with every other Galaxy call in
  the app and to not depend on the target Galaxy being configured
  permissively.
- **`Tool.import_tools()` re-fetched and re-saved every `Citation` row on
  every `importtools` run**, even for tools that already existed — no dedup
  on that table, so re-running the command (which happens on every
  container restart via `docker/init.sh`) duplicated citations
  indefinitely. Fixed by only fetching/saving citations when a tool is
  newly created or force-reimported — same idempotency class as the
  `importworkflows`/`addgalaxykey` fixes below.
- **Deliberately NOT ported: `master`'s `importworkflows` matching regex
  change** (`'oneclick'` → `'/oneclick'`, slash-prefixed) and its plain
  get-or-create dedup. `upgrade`'s own `importworkflows`/`addgalaxykey`
  already use `update_or_create` (galaxy re-imports its bundled workflows
  with a fresh `id_galaxy` on every restart — `master`'s `id_galaxy`-keyed
  existence check doesn't actually solve that duplication, `upgrade`'s
  `galaxy_server`+`slug`-keyed one does). Adopting the slash-prefixed regex
  would also be a breaking change against any already-deployed instance
  whose imported workflows are named `"<Tool> OneClick"` (space) rather
  than `"<Tool>/OneClick"` (slash) — check the actual imported workflow
  names on the target Galaxy before ever changing this regex.
- Everything else `master`-only (`startup.sh`/old `Dockerfile`/old
  `.dockerignore` changes, the `addgalaxykey.py` get-or-create dedup, a
  couple of dependabot bumps) was superseded by this branch's own
  independent fixes, or made moot by the docker-compose rewrite.
