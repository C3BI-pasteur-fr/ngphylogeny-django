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

**"BLAST" is a 5th category, merged in from a different app/table
entirely** — `blast.BlastRun` (see "App responsibilities") has no
`WorkspaceHistory` row and no `workflow_category`, so it's not a real
value of that field. `gather_last_7_days()`/`gather_all_time()` each
separately query `BlastRun` (grouped by `date`/counted overall) and merge
the result into the same `by_day_category`/`by_category` dicts under a
synthetic `'blast'` key, alongside `CATEGORY_LABELS['blast'] = 'BLAST'`/
`CATEGORY_COLORS['blast']`. Both BLAST servers (NCBI and Pasteur) are
lumped into this one category, same as how OneClick/Advanced/A La Carte
are each already a single category regardless of which specific tool
ran. Every category-consuming function downstream (the daily/all-time
charts, `category_columns`/`per_category`/`alltime_category_table` in
the template) is already fully generic over whatever keys are in
`CATEGORY_LABELS`, so this needed no template changes at all — only
`reports.py`'s two `gather_*` functions and the label/color dicts. Counts
`deleted=True` `BlastRun` rows too, same "usage report, not a
what's-still-retained report" reasoning as `WorkspaceHistory` (see this
section's first paragraph). Alongside this, `BlastRun.date`'s field
default (`blast/models.py`) and `deleteoldblastruns()`'s cutoff
(`blast/tasks.py`) were switched from naive `datetime.now()` to
`django.utils.timezone.now()` — same bug class as `Workflow.date` above,
just not yet fixed for the `blast` app specifically, and would have
made the new day-bucketed BLAST query's day boundaries unreliable.

**Restoring historical `workspace_workspacehistory` data** (e.g. into a
fresh deployment like `ngphylogenyfr-dev`, so the report reflects real
usage instead of just a handful of test submissions) only needs that one
table from a production dump — `reports.py` reads `created_date`,
`workflow_category`, `workflow_steps`, and `workflow__name` (the last via
a LEFT JOIN through the `workflow` FK, used only as a display fallback
and gracefully `None` if it doesn't resolve). It deliberately does **not**
need `workflows_workflow` restored alongside it (a real production dump
can have 800,000+ rows there from years of per-run duplicates — see
"Workflow duplicates" above — reimporting it would undo any cleanup done
on the target Galaxy for no report benefit) or `auth_user` (the report
doesn't group by user). Restore by extracting just that table's `COPY`
block from the dump (`grep -n "^COPY workspace_workspacehistory "
dump.sql` to find it, then the matching `\.` terminator), with three
columns rewritten before import:
- `workflow_id` / `user_id` → `\N` (both nullable - drops the FK
  references to data deliberately not being restored).
- `id` → dropped from the `COPY` column list entirely, letting Postgres
  assign fresh ids via the table's own sequence - avoids any collision
  with rows already created by real usage/testing on the target
  deployment (a production dump's ids can be arbitrarily large/small) and
  needs no manual sequence resync afterward.
- `galaxy_server_id` is **not** nullable, so it can't just be dropped -
  wrap the `COPY` in `BEGIN; ALTER TABLE workspace_workspacehistory
  DISABLE TRIGGER ALL; COPY ...; UPDATE workspace_workspacehistory SET
  galaxy_server_id = (SELECT id FROM galaxy_server WHERE url LIKE
  '%<real galaxy host>%') WHERE galaxy_server_id = <dump's original
  value>; ALTER TABLE workspace_workspacehistory ENABLE TRIGGER ALL;
  COMMIT;` - Django's FK constraints aren't `DEFERRABLE`, so temporarily
  disabling the table's FK-enforcement triggers is the only way to let
  the value briefly not resolve to a real row while the fixup `UPDATE`
  runs in the same transaction.

**Both restore scripts (`prepare_legacy_dump_import.sh`/
`prepare_ifbcloud_delta_import.sh`, both untracked) now also restore
`blast_blastrun`/`blast_blastsubject`, alongside
`workspace_workspacehistory` as before**, so the "BLAST" category (see
above) can reflect real pre-existing historical volume too, not just
usage since the restore. `blast_blastrun`/`blast_blastsubject` needed
a different id/FK strategy than `workspace_workspacehistory` - verified
directly against a real `manage.py migrate`-generated schema (a
throwaway Postgres via `docker run postgres:15`), not assumed from the
old dump's own (not necessarily current) DDL:
- `blast_blastrun.id` is a `uuid` with **no DB-level default** (Django
  assigns it client-side via `uuid.uuid4()` on save, not a Postgres
  sequence/`gen_random_uuid()`) and the row has no FK to anything - so
  its original id is restored completely as-is rather than
  dropped-and-reassigned, which is also collision-safe (uuid). One
  consequence: unlike `workspace_workspacehistory` (fresh
  sequence-assigned id every run), re-running `prepare_legacy_dump_import.sh`'s
  *output* a second time against the same target fails outright on a
  primary key violation (transaction rolled back cleanly, verified) -
  this is a one-shot historical import, not meant to be re-run, whereas
  the IFB Cloud script's delta approach naturally dedupes by the same
  preserved `id` directly (simpler than `workspace_workspacehistory`'s
  `history`-column indirection, needed there only because its own id
  *is* reassigned).
- `blast_blastsubject.id` *is* a real Postgres identity column with
  nothing else referencing it, so it's dropped/reassigned - same
  reasoning as `workspace_workspacehistory.id`. Its `blastrun_id` is
  kept as-is (referencing the preserved `blast_blastrun` ids above).
  Every restored `blast_blastrun` row has `deleted` forced to `true`,
  same reasoning as `workspace_workspacehistory`'s rows: so
  `blast.tasks.deleteoldblastruns()`'s daily cleanup (`deleted=false`
  filter) never tries to call Galaxy's real delete-history API against
  years-old, near-certainly-defunct Pasteur BLAST histories.
- The `blast_blastrun`<->`blast_blastsubject` FK is `DEFERRABLE INITIALLY
  DEFERRED` on the real deployed schema (also verified directly, not
  assumed) - unlike `workspace_workspacehistory`'s FKs, which aren't -
  so no trigger-disabling is needed for these two tables, just COPYing
  `blast_blastrun` before `blast_blastsubject` in the same transaction.
- The IFB Cloud delta script now takes a 4th argument, a separate
  tracking file for already-imported `blast_blastrun` ids (distinct
  from the existing history-id tracking file - different table, different
  natural key) - `touch` it once before the first run, same as the
  existing tracking file. A `blast_blastsubject` row is only ever
  imported alongside a `blast_blastrun` row that's new in that same
  run, so it needs no tracking file of its own.

Both extended scripts were verified end to end against a real Postgres
(not just read over) - real sample rows from `ngphylo_dump.sql`,
migrated with the actual `manage.py migrate` schema, imported, and
queried back; the delta script was also run through three rounds
(initial import, an unchanged re-run producing zero new rows, and a
real incremental delta) to confirm dedup and the tracking-file updates
behave correctly. One thing that surfaced during that verification, true
of the original script too and not something either script guards
against: the tracking file(s) are updated as soon as the SQL is
*written*, regardless of whether it was actually successfully applied
to the target Postgres - if a generated import fails partway (transaction
rolled back), its ids are still recorded as "imported," and the next
delta will skip them. Rerunning after a failed apply needs the affected
line(s) manually removed from the relevant tracking file first.

**`prepare_legacy_dump_import.sh` writes two separate output files, not
one** (`<dump> <workspacehistory_output.sql> <blast_output.sql>`) - each
wrapped in its own `BEGIN`/`COMMIT`, so either can be applied without the
other (e.g. restoring only the BLAST history into a deployment that
already has its workflow history, or vice versa). `prepare_ifbcloud_delta_import.sh`
still writes one combined output - not split, since it wasn't asked for
there, but the same approach would apply if it ever is.

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

**`deploy-dev`'s cluster access is real and confirmed working end to end**
(since 2026-09): namespace `ngphylogenyfr-dev`, GitLab Environment
`k8sdev-ngphylogenyfr-dev`, domain `ngphylogenyfr.dev.pasteur.cloud`,
Deploy Token + kubectl context wired into a runner — a real OneClick
workflow has been submitted and completed successfully against it. **See
"Getting `deploy-dev` from green pipeline to actually working" below for
everything that took to get there** — none of it was a manifest/pipeline
authoring mistake caught by review, all of it only surfaced by actually
running a real deploy against a real cluster. **`deploy-prod`'s
namespace/environment/domain are still a placeholder guess**
(`ngphylogeny-prod`/`k8sprod-ngphylogeny`/`ngphylogeny.pasteur.cloud`,
following drmab-web's own naming convention) — nothing prod-side is
provisioned yet; confirm the real values once it is (same process as
dev: update `.gitlab-ci.yml`'s `deploy-prod` block and, for the public
hostname, manifest.yaml's `NGPHYLO_HTTPS_HOST`/Ingress `host`), and
expect to hit the same class of first-real-deploy issues listed below
again against whatever cluster/namespace prod actually turns out to be.
`deploy-dev` triggers on every push to `upgrade`; `deploy-prod` requires
a manual trigger from the pipeline page even then, on purpose — nothing
rolls out to production automatically.

**Getting `deploy-dev` from green pipeline to actually working** took
several rounds of real-cluster-only issues, worth knowing about before
repeating this for prod:
- **RBAC**: the GitLab runner's ServiceAccount
  (`system:serviceaccount:gitlab-runner:default`) has no write access to
  a namespace by default — `kubectl delete/create/apply` all fail with
  `Forbidden`. `k8s-rbac-gitlab-runner.yaml` (checked into this repo, not
  applied by the pipeline itself — the runner's SA can't grant itself
  more access) is the `Role`/`RoleBinding` a cluster admin needs to apply
  once per namespace, granting CRUD on Secrets/Deployments/Jobs/Services/
  Ingresses/PVCs/Pods.
- **Runner selection**: this project has multiple runners (some `k8s`-
  tagged, presumably others not), but the `.deploy` job currently sets no
  `tags:` at all — it relies on whichever runner(s) accept untagged jobs
  actually being the right one for the cluster this repo targets. If a
  redeploy ever lands on the wrong runner, this looks identical to a
  missing RBAC grant (same generic `gitlab-runner:default` identity in
  the resulting `Forbidden` error either way) — check the job's actual
  runner on its GitLab page, and this project's registered runners under
  Settings → CI/CD → Runners, before assuming it's RBAC again.
- **The global `before_script` breaks the deploy job**: `pip install -r
  requirement.txt` runs by default for every job, but `.deploy`'s image
  is kubectl/yum-based with no Python at all — needs its own empty
  `before_script: []` override (same fix `build` already needed for its
  own reasons).
- **Migrations vs. a *persistent* database**: `makemigrations`-fresh-
  every-start (see "Migrations are not committed" above) only produces a
  correct schema the *first* time it runs against a given database.
  `ngphylogenyfr-dev`'s Postgres is a PVC-backed Deployment, not wiped
  between deploys — its `django_migrations` table already has
  `<app>.0001_initial` recorded as applied from the very first `init`
  run. On every later redeploy, `makemigrations` regenerates a migration
  file with the *same* name (Django numbers from scratch since there's no
  history file to build on) but whatever the *current* `models.py` says,
  and `migrate` matches purely by app+name — sees that name already
  applied, and silently no-ops the real `ALTER TABLE`, even though the
  regenerated file's actual operations changed. Any field change (e.g.
  `Citation.reference` going from `CharField(1000)` to `TextField` this
  session) needs a manual `ALTER TABLE ... ALTER COLUMN ... TYPE ...`
  run directly against the live Postgres pod to actually take effect —
  `migrate`'s own success/failure tells you nothing about whether it did.
  This will recur for every future model change deployed here; the
  durable fix (committing migration files to git, so `migrate` can apply
  real incremental changes) was discussed and deliberately deferred in
  favor of handling drift manually case by case.

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

**The init Job's `GALAXY_WORKFLOW_IDS` CI/CD variable (plain, optional)
should be set to the 4 known galaxy.pasteur.fr base-workflow ids**:
`0c0a83400cbba3e9` (FastME/OneClick), `7e182aa7ef0fb860`
(FastTree/OneClick), `617de6dd70aae83a` (PhyML/OneClick),
`6f4b7c17419da3e5` (PhyML+SMS/OneClick) — long-lived ids, created once on
2019-01-16 and never touched since (found via the real production
database dump, `ngphylo_dump.sql`, not the k8s namespace's own near-empty
dev DB). Passed as `NGPHYLO_GALAXY_WORKFLOW_IDS` to `docker/init.sh`,
which uses `importworkflows --wfids=...` (fetches those specific ids
directly, one GET per id) instead of `--wfnamefile=wfnames.txt` (lists
and filters Galaxy's whole workflow collection by name) whenever this
variable is set — see `tools/management/commands/importworkflows.py`.
Left unset, `docker/init.sh` falls back to the name-based listing
behavior unchanged, which is what `docker-compose.yml`/
`docker-compose.standalone.yml` still use (a different, freshly-created
Galaxy instance on each of those has no such stable pre-known ids).

**Maintenance mode**: the `web` Deployment's (only - not `init`/
`celery-worker`/`celery-beat`) `NGPHYLO_MAINTENANCE_MODE` env var, driven
by the `MAINTENANCE` GitLab CI/CD variable ("true" to enable, anything
else/unset to leave it off), makes `NGPhylogeny_fr.middleware.
MaintenanceModeMiddleware` serve `templates/maintenance.html` (503) for
every request instead of routing normally - a previously-dead, hardcoded-
stale-date template that used to be wired in via a commented-out catch-
all URL pattern (`re_path(r'.*', ...)`, now removed) rather than a
setting. The message is deliberately generic ("currently under
maintenance... we'll be back shortly") rather than naming a specific
reason/date, since it's meant to be reusable every time this gets flipped
on, not rewritten per-incident. Takes effect on the next deploy, same as
every other CI/CD variable here - re-run/retry the deploy job to flip it
either way, changing the GitLab variable alone doesn't touch an
already-running pod.

**Pasteur BLAST activation**: `settings.BLASTS['pasteur']['activated']`
(gates the Pasteur Galaxy BLAST server option throughout `blast/models.py`
- 6 separate checks) is env-var driven the same way, via
`NGPHYLO_PASTEUR_BLAST_ENABLED`/the `PASTEUR_BLAST_ENABLED` CI/CD
variable, set on both `web` and `celery-worker` (`launch_pasteur_blast`/
`checkblastruns` in `blast/tasks.py` read it too, not just the form).
Off by default, matching the hardcoded `False` it replaced. **Mutually
exclusive with `BLASTS['ncbi']['activated']`, not independently
toggleable** - enabling Pasteur deactivates NCBI's public server option
(previously unconditionally `True`), and vice versa: prefer Pasteur's own
controlled Galaxy BLAST server exclusively once it's available, rather
than also still offering NCBI's shared, rate-limited public
infrastructure. Both read the same `_PASTEUR_BLAST_ENABLED` module-level
variable in `settings/base.py` (Python dict literals can't
cross-reference each other's values directly).

**`launch_pasteur_blast()` had its own Python 2->3
`NamedTemporaryFile()` binary-mode bug** (`blast/tasks.py`), only
surfaced once Pasteur BLAST activation above actually let a real
submission reach it: `tmp_file.write(sequence)` wrote the (plain `str`)
query sequence into a `NamedTemporaryFile()`, which defaults to binary
mode - `TypeError: a bytes-like object is required, not 'str'`. Same bug
class as `workflows.tests.ProcessFileToUploadTest`'s
`process_file_to_upload()` fix. Fixed with
`tempfile.NamedTemporaryFile(mode='w')`; regression test
`blast.tests.LaunchPasteurBlastTest`.

**`checkblastruns()` had a real race with `launch_pasteur_blast()`,
only surfaced once a real Pasteur submission actually raced its own
1-minute Celery-beat check**: `launch_pasteur_blast()` saves the run as
`PENDING` right after creating its Galaxy history, then only sets
`history_fileid` afterwards, once the (network-bound) file upload + tool
run calls finish. If `checkblastruns()` polls in that window, it calls
`galaxycon.histories.show_dataset(b.history, '')` - the empty dataset id
turns the URL into Galaxy's history *contents list* endpoint instead of a
single dataset, returning a `list`, not a `dict`: `AttributeError: 'list'
object has no attribute 'get'`. Worse, the whole per-run loop used to
share one `try/except`, so this (or any other single run's failure)
silently aborted checking of every other pending/running Pasteur run in
that same pass too. Fixed by excluding `history_fileid=''` from the
polled queryset (that run is picked up again once it's set) and giving
each run its own `try/except` so one failure can't starve the rest;
regression tests `blast.tests.CheckBlastRunsTest`.

**`checkblastruns()` had no timeout at all on the actual Galaxy-side
blast computation**, only surfaced by a real Pasteur run against `nt`
that sat showing `Running` for over an hour with `message` empty and no
way to tell a genuinely slow search apart from one stuck on Galaxy's/the
cluster's side. This is a different phase than `launch_ncbi_blast`'s/
`launch_pasteur_blast`'s own `soft_time_limit`/`time_limit` (both only
bound *submitting* the job - the actual computation runs async on Galaxy
afterwards, polled here). `PASTEUR_RUN_STALE_AFTER` (`blast/tasks.py`,
currently 3 hours, same "generous enough for a real search, bounded
enough to recover" guess as the submission timeouts) is checked against
`BlastRun.date` at the top of each run's own per-run `try` (see the
`checkblastruns()` race fix above for why each run already has one) -
past the cutoff, `show_dataset` is skipped entirely (no point asking
Galaxy about a run already being abandoned), the run is marked `ERROR`
with a clear message, and its Galaxy history is queued for deletion
(`deletegalaxyhistory.delay(...)`, not called directly - same "don't
block this pass on one more Galaxy call" reasoning as the submission
timeout's own cleanup). Regression test:
`blast.tests.CheckBlastRunsTest.test_gives_up_on_runs_stuck_past_the_staleness_cutoff`.

**`deleteoldblastruns()` (the daily 2am, 14-day-cutoff cleanup) now
queues `deletegalaxyhistory` instead of calling it directly** - same
"don't block this batch on one Galaxy call" reasoning as the two
timeouts above, previously not applied here even though this is the
oldest of the three cleanup paths. It also **clears `query_seq`/`tree`**
on every run it cleans up, to free space on rows old enough to be
deleted anyway - confirmed safe for the daily report
(`workspace/reports.py` only ever reads `BlastRun`'s `date`/`deleted`/
`id`/`query_length`, never `query_seq`/`tree`). The redundant `e.save()`
right after `e.soft_delete()` (which already calls `self.save()`
internally) was also dropped. Regression tests:
`blast.tests.DeleteOldBlastRunsTest`.

**New field `BlastRun.query_length`** (`PositiveIntegerField(null=True,
blank=True)`) exists specifically so the sequence length survives
`deleteoldblastruns()` clearing `query_seq` - set at submission time
(`launch_ncbi_blast`/`launch_pasteur_blast`, right alongside `query_seq`,
before anything else - including the alphabet check - can short-circuit
the run into `ERROR`) and re-derived defensively at cleanup time for any
row where it's still `NULL` (predates the field, or some other path
never set it), right before `query_seq` is cleared. A row already
cleaned up *before* this field existed has no way to recover its length
- `query_seq` is already gone by then - and stays `NULL` permanently;
nothing to be done about that historical gap. `null=True` since existing
rows aren't backfilled (migrations aren't data-migrated in this
project). Needs a manual `ALTER TABLE` on any already-deployed Postgres,
same "migrate alone won't apply this" reasoning as `history`/
`history_fileid` above - verified against the real migration's own
generated SQL (`manage.py sqlmigrate blast 0003`), not guessed (Django's
`PositiveIntegerField` also emits a `CHECK` constraint):
```sql
ALTER TABLE blast_blastrun ADD COLUMN query_length integer NULL CHECK (query_length >= 0);
UPDATE blast_blastrun SET query_length = LENGTH(query_seq)
  WHERE query_length IS NULL AND query_seq IS NOT NULL AND query_seq <> '';
```
The `UPDATE` only recovers length for rows whose `query_seq` hasn't
already been cleared by a previous cleanup run - same unavoidable gap as
above for anything already cleaned up. Regression tests:
`blast.tests.LaunchNcbiBlastTest`/`LaunchPasteurBlastTest`/
`DeleteOldBlastRunsTest`.

**`BlastRun.history`/`history_fileid` were `CharField(max_length=20)` -
too narrow for this Galaxy server's real encoded ids.** Once the two
bugs above were fixed, a real Pasteur submission ran for real on Galaxy
(job genuinely executing) but crashed saving that fact back:
`django.db.utils.DataError: value too long for type character
varying(20)` on `history_fileid` in `launch_pasteur_blast()`'s final
`b.save()` - the run was left showing `PENDING` in NGPhylogeny
indefinitely while actually running/finishing on Galaxy. Bumped both to
`max_length=250`, matching `workflows.Workflow.id_galaxy`'s existing
convention for the same kind of value (an opaque Galaxy-provided encoded
id) elsewhere in this codebase. `workspace.WorkspaceHistory.history` has
the exact same `max_length=20` shape and hasn't hit this yet - only
because OneClick/Advanced history ids have happened to fit so far - and
would need the same fix if it ever doesn't.

Can't be caught by `manage.py test`: CI's test DB is sqlite (no
`NGPHYLO_DATABASE_HOST` set - see `.gitlab-ci.yml`'s `test` job), which
doesn't enforce `CharField` `max_length` at the DB layer the way
Postgres does; `blast.tests.LaunchPasteurBlastTest`'s regression test
checks the field's `max_length` directly instead of reproducing the
crash. **And per "Migrations vs. a persistent database" above, this
needs the same manual fixup on any already-deployed Postgres**:
`makemigrations`+`migrate` regenerating a same-named `0001_initial` a
second time silently no-ops the real `ALTER TABLE` on a database that
already has that migration recorded as applied - run this by hand
against the live Postgres pod for any deployment that already has BLAST
data:
```sql
ALTER TABLE blast_blastrun ALTER COLUMN history TYPE varchar(250);
ALTER TABLE blast_blastrun ALTER COLUMN history_fileid TYPE varchar(250);
```

**BLAST completion email now reuses the workflow job-completion email's
branded HTML template**, rather than the hand-built plain-text
`send_mail()` call it used before. `workspace/emails.py`'s MIME/logo
wiring was split out into `build_branded_html_email(subject, plain_text,
html, recipient)`, and its `site_url()` helper (https if
`NGPHYLO_HTTPS_HOST` is set, http otherwise) was made non-private
(`site_url`, not `_site_url`) since it's now used across apps - both
shared by the new `blast/emails.py`
(`build_blast_completion_email()`/`send_blast_completion_email()`), which
`blast/tasks.py`'s three notification sites (`launch_ncbi_blast`,
`launch_pasteur_blast`'s counterpart in `checkblastruns`) now call
instead of building `send_mail()` messages by hand. Both notification
types render the *same* `templates/workspace/job_completion_email.html`
file directly (not a copy) - `success_message`/`error_message` context
variables are the only notification-specific text, so the two stay
visually identical rather than drifting apart over time. See
`blast.tests.BlastCompletionEmailTest`, which mirrors
`workspace.tests.JobCompletionEmailTest`'s coverage of the same shared
machinery.

**`launch_pasteur_blast()` got the same 10-minute
`soft_time_limit`/`time_limit` as `launch_ncbi_blast`**, plus cleanup of
any Galaxy history it already created before timing out. Unlike NCBI,
the actual blast computation runs asynchronously on Galaxy once
submitted and is separately monitored (still with no timeout of its
own) by `checkblastruns()` - a hang inside `launch_pasteur_blast` itself
is most likely in the (network-bound) `create_history`/`upload_file`/
`run_tool` calls that submit the job in the first place. On
`SoftTimeLimitExceeded`, if `b.history` was already set (the Galaxy
history got created before the timeout fired), `deletegalaxyhistory` is
queued (`.delay()`, not called directly - this task already blew its
own time budget) to clean it up rather than leaving an orphaned history
nothing will ever reference again; previously an orphaned history like
this would just sit until `deleteoldblastruns()`'s 14-day cutoff.
Regression test: `blast.tests.LaunchPasteurBlastTest
.test_timeout_marks_error_and_cleans_up_galaxy_history`.

**Deleting a BLAST run (`DeleteBlastRunView`) crashed with
`TemplateDoesNotExist: blast/blastrun_confirm_delete.html`** - a real
Django 4.x breaking change, not something introduced by this project's
own rewrite. `BaseDeleteView.post()` was rewritten to go through
`FormMixin` (`get_form()`/`form_valid()`/`form_invalid()`) and no longer
calls `self.delete()` at all (the class docstring's own
`DeleteViewCustomDeleteWarning`, visible in the worker/web log right
before the crash, says exactly this). The view's `get()` used to forward
to `self.post()` to skip `DeleteView`'s confirmation page (this view
never had one) - but `get_form_kwargs()` only binds
`request.POST`/`request.FILES` onto the form when `self.request.method`
is actually `'POST'`, which it still isn't when `post()` is called by
hand from inside a GET request. The resulting unbound form is always
invalid, so it fell through to `form_invalid()`'s default: render the
confirmation template - which was never created because this view never
wanted one. **A real POST request had a second, more severe, silent bug
for the same reason**: `BaseDeleteView.form_valid()` calls
`self.object.delete()` directly - Django's real hard delete - so the
custom `delete()` override (meant to soft-delete via
`BlastRun.soft_delete()`) was already dead code on the POST path too,
just without a visible crash. Fixed by dropping the `delete()` override
entirely and handling both `get()`/`post()` directly with the same
soft-delete-and-redirect logic, bypassing `FormMixin`'s form machinery
altogether - this view never validated an actual form. Regression tests:
`blast.tests.DeleteBlastRunViewTest` (both confirmed to fail - one by
crashing, one by hard-deleting the row - against the pre-fix code).

**BLAST analysis was briefly, temporarily disabled** (code-level, not
via a CI/CD variable) right after real usage surfaced two open issues:
`launch_ncbi_blast`'s NCBI client could hang indefinitely with no
timeout (`blast/tasks.py`, now bounded to 10 minutes via
`soft_time_limit`/`time_limit`), and the Pasteur BLAST server option had
never actually been activated at all (now env-var driven, see above).
`BlastView.dispatch()` (`blast/views.py`) served
`templates/blast/blast_disabled.html` (503) instead of the real form,
and the "Blast Analysis" menu link was commented out in
`templates/base.html`. Both are reverted now that the underlying issues
are fixed - `templates/blast/blast_disabled.html` is left in the repo,
unreferenced, in case a future incident needs the same quick disable
again (`git log` for `blast/views.py`/`templates/base.html` has the
exact diff to reapply).

**`/status` is deliberately exempt from maintenance mode** - both
`readinessProbe` and `livenessProbe` on the `web` Deployment point at it,
not `/`. First real use of `MAINTENANCE=true` (2026-09-15, `deploy-dev`)
took the site *fully* down rather than showing the maintenance page: the
middleware's 503 on every path included whatever the probes hit, and
Kubernetes reads a failing `httpGet` probe as "the container is broken",
not "intentionally in maintenance" - repeated liveness failures
restarted the pod endlessly, and the Service ended up with zero ready
endpoints (the client-visible symptom was an infrastructure-level "no
available server", not even an app-level error). Immediate recovery was
`kubectl set env deployment/ngphylogeny-web -n <namespace>
NGPHYLO_MAINTENANCE_MODE=False` (patches the running Deployment directly,
faster than waiting on a pipeline run) plus flipping the `MAINTENANCE`
CI/CD variable back so the next deploy didn't reintroduce it. If any
future page/endpoint also needs to stay reachable during maintenance,
add it to `MaintenanceModeMiddleware.EXEMPT_PATHS`, not just to whatever
the probes happen to check.

**Every container (`init` Job, `web`, `celery-worker`, `celery-beat`) sets
`DJANGO_SETTINGS_MODULE` explicitly** to `NGPhylogeny_fr.settings.prod` in
its own `env:` list (no anchor/YAML-reuse across them — anchors don't
survive across `---`-separated documents in the same file, only within
one; `kubectl apply --dry-run` catching an unresolvable anchor reference is
exactly how this was found while writing these manifests). This is the
same variable-name trap documented in "`NGPHYLO_SETTINGS_MODULE` vs
`DJANGO_SETTINGS_MODULE`" above for the IFB Cloud deployment — getting it
wrong doesn't error, it just silently serves `DEBUG=True`.

### Two long-standing history/upload bugs only a real submission surfaced

Getting the very first real OneClick workflow to actually complete
against `deploy-dev` (see above) surfaced two bugs that predate this
session's Python 2→3 rewrite entirely — neither one is Kubernetes-
specific, both apply to every deployment, they were just never exercised
by a genuinely fresh session hitting these exact code paths before:

- **`workspace/views.py`'s `get_or_create_history()` returned the wrong
  type on a session's first history.** `create_history()` returns the
  full `WorkspaceHistory` model instance (other callers, e.g.
  `tools/views.py`, need it for its FKs/`wf_category`/`wf_steps`) — a
  return-type change from commit `488070129` (2018!) that
  `get_or_create_history()` never picked up, despite its own docstring
  saying `:return: history_id`. `data/views.py`'s `UploadMixin` (which
  `WorkflowOneClickListView` inherits `form_valid` from unmodified, so
  this is OneClick's actual paste/upload step) passes that straight into
  bioblend's `paste_content`/`upload_file` as `history_id`, which then
  crashed trying to JSON-encode it: `TypeError: Object of type
  WorkspaceHistory is not JSON serializable`. Only fires on a session's
  *first* history (later requests reuse the correct string id already
  cached in `request.session['last_history']`), which is exactly what a
  brand-new namespace's very first real submission guarantees and normal
  repeated testing in the same browser session almost never does. Fixed
  by unwrapping `.history` in `get_or_create_history()`.
- **`data/views.py`'s `UploadView.form_valid()` called `super().form_valid()`
  with no arguments** — Django's `FormMixin.form_valid(self, form)` has
  always required `form`; this specific call has been missing it since
  commit `e97ad0bc` (2018-03-14), predating even the bug above. The
  Galaxy upload itself (`self.upload_content()`/`self.upload_file()`)
  already succeeds by the time this runs, so the practical effect was
  always a 500 error page immediately *after* a technically-successful
  submission rather than a clean redirect to the history page — easy to
  mistake for a harmless glitch since the user's actual data/analysis
  still went through, which is the likely reason this went unnoticed for
  8 years. Fixed by passing `form` through.

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
