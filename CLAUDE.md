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

**24h staleness cutoff for regular workflow/tool runs**
(`workspace.tasks.WORKFLOW_RUN_STALE_AFTER`) - the same class of gap
`blast.tasks.checkblastruns()`'s own `PASTEUR_RUN_STALE_AFTER` already
covers for BLAST, just never ported over to `updateworkspacestatus()`: a
`WorkspaceHistory` used to stay `finished=False` forever if its Galaxy
jobs got stuck (a hung cluster node, a tool that never returns) -
`launchmonitorworkspaces` just keeps polling it every minute, and
`deleteoldgalaxyhistory` only ever looks at `finished=True` rows (see
above), so nothing ever cleaned it up either. Reference point is
`WorkspaceHistory.created_date`, same choice as `BlastRun.date` for
BLAST.

Once a history's still running/queued past that cutoff,
`cancel_stale_jobs()` calls bioblend's `gi.jobs.get_jobs(history_id=...)`
(one call, not filtered per-state server-side, to keep this endpoint's
own Galaxy load down - see the step-chain section below for why that
matters here) and `cancel_job()` on whatever it finds still
new/queued/running. `cancel_job()` is explicitly a *per-job* operation -
it deletes that job's own not-yet-materialized output dataset, but
every other, already-completed dataset in the same history (earlier
steps that finished fine) is untouched and stays downloadable - a
deliberate choice (confirmed with the user before implementing): a
stuck step shouldn't cost the user results they already have. The
newly-error'd datasets are picked up by `deleteoldgalaxyhistory`'s
normal 14-day cleanup from there, same as any other finished-with-error
history - no separate cleanup path needed.

The still-running/queued dataset *states* are forced to `'error'`
directly in the fetched `history_content` (not inferred by re-fetching
to see how Galaxy itself now represents a cancelled job's dataset) right
before the existing finished/error-detection loop runs - that loop then
naturally computes `finished=True`/`error=True` and fires the existing
job-completion email path unchanged, with no duplicated logic.

Test coverage: `workspace.tests.UpdateWorkspaceStatusStaleRunTest` -
goes through `galaxy_connection()` for real (this is a Celery task, not
a view - `request.galaxy`/`connection_galaxy` don't apply), same
`Server` + anonymous `GalaxyUser` DB fixture pattern established
elsewhere this session, with `HistoryClient.show_history`/
`JobsClient.get_jobs`/`JobsClient.cancel_job` patched directly. Confirms
only the actually-stuck job gets cancelled (not an already-`ok` one in
the same history), the already-`ok` dataset's own state is left alone,
and a history that hasn't hit the cutoff yet is left running untouched.

**`workspace.views.running_jobs_view` (`/workspace/running`, URL name
`running_jobs`) is a live, admin/staff-only page listing everything
currently running** - both `WorkspaceHistory` rows (`monitored=True`,
`finished=False`, `deleted=False`) and `BlastRun` rows (`status` still
`PENDING`/`RUNNING`, `deleted=False`), merged into one list and sorted
oldest-first, so a run approaching or past one of the two staleness
cutoffs above (`WORKFLOW_RUN_STALE_AFTER`/`PASTEUR_RUN_STALE_AFTER`)
surfaces at the top instead of being buried under newer ones. Per row:
name (linking to the history/BLAST detail page), type (`workflow_category`
mapped through `workspace.reports.CATEGORY_LABELS` - the same mapping
the daily report already uses - or, for BLAST, the server+program+status),
created date, runtime (just `now - created_date`, rendered via Python's
own `timedelta.__str__`, no extra dependency), and steps done/running -
counted from `history_content` for workflows (a dataset's own `state`),
or a fixed 0/1 for BLAST (it's not a multi-step pipeline the same way).
Deliberately **not cached** (unlike `daily_report_view`'s 15-minute
cache) - the whole point of this page is a live, up-to-the-moment
picture, not a stale snapshot of what was running a quarter-hour ago.
`history_content_json` is defensively filtered to dict-shaped entries
only, same reasoning as `WorkspaceHistoryObjectMixin.get_context_data`'s
own note on this (a genuinely malformed one shouldn't crash a page that
reads it).

The dataset-table look (card container, color-coded status pills,
line-art SVG icon buttons, fixed action-slot alignment - see the
step-chain section below) moved from being inline in
`history_contents_provenance_ajax.html`'s own `{% block stylesheet %}`
into `assets/css/custom.css` (global, loaded on every page via
`base.html`) once this page needed the exact same look for an unrelated
table, rather than duplicating that CSS a second time - `@keyframes
workflow-step-pulse` moved with it, since `.status-pill-running`'s own
animation depends on it regardless of which page it's used from.

Test coverage: `workspace.tests.RunningJobsViewTest` - same
access-control pattern as `DailyReportViewTest` (anonymous/non-staff
redirected to admin login; a staff user sees the page), plus content
checks that a finished or deleted `WorkspaceHistory`/`BlastRun` never
shows up, that mixed workflow+BLAST rows sort oldest-first together,
and that the done/total step count is correct for a partially-completed
run.

### Graphical step-chain on the history detail page

`templates/workspace/include/history_contents_provenance_ajax.html`
(included by `workspace/history.html`, the `HistoryDetailView` template)
shows a horizontal chain of colored boxes above the existing dataset
table - one box per tool step, green/blue/gray/red for
finished/running/pending/error. Feasible cheaply because the data was
already there: `object.history_content` (fetched fresh on every page
load by `WorkspaceHistoryObjectMixin.get_object()`/
`updateworkspacestatus`) already carries each Galaxy dataset's real-time
`state`, and `WorkflowGalaxyFactory` chains tools by matching EDAM
formats step-to-step, so these workflows are strictly linear (no
branches) - a straight left-to-right chain is an accurate
representation, not a simplification. Real Galaxy workflow invocations
too (`gi.workflows.invoke_workflow()`, `workflows/views/generic.py`/
`wkadvanced.py`), not tool-by-tool orchestration - which is what makes
"pending" steps visible at all: Galaxy creates every step's output
dataset (state `new`) up front on invocation, not just the one
currently running, so a not-yet-reached step already exists in
`history_content` before it starts.

Built as an entirely separate JS pass from the table's own tool-name
resolution (same `get_dataset_tool`/`get_tool_name` AJAX endpoints,
called again independently rather than sharing the table's DOM/AJAX
calls) - deliberately not refactored to share one round of calls,
so this addition can't regress the existing (working, somewhat fragile)
table-grouping JS. One real difference from the table: the table sorts
newest-first (`dictsortreversed:"hid"`) to show recent activity at the
top; the chain sorts oldest-first (`dictsort:"hid"`) since a
left-to-right sequence of boxes reads naturally as progress through the
pipeline. Colors reuse `templates/workspace/job_completion_email.html`'s
existing `.status-ok`/`.status-error` palette for consistency.

No test goes through `HistoryDetailView` itself - `galaxy.decorator.
connection_galaxy` (the decorator gating that view) has no mocking
pattern anywhere in this codebase to build a real request/response test
on, and building one from scratch was out of scope for what's
fundamentally a template/JS change. `workspace.tests.
HistoryStepChainTemplateTest` instead renders `workspace/history.html`
directly via `render_to_string()` with a fake object - same approach
`DailyReportTest` already uses to test `render_report_html()` without
going through its own view.

**Partial refresh instead of a full `location.reload()` every 10s.** The
step chain + table (everything that actually changes as a run
progresses) now lives in its own template,
`templates/workspace/include/history_contents_refreshable.html`, `{%
include %}`d once into `history_contents_provenance_ajax.html` inside
`<div id="history-refreshable-region">`. `workspace.views.
HistoryContentRefreshView` (`GET /workspace/history/<id>/refresh`, urls.py
name `history_content_refresh`) renders just that fragment; client-side,
`refresh()` calls jQuery's `$('#history-refreshable-region').load(url)`
on the same 10s timer that used to call `location.reload()`. `.load()`
specifically executes `<script>` tags in the fetched HTML (unlike a
plain `$.get()` swapped in via `.html()`, which doesn't reliably run
embedded scripts) - so the fragment's own scripts (step-chain building,
table tool-name grouping, tooltip init, the "stop polling once finished"
check) re-run on every poll exactly as they do on first load, with no
separate client-side re-init path to maintain. `window.historyRefreshTimer`
is a global (not a local `var`) specifically so the refreshable
fragment's own finished-check can `clearInterval()` it from inside
`.load()`-injected script, on every reload, not just the first.

Splitting out the fragment means anything whose *filling* script only
ever runs once at the outer page's `$(document).ready` - table state,
the step chain - has to actually re-run on every poll, not just render
once and sit there stale; everything that needs that moved into the
fragment. The `#pub-container` citations list is the opposite split:
the `<ul id="pub-container">` element itself stays in the outer shell
(it isn't duplicated between the live/staging containers the way
`#workflow-step-chain`/`#myTable` are, so nothing about the staging
swap needs it to move), but its `get_dataset_citations` AJAX *call* was
moved into the fragment anyway - it used to fire once, from the outer
shell's own `$(document).ready`, which meant a page loaded during the
"please wait" state (before any tool had actually run yet) never
fetched citations again once the run produced some, even after
switching to the real step-chain/table view. Now it's part of the
fragment's own per-poll script (a plain, unscoped `$('#pub-container')`
- safe regardless of staging, precisely because that element is never
duplicated), so it picks up new citations on the very next poll instead
of never.

**Background-build the refreshed fragment instead of reloading it
straight into the visible region.** A plain `.load()` into
`#history-refreshable-region` on every poll (the version above) still
had a visible flaw: `.load()` blanks the target immediately, and the
fragment's own tool-name resolution is several sequential/parallel AJAX
round trips (`get_dataset_tool` per dataset, `get_tool_name` per tool
group, for both the step chain and the table independently - see
above) that can visibly take a moment, so the whole region flashed empty
and slowly rebuilt on every single tick. Fixed by adding a second,
permanently-hidden container, `#history-refreshable-staging`
(`display:none`, sits right next to `#history-refreshable-region` in the
outer shell), and having `refresh()` `.load()` into *that* instead, with
a `?staging=1` query param (`HistoryContentRefreshView.get_context_data()`
reads it into a `staging` context flag). The fragment template computes
`var $root = ...#history-refreshable-staging or #history-refreshable-region...`
right from that flag, and scopes every single selector inside it through
`$root.find(...)` instead of a bare `$('#id')` - necessary because, for
however long a staged build takes, `#workflow-step-chain`/`#myTable`/a
dataset's own id exist twice in the live document at once (the old,
still-visible copy, and the new one quietly being built) - `$root.find
('#id')` stays correct even then: jQuery's `getElementById`-based fast
path for id selectors only kicks in when the search context is the
`document` itself, so scoping to a plain element context (`$root` here)
falls through to a real, properly-scoped `querySelectorAll` search
instead (checked directly in Sizzle's own source, not assumed). Once
*both* the step chain and the table report fully resolved - via two new
`$.Deferred()`s, `window.__historyStepChainReady`/`__historyTableReady`,
only now resolved after the per-group/per-dataset tool-*name* calls
finish too, not just the earlier tool-*id* wave the old code already
waited on - a final `{% if staging %}`-only script moves the staged
content over wholesale (`$('#history-refreshable-region').empty()
.append($root.contents())`) and clears `window.historyRefreshInFlight`.
That flag (plus a 30s watchdog `setTimeout` as a fallback, in case an
AJAX failure ever left a deferred permanently unresolved) is what
`refresh()` checks before starting the next tick, so an overlapping
build can never start mid-swap. The very first render (direct `{%
include %}`, no query param, `staging` unset/`False`) skips all of this
by construction: `$root` resolves straight to `#history-refreshable-region`
and the swap script isn't even in the rendered output, so there's no risk
of it ever running against an untouched, empty `#history-refreshable-staging`
and wiping the just-rendered live page.

**The "please wait, analysis initializing" state was unified into the
same fragment/shell instead of staying a second, separate template.**
`templates/workspace/include/history_wait.html` used to be shown by
`workspace/history.html` in place of `history_contents_provenance_ajax.html`
whenever a history had fewer than 2 datasets yet (a run just starting) -
its own copy of the name/email inputs (differently styled, plus a Url
field the real view no longer has at all - see the panel-restyle work
above), a full `location.reload()` every 10s (none of the staged-
background-refresh work above applied to it), and its own copy of the
"This page is will be refreshed in N sec." banner. All of that drifting
out of sync with the real view was the actual bug being fixed here, not
just the individual symptoms - so rather than patch `history_wait.html`
to separately re-implement the same panel styling/staged-refresh a
second time, it's deleted outright and that state is now just a branch
inside `history_contents_refreshable.html` itself: `{% if
object.history_content and object.history_content|dictsortreversed:
"hid"|length > 1 %}` picks between the step-chain/table markup (unchanged)
and a plain "please wait" message + spinner, with `history.html` now
unconditionally including `history_contents_provenance_ajax.html`
(previously also `{% if %}`-gated on the same dataset-count check).
Both `__historyStepChainReady`/`__historyTableReady` resolve immediately
in the wait branch (nothing async happens there) so a staging swap never
waits on them. The info-refresh banner itself - in both the old
`history_wait.html` and the real view - is gone entirely, not just moved:
it described the old full-page-reload behavior and was never updated for
the background-staged refresh, which doesn't visibly reload anything for
a banner to announce.

**This feature's 10s-polling frequency against Galaxy caused a real
production incident: the pod itself got killed and restarted (observed:
5 restarts, exit code 137/SIGKILL), not just one endpoint erroring.**
First surfaced as a `bioblend.ConnectionError` (a `galaxy.pasteur.fr`
502, "Proxy Error... Error reading from remote server") logged as an
unhandled Django 500 from `workspace.views.get_dataset_toolprovenance` -
but an unhandled exception in one Django view can't crash/restart a pod
by itself (Django's exception middleware catches it, logs a 500, keeps
serving other requests) - so that traceback was a symptom of the same
underlying Galaxy outage, not the actual crash mechanism. The real
mechanism: `GalaxyUser.get_galaxy_instance` (`galaxy/models.py`)
constructed its bioblend `GalaxyInstance` with bioblend's own default,
`timeout=None` (wait forever) - and the k8s deployment's `web` container
runs uwsgi with only `--processes 4 --threads 2` (manifest.yaml, 8
request-handling slots total) plus a 120s `--harakiri`. The step chain/
table polls Galaxy-backed AJAX endpoints (`get_dataset_toolprovenance`
once per dataset, `get_tool_name` once per tool group, from *both* the
step chain and table scripts independently - see above) every 10s, for
every open history page - at that call volume, once `galaxy.pasteur.fr`
started degrading, enough of those 8 slots ended up stuck waiting
(indefinitely, since nothing bounded them) on a slow/unresponsive Galaxy
that new requests - *including the liveness probe's own plain GET
/status* - couldn't get a free slot in time. `failureThreshold: 3` failed
probes later (manifest.yaml, `periodSeconds: 10`, default 1s probe
timeout), Kubernetes decided the container was unhealthy and killed it -
which just repeated on restart, since the same polling resumed against a
Galaxy that hadn't actually recovered yet.

Fixed with two changes together, since neither alone is sufficient:
- **`GalaxyUser.GALAXY_REQUEST_TIMEOUT = 30`**, set as a plain attribute
  on the constructed `GalaxyInstance` after the fact (bioblend's own
  `GalaxyInstance.__init__`, unlike the lower-level `GalaxyClient` it
  wraps, doesn't actually accept/forward a `timeout` kwarg at all in the
  installed version - verified directly, not assumed - so it has to be
  set as `gi.timeout = ...` post-construction; `make_get_request`/
  `make_post_request` just read `self.timeout` off the instance on every
  call, so this takes effect immediately). Bounds how long any single
  Galaxy call can occupy one of the 8 worker slots, well under uwsgi's
  own 120s harakiri, so a stalled/degraded Galaxy frees workers back up
  quickly instead of the pool draining to zero.
- **The three endpoints' exception handling was broadened from just
  `except ConnectionError` to `except (ConnectionError,
  requests.exceptions.RequestException)`.** Verified directly (not
  assumed) that this distinction matters: bioblend's own retry loop
  (`bioblend.galaxy.client.Client._get`) only catches
  `requests.exceptions.ConnectionError` itself (converting it to
  bioblend's own `ConnectionError` after `max_get_attempts` - 1 by
  default, matching the observed traceback's "0 attempts left") - a
  plain `requests.exceptions.ReadTimeout` (the likely shape of a slow/
  overloaded-but-still-connected Galaxy, and exactly what the new
  `timeout=30` above is meant to turn a stuck call into) is a
  `requests.exceptions.Timeout` subclass, *not* a `ConnectionError`
  subclass, and would otherwise propagate straight through both
  bioblend's retry logic and a narrower `except ConnectionError` alike.
  All three (`get_dataset_toolprovenance`/`get_dataset_citations`,
  `get_tool_name`) return a clean non-2xx JSON response on either
  exception type instead of an unhandled Django 500 - not itself what
  prevented the pod restarts (see above), but still worth fixing on its
  own merits: `get_dataset_citations` loops calling
  `show_dataset_provenance` once per dataset in the *whole* history
  within a single request, so it's even more exposed than the other two,
  and now `continue`s past just the one bad dataset instead of failing
  the whole citations fetch. The client side already degrades gracefully
  from a failed request without any changes needed: the step-chain/
  table's own `$.when.apply($, deferreds).then(success, fail)`/
  `.always()` handling (see above) already resolves
  `__historyStepChainReady`/`__historyTableReady` either way, so a poll
  that hits this just leaves that one box/cell with its placeholder text.

`galaxy.tests.GalaxyUserGetGalaxyInstanceTest` covers the timeout fix
directly (a plain model property, no `connection_galaxy` decorator
involved) - the existing `Server.save()`-makes-a-real-HTTP-call-on-create
wrinkle (see "Restoring historical..."/elsewhere) has an established
mocking pattern already (`patch('galaxy.models.requests.get', ...)`),
used here too. At the time this was written there was no regression
test for the broadened exception handling itself - reaching that
specific path means driving a *decorated*
(`galaxy.decorator.connection_galaxy`) view through Django's test
client with a mocked bioblend call - verified instead by inspection
(including live-checking bioblend's actual exception hierarchy/retry
behavior in a real Python shell, not guessed) and the existing suite/
flake8 staying green. That pattern turned out to already exist, though
(`tools.tests.GetToolNameViewTest`, predating this session's work) and
got reused for real by the load-reduction work below
(`workspace.tests.HistoryContentRefreshViewTest`) - a real
`Server`/anonymous `GalaxyUser` DB fixture, then `self.client.get(...)`
with the specific bioblend client method patched
(`bioblend.galaxy.histories.HistoryClient.show_dataset_provenance` etc.)
rather than mocking HTTP directly.

**Further cut the Galaxy call volume behind this feature - the actual
root cause of the incident above, not just its two symptoms (the 500s,
the pod restarts) - since step chain/table/citations were each
independently re-deriving the exact same dataset->tool_id/tool_id->name
mappings every 10s poll.** Before this: roughly `3N + 2K` Galaxy calls
per poll (N datasets, K distinct tools) - `get_dataset_toolprovenance`
once per dataset from *both* the step chain and the table, `get_tool_name`
once per tool group from *both* of those too, and `get_dataset_citations`
re-deriving the same per-dataset provenance a third time on top. Two
changes, both in `workspace/views.py`:
- **`resolve_dataset_tools(gi, galaxy_server, history_id, dataset_ids)`**
  - the one piece that still needs a real Galaxy call (bioblend has no
  bulk provenance-by-dataset-ids API) - now called *once* per poll, by
  `HistoryContentRefreshView.get_context_data()`, instead of by three
  separate call sites. Tool *names*, though, are resolved from this
  app's own `Tool` model (mirrors every tool NGPhylogeny actually runs -
  see "App responsibilities") before ever calling Galaxy's `show_tool` -
  since NGPhylogeny only ever submits its own preconfigured/imported
  workflows, the tool that produced any given dataset is almost always
  already known locally, cutting the `K` term to near zero in the
  common case (a DB lookup, not a Galaxy call, and one that doesn't
  depend on Galaxy being reachable at all). `tools.views.get_tool_name`
  itself got the same local-first lookup for consistency, even though
  nothing polls it anymore post-refactor (kept as a still-valid,
  standalone endpoint rather than deleted, in case anything else ever
  wants a one-off tool-name lookup).
- **`build_citations(dataset_tool_ids)`** - reuses the same resolved
  mapping instead of `get_dataset_citations` re-fetching provenance for
  every dataset a third time.

Both are computed once, server-side, and passed into
`history_contents_refreshable.html`'s context as plain dicts/a list
(`dataset_tool_ids`/`tool_names`/`citations`) - the step chain, the
table, and the citations list (`#pub-container`) all just read from
them directly now (embedded as JS objects/an array via `escapejs`, same
escaping discipline as `historySteps` already used), with **no client-
side AJAX calls left in this fragment at all**. This also simplified
the staging-swap mechanics: `__historyStepChainReady`/`__historyTableReady`
still exist (the staging-swap script - see above - still waits on them),
but now resolve synchronously since there's no more async tool-name
resolution to wait for. Net Galaxy load per poll: `N` calls (down from
`3N + 2K`) - and even those only run once at all, not from two
independent client scripts. `get_dataset_toolprovenance`/`get_tool_name`/
`get_dataset_citations` (the original standalone AJAX endpoints) are
kept, untouched otherwise, in case anything else ever calls them - just
no longer polled by this page.

Test coverage: `workspace.tests.ResolveDatasetToolsTest` unit-tests
`resolve_dataset_tools`/`build_citations` directly against a mocked
`gi` (local-tool-preferred, Galaxy fallback, a failed dataset left out
rather than raising); `workspace.tests.HistoryContentRefreshViewTest`
goes through the real decorated view end to end (see the test-pattern
note above) and asserts the rendered fragment contains zero occurrences
of `get_dataset_tool`/`get_tool_name`/`get_dataset_citations` - the
actual thing this change set out to prove, not just that the helper
function works in isolation.

**Real bug from the change above: `dataset_tool_ids`/`tool_names`/
`citations` were only ever computed by `HistoryContentRefreshView.
get_context_data()` - `HistoryDetailView` (the very first render, before
any poll) never set them at all**, so its own inclusion of
`history_contents_refreshable.html` rendered the Tool column and
citations list permanently unresolved. Harmless for a still-running
history (the first poll, 10s later, backfills it) but permanent for an
*already-finished* one: `object.finished` being true clears the
client's own polling timer before it ever fires once (see
`history_contents_provenance_ajax.html`'s `{% if object.finished %}`),
so `/refresh` is never called at all. Caught live on a real completed
`deploy-dev` run (`workspace/history/<id>`) - fetched the actual page
(publicly viewable, no auth needed for an anonymous-session history)
and read the embedded `datasetToolIds`/`toolNames` JS literals directly
out of the raw HTML (this only works because they're rendered server-
side as literal JS source, not fetched at runtime - the same property
that made this change possible in the first place) to confirm both were
empty `{}` despite 13 real, fully-`ok` datasets - not a hunch, verified
against the real deployed data. **Not** a dataset-collection issue
(checked first, since `show_history(contents=True)` mixes datasets and
collections with no code anywhere handling the latter's different shape
- ruled out because the affected row was confirmed to be a plain
dataset).

Fixed by moving the whole computation up into
`WorkspaceHistoryObjectMixin.get_context_data()` - both `HistoryDetailView`
and `HistoryContentRefreshView` inherit `(WorkspaceHistoryObjectMixin,
DetailView)` in that MRO order, so defining it once on the mixin (calling
`super().get_context_data()` to still reach `DetailView`'s own) reaches
both automatically; `HistoryContentRefreshView`'s own override shrinks to
just adding `staging` on top. Regression test:
`workspace.tests.HistoryContentRefreshViewTest.
test_history_detail_view_also_resolves_tool_names_on_first_load` - marks
the fixture history `finished=True` and hits `history_detail` (not
`history_content_refresh`) directly, exactly reproducing the real
failure mode.

**A second, compounding bug found while diagnosing the above: the
table's own grouping script collapses every row into one whenever
nothing resolves, not just leaving each with a blank Tool label.** Its
loop variable `id_tool` starts at `''`, and an unresolved dataset's own
cell also stays at its default `data=''` (nothing overwrites it - there
was nothing to set). A bare `id_tool == toolId` comparison then matches
`'' == ''` starting from the very first row, so every such row gets
folded into "the same tool" and `$(el).remove()`d except the first -
this is *worse* than just a blank label, it makes most of the table
vanish outright. This exact scenario (`dataset_tool_ids` completely
empty) is precisely what the `HistoryDetailView` bug above produced on
a real finished run, so both bugs compounded on the same real page.
Fixed with an explicit `toolId !== ''` guard alongside the equality
check, so an unresolved row is never treated as matching anything -
including another unresolved row - and keeps its own place in the table
with just a blank Tool cell, the intended degraded state. This is
defense in depth independent of the fix above: an individual dataset's
provenance call can still transiently fail on its own (see
`resolve_dataset_tools`'s own `continue`-on-failure), so a partially-
empty `dataset_tool_ids` is a real, expected case this has to handle
correctly regardless.

**Dataset table redesign**: `#myTable` (in `history_contents_refreshable.html`)
moved from a plain Bootstrap-striped table with bare glyphicon buttons to
a card container (`.history-table-card`), color-coded status pills
reusing the step-chain's own finished/running/pending/error palette
(`.status-pill-*`), and line-art SVG icon buttons (`.icon-btn`, inline
SVGs - not glyphicons, not emoji) instead of a row of identical small
buttons. `table-layout: fixed` with an explicit `<colgroup>` (Tool 15%,
Step 6%, File name 44%, Status 13%, Actions 22%) replaces the previous
`style="width: ..."` guesses on two of five `<th>`s. Iterated as an
Artifact mockup built from a real production history's actual rendered
data (fetched directly via `curl` - real tool names/groupings, not
invented ones) before touching the template, per the session's usual
"propose visually first" pattern for UI changes.

**Fixed action-slot alignment**: every row now shows the same action
icons in the same column position, whether or not that row actually has
each one - previously a variable-length list of only-the-applicable
buttons meant e.g. Download could be the 3rd icon in one row and the 2nd
in the next, depending on what else that row had. Canonical order:
session toggle, parameters, stdout/messages, download, display, a
*shared* tree-viewer/MSAViewer slot (never both on the same row, since a
dataset is never simultaneously a tree output and a fasta alignment),
then iTOL. A row missing a given action gets an invisible `.action-slot`
(same width as `.icon-btn`, or `.action-slot.wide` matching the 36px
`.itol-badge`) there instead of the next real button sliding left into
that column.

**Two real buttons the approved mockup didn't show, restored in the
actual implementation**: the mockup (being a quick visual proposal, not
a full behavior spec) dropped the "Show stdout" button entirely and
never depicted the `error`-state row's "Show messages" button at all -
implementing it directly from the mockup would have been a silent
functional regression, not just a restyle. Both are back: stdout uses a
plain chevron SVG in the shared 3rd slot; an `error` row uses that exact
same slot for "Show messages" instead (a row is never both `ok` and
`error`), styled `.icon-btn.danger` (red, reusing the status pill's own
error color) rather than `.warn` (amber - already used for iTOL, a
different, unrelated meaning) to keep the color coding consistent with
what "error" already means everywhere else on this page.

**The real app's own `static/images/ptree.svg` (used by the old "Viewer"
button) was deliberately not reused for the new tree-viewer icon.** It's
a solid white-filled mark meant for a colored button background: this
redesign's icon buttons are the opposite (light tinted background,
`currentColor`-stroked icon), so the white fill would simply be
invisible on them. Kept the inline cladogram SVG from the mockup instead
(a real branching-tree diagram, `stroke="currentColor"`, adapts to the
button's own color automatically) - already verified visually via the
Artifact preview before implementing, not a first attempt here.

Test coverage: `workspace.tests.HistoryTableRedesignTest` (an `ok` row
gets a `status-pill-ok`/"Done" pill and none of the old glyphicon
classes; a plain dataset with no special extension still gets its
trailing empty action slots so alignment holds; an `error` row gets a
`status-pill-error` pill and the `.icon-btn.danger` "Show messages"
button) - rendered directly via `render_to_string`, the same
no-view-mock pattern as this section's other template tests.

**The Actions column's initial percentage width (22%) clipped the last
icon(s) - iTOL first, being rightmost - on real deployments.** Caught
live (`curl`-fetched the actual page and confirmed the button really
was present in the server-rendered HTML - a CSS layout bug, not a
missing-button one) once real data was viewed at real page widths.
`table-layout: fixed` locks a percentage column to exactly that share of
the table's width regardless of content, and 7 icon slots (up to 36px
each, the `.itol-badge`) need roughly 270px+ - on any page width where
22% comes out narrower than that, the overflow got clipped by
`.history-table-card`'s own `overflow: hidden` (needed for the rounded
corners). Fixed by giving Actions (and Tool/Step/Status) fixed pixel
widths in the `<colgroup>` (140px/50px/-/90px/280px) instead of
percentages, and leaving File name's own `<col>` with no explicit width
at all - in `table-layout: fixed`, a column with no specified width
takes whatever's left after the fixed-width columns are subtracted, so
Actions is now guaranteed enough room regardless of table width, at the
cost of File name shrinking (it already truncates with an ellipsis + a
tooltip) on narrow pages instead.

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

**All-time BLAST query-length histogram** — a 6th chart, in the "Since
the beginning" card, unrelated to the category breakdown above (it's
not bucketed by day or category at all). `gather_blast_query_lengths()`
returns every non-`NULL` `BlastRun.query_length` (all-time, `deleted=True`
included, same reasoning throughout this section);
`render_blast_query_length_histogram()` renders it via `ax.hist(...)`,
capped at one bin per distinct value so a handful of searches doesn't
get smeared across 30 mostly-empty bins. Rows with `query_length` still
`NULL` are excluded from the histogram entirely (not counted as
zero-length) - see `BlastRun.query_length`'s own note above for why a
row can still be `NULL`. `CID_BLAST_LENGTH_HISTOGRAM`/
`blast_length_chart`/`blast_length_count` follow the exact same
cid-name/context-key/images-dict wiring as every other chart in this
module - the shared `CHART_CONTEXT_KEYS` swap (cid: for email,
data: URI for the web page) and `send_daily_report()`'s generic
`for cid, png_bytes in images.items()` attachment loop needed no changes
at all to pick this one up.

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
running a real deploy against a real cluster. **`deploy-prod`'s real
values**: namespace `ngphylogenyfr-prod`, GitLab Environment
`k8sprod-02-ngphylogenyfr`, public URL `ngphylogeny.fr` (external
ingress class - was briefly `ngphylogeny.pasteur.cloud`/`internal`
before the real production domain was confirmed and wired in), set via
`.gitlab-ci.yml`'s `deploy-prod` block's `PUBLIC_URL`/`INGRESS_CLASS`
variables (consumed by `manifest.yaml`'s `${PUBLIC_URL}`/
`${INGRESS_CLASS}` - the app's own `NGPHYLO_HTTPS_HOST` env var and the
Ingress `host`/`ingress.class` all derive from the same two variables,
not set independently). Expect to hit the same class of
first-real-deploy issues documented below the first time `deploy-prod`
is actually triggered against this domain, the same way `deploy-dev`
did. `deploy-dev` triggers on every push to `upgrade`; `deploy-prod`
requires a manual trigger from the pipeline page even then, on purpose
— nothing rolls out to production automatically.

**`deploy-prod`'s TLS is via cert-manager, patched in separately from
the shared `manifest.yaml`.** First attempt was a pre-existing
Pasteur-wide `wildcard-pasteur-tls` secret an operator suggested
(`*.pasteur.fr` wildcard) - didn't work, since `ngphylogeny.fr` isn't a
`pasteur.fr` subdomain at all (a separate domain registration, not
covered by that wildcard's SAN list). A manual `certbot`/DNS-01 attempt
was also started (independent of the cluster entirely - would have
needed a `_acme-challenge.ngphylogeny.fr` TXT record and a recurring
manual renewal, since nothing would auto-renew it) but abandoned once
the operator confirmed cert-manager is already installed cluster-wide
with a `ClusterIssuer` - the standard, auto-renewing path. **The
`ClusterIssuer`'s real name is `letsencrypt`**, not `letsencrypt-prod` -
that first guess (a plausible-looking name, not actually confirmed)
went in first and failed with `clusterissuer.cert-manager.io
"letsencrypt-prod" not found`, caught via `kubectl describe
certificaterequest` (the `Certificate`/`CertificateRequest` objects
were created fine - `ClusterIssuer` resolution happens one step later,
inside the `CertificateRequest`'s own `IssuerNotFound` condition, so no
`Order`/`Challenge` objects ever got created at all while this was
wrong). Confirming the real name needed the operator - `ClusterIssuer`
is cluster-scoped, and the deploying user's own kubectl access got
`Forbidden` (not `NotFound`) on both `list` and `get` for it, so it was
never independently checkable from this project's own access level.
Since `manifest.yaml` is shared with `deploy-dev` through plain
`envsubst` (no conditionals), the `cert-manager.io/cluster-issuer`
annotation and `tls:` block aren't in `manifest.yaml` itself - baking
them in unconditionally would've also applied to dev's Ingress (with an
empty `secretName` there, untested, and dev's `internal` ingress class
may already handle TLS some other way outside this Ingress object
entirely). Instead, the `.deploy` script's own `kubectl patch ingress`
step (right after the deployment-rollout patches) applies both, gated
on `INGRESS_TLS_SECRET_NAME` being set - only `deploy-prod`'s own
`variables:` block sets `INGRESS_CLUSTER_ISSUER`/
`INGRESS_TLS_SECRET_NAME`, so the patch is a no-op for dev. `secretName`
(`ngphylogeny-fr-tls`) doesn't need to exist beforehand - cert-manager's
ingress-shim watches for the annotation and creates/populates/renews
that Secret itself once the annotation lands.

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

**Pasteur's wrapper tool ids changed** (`settings.BLASTS['pasteur']['progs']`,
`toolshed.pasteur.fr/repos/fmareuil/ncbi_blast_plus/...` ->
`toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/.../2.14.1+galaxy2`),
independently of anything in this session - see git history for the
exact diff. This surfaced two real bugs, both about that literal `+` in
the new ids' version suffix, not the id change itself:
- **`blast/urls.py`'s `available_blasts_dbs`/`blast_example` routes
  404ed for every Pasteur program.** Both take the tool id as a path
  segment, inserted as-is (not `encodeURIComponent`-escaped) by
  `templates/blast/blast.html`'s JS - their `prog` regex (`[\w/\.]+`,
  `\w` = letters/digits/`_`) didn't allow `+`. Symptom: the BLAST
  database dropdown came back empty and the example-sequence button did
  nothing, for every Pasteur program (NCBI's own progs are plain words
  like `blastn` - never hit this). Fixed by widening the regex to
  `[\w/.+-]+` (`-` added too, pre-emptively - not currently in any id
  here, but a very plausible one in a toolshed owner/repo name).
  Regression tests: `blast.tests.BlastAjaxEndpointsTest`.
- **The `blastx` entry's id had a typo**: `toolshedtoolshed.g2.bx.psu.edu/...`
  (duplicated prefix) - Galaxy wouldn't recognize that tool id at all,
  breaking `blastx` submissions on Pasteur specifically. Not caught by
  any test (nothing exercises the actual tool id strings against a real
  Galaxy) - found by inspection when the other two bugs above were being
  diagnosed.
- Separately, `blast/tests.py` hardcoded the *old* blastn wrapper id as
  a literal string - broke this test file's own Pasteur-path tests the
  moment the real id changed, unrelated to either bug above. Now derived
  from `settings.BLASTS['pasteur']['progs']` by looking up the `'blastn'`
  type (`_pasteur_blastn_prog()`), so it can't go stale the same way
  again.

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

**`BlastRun.BLASTSERVERS` labeled NCBI runs as "Pasteur" too**
(`(NCBI, 'Pasteur')` instead of `(NCBI, 'NCBI')`) - found while building
`workspace.views.running_jobs_view` (its own `dict(BlastRun.BLASTSERVERS)`
lookup would have silently inherited this). No DB/migration impact -
`choices=` is display-only, the stored `server` value itself
('pasteur'/'ncbi') was never wrong. Fixing it also exposed a second,
previously-inert bug right next to it: `BlastRun.server_str()` compared
`self.status` (a `RUNSTATUS` code, e.g. `'P'`/`'R'`) against
`BLASTSERVERS` codes (`'pasteur'`/`'ncbi'`) - two code spaces that can
never match, so it always fell through to returning `'Error'`
regardless of the actual server. Nothing currently calls `server_str()`,
so this had zero live impact, but it was broken and sitting right next
to the fix - corrected alongside rather than left broken. Regression
tests: `blast.tests.BlastRunServerStrTest` (both server values return
their real label, not `'Error'`, on a plain unsaved instance);
`workspace.tests.RunningJobsViewTest` now also asserts an NCBI run
actually shows `"NCBI BLAST"` on the running-jobs page.

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

### Contact form now also emails staff, not just the DB row

`surveys.views.FeedbackCreateView.form_valid()` (previously the plain,
unmodified `CreateView.form_valid()` — no hook existed here at all)
now also sends a plain-text notification to
`settings.NGPHYLO_CONTACT_FORM_RECIPIENTS` (comma-separated, envsubst'd
from the `CONTACT_FORM_RECIPIENTS` GitLab CI/CD variable — see
`.gitlab-ci.yml`'s `.deploy` comment and `manifest.yaml`) *in addition to*
saving the `Feedback` row — the DB save was and remains the primary,
unconditional effect. Left unset (the default), `surveys/emails.py`'s
`send_feedback_notification_email()` no-ops before building anything,
same "safe to leave enabled everywhere" convention as
`NGPHYLO_REPORT_RECIPIENTS`/`send_daily_report` (see "Daily
workflow-usage report" above) — this setting is genuinely independent
of that one (a separate GitLab variable, a separate recipient list),
just following the same shape.

Deliberately plain text via a bare `EmailMessage`, not
`workspace/emails.py`'s branded `build_branded_html_email()` template
(unlike every other notification in this codebase — job completion,
BLAST completion): those exist to look good to an *end user* clicking
through to their results; this is an internal operational notice to
staff, carrying the submitter's own message verbatim, with `reply_to`
set to the submitter's email so a staff member can just hit reply
instead of copy-pasting their address out of the body — no logo/CID
MIME wiring needed for that.

Sending is wrapped in the exact same `try/except SMTPException` /
`except Exception` / `logging.warning` pattern already used by
`workspace.tasks.updateworkspacestatus` and `blast/tasks.py`'s own
completion-email sends — but for a different reason: those run inside a
Celery task, where a swallowed failure only affects that one async job.
This one runs inline in the submitter's own request/response cycle, so
without the same swallowing, an SMTP hiccup would turn an
already-successful submission (the `Feedback` row is saved *before* the
send is attempted, via `super().form_valid(form)`) into a 500 error page
for someone who just successfully contacted support.

`docker-compose.yml`/`docker-compose.standalone.yml` pass
`NGPHYLO_CONTACT_FORM_RECIPIENTS` through with the same
`${VAR:-}`-empty-default convention as every other optional email
setting; README.md's "Email" section and `manifest.yaml`'s three
container env blocks (`init`/`web`/`celery-worker`/`celery-beat` — added
to all four for consistency with `NGPHYLO_REPORT_RECIPIENTS`, even
though this send currently only ever happens on `web`, which is the one
that actually serves `FeedbackCreateView`) were updated to match.
Regression tests: `surveys.tests.FeedbackNotificationEmailTest` (email
content/recipients/no-op-when-unset, mirroring
`workspace.tests.JobCompletionEmailTest`'s style) and
`FeedbackCreateViewTest` (the view saves the row and calls the send
exactly once, and a raised `SMTPException` from the send doesn't
propagate — both driven through `form_valid()` directly with a mocked
form, not a real captcha-validated POST, since nothing in this codebase
has a pattern for bypassing django-simple-captcha in a test post; see
`data.tests.StaticPagesSmokeTest` for the existing real-render coverage
of that part of the form instead).

### Header Pasteur logo rendering full-size in Safari

`.pagehead-logo` (see "Can you add the white pasteur logo..." above) was
only constrained via CSS (`height: 44px; width: auto`) - no size hint
existed anywhere in the actual HTML. The source image
(`logo_institut_pasteur_white.png`) is 1793x661px; reported in Safari as
the logo rendering at roughly that native size, overwhelming the page -
Safari is known to be more visibly aggressive than Chrome/Firefox about
briefly (or, in some caching/load-ordering edge case) rendering an `img`
at its intrinsic size before/without the constraining stylesheet taking
effect, since nothing gave it an explicit box to reserve. Fixed by adding
real `width="119" height="44"` HTML attributes directly on the `<img>`
tag in `templates/base.html` (119 = the actual 1793:661 aspect ratio at
44px tall, computed via Pillow, not guessed) - HTML width/height
attributes are honored by every browser immediately, independent of
whether/when CSS has loaded, which is the standard fix for this class of
bug (also incidentally prevents layout shift). `.pagehead-logo` in
`assets/css/custom.css` additionally got `max-width: 160px; max-height:
44px;` as a belt-and-braces cap, in case anything else about the sizing
ever regresses again.

### Workflow list/form pages 500'd when Galaxy itself is unreachable

A real production incident: `galaxy.pasteur.fr` went into its own
maintenance mode (a plain 503 HTML page, not a Galaxy-shaped error) while
NGPhylogeny.fr was still up, and every workflow list/form page - OneClick
and Advanced both - 500'd instead of degrading gracefully. Real
traceback: `bioblend.ConnectionError: GET: error 503: ...Site
Maintenance...` raised straight out of `workflows/models.py`'s
`Workflow.fetch_details()` (`galaxyinstance.workflows.show_workflow(...)`),
via `workflows/views/generic.py`'s `WorkflowListView.workflow_list`
(a `cached_property` that calls `fetch_details()` once per base
workflow, completely uncaught) from a plain `GET /workflows/advanced/`.
Same uncaught call, same crash shape, in two more places reachable by
just browsing (not submitting anything): `WorkflowFormView.
get_context_data()` (the single-workflow OneClick page,
`workflow_oneclick_form`) and `WorkflowAdvancedFormView.
get_context_data()` (`workflows/views/wkadvanced.py`, the single-workflow
Advanced page, `workflows_advanced_fullsteps`). Since this app does no
computation of its own and only orchestrates Galaxy (see "What this is"),
a real Galaxy outage genuinely means these three pages can't be shown -
the fix isn't to route around Galaxy, just to fail predictably with a
clear message instead of Django's generic 500.

`workflows/views/generic.py` now defines `GALAXY_UNREACHABLE_EXCEPTIONS`
(`bioblend.galaxy.client.ConnectionError` + `requests.exceptions.
RequestException` - the same broadened pair already used throughout
`workspace/views.py`, see the "pod restart" section above for why
`ConnectionError` alone isn't enough - a `ReadTimeout` from
`GalaxyUser.GALAXY_REQUEST_TIMEOUT` is a `requests.exceptions.Timeout`
subclass, not a `ConnectionError` subclass) and
`galaxy_unavailable_response(request)`, which logs a warning and renders
the new `templates/workflows/galaxy_unavailable.html` (extends the
existing `error.html` shell - same pattern as `404.html`/`500.html` - with
a 503 status) instead. `WorkflowListView.get()` and `WorkflowFormView.
get()` (both in `generic.py`, covering OneClick+Advanced's list pages and
OneClick's single-workflow page) and `WorkflowAdvancedFormView.get()`
(`wkadvanced.py`, imports the same two names) each just wrap their
existing `super().get(...)`/`self.get_context_data(...)` call in a
`try/except GALAXY_UNREACHABLE_EXCEPTIONS` and delegate to the shared
helper - no other behavior changed.

**Deliberately not extended to every other Galaxy call site in this
app** (the workflow-maker form, rerun, the A La Carte POST path, actual
job submission) - those are either POST-triggered (a user actively
submitting work already expects some risk of failure, a materially
different situation from just browsing a list) or weren't part of the
reported/reproduced failure, and blanket-wrapping every Galaxy call
across the app was out of scope for this fix. `workflows_alacarte_build`
specifically wasn't touched because its own `GET` (the actual "A La
Carte" page load) doesn't call Galaxy at all - only building/submitting
a custom workflow on `POST` does.

Regression tests: `workflows.tests.GalaxyUnavailableTest` - real
`Server`/anonymous `GalaxyUser` DB fixture (the established
`connection_galaxy`-decorated-view test pattern, see the pod-restart
section above), `bioblend.galaxy.workflows.WorkflowClient.show_workflow`
patched to raise a real `ConnectionError` shaped like the actual incident
("Site Maintenance", `status_code=503`), asserting all four affected URLs
(`workflow_oneclick_list`, `workflows_advanced`, `workflow_oneclick_form`,
`workflows_advanced_fullsteps`) return `503` and render
`workflows/galaxy_unavailable.html` rather than raising.

### Permalink to "Workspace" (the previous-analyses list)

`request.session['histories']` (a plain list of Galaxy history ids -
`workspace.views.create_history`/`PreviousHistoryListView`) is the
*only* thing "Workspace" (`/workspace/histories`,
`previous_analyses.html`) is built from - clear cookies, switch
browsers/devices, or just come back much later, and the page shows
nothing, even though the underlying `WorkspaceHistory` rows (and their
Galaxy data, until the 14-day cleanup) are still there. Added a
permalink: `PreviousHistoryListView.get_context_data()` now also builds
`permalink_url` - `site_url(reverse('workspace_permalink', kwargs=
{'token': signing.dumps(history_ids, salt=PERMALINK_SALT)}))` - shown
on the page as a read-only, click-to-select input with a "Copy" button
(plain `document.execCommand('copy')`, matching this app's existing
jQuery-2.1.4/Bootstrap-3 stack rather than reaching for the newer
`navigator.clipboard` API). Built from `self.request.session['histories']`
*after* `get_queryset()` has already run (it prunes deleted/no-longer-
matching ids first - `BaseListView.get()` always calls `get_queryset()`
before `get_context_data()`), so the token never encodes a stale id
`get_queryset()` would just drop again anyway.

`workspace.views.WorkspacePermalinkView` (`GET
/workspace/permalink/<token>`, no `connection_galaxy` decorator needed -
it only touches the session, no Galaxy call) decodes the token with
`signing.loads(..., salt=PERMALINK_SALT)`, **merges** (set union, not
replace) the decoded ids into whatever's already in *that* browser's own
`request.session['histories']`, and redirects to `previous_analyses` -
so following a permalink in a browser that already has its own,
different analyses adds to that list instead of clobbering it. A
tampered/invalid token (`signing.BadSignature`) shows a plain error
message (`django.contrib.messages`, same `{% if "error" in message.tags
%}` block already established in `workflows/workflows_alacarte.html`,
now duplicated into `previous_analyses.html` too) and redirects back
rather than raising - it's user-shareable input, so it has to degrade
cleanly.

**Deliberately no encryption, no expiry, no server-side revocation, and
no new model/migration** - `django.core.signing` only guarantees the
token wasn't *tampered with*, not that it's secret, and that's fine
here: an individual history is already reachable the exact same way with
no ownership check at all -
`WorkspaceHistoryObjectMixin.get_object()` (`history_detail`,
`/workspace/history/<id>`) looks a history up by its Galaxy id with no
session/ownership gate whatsoever, so anyone who already knows (or
receives) one of these 16-character ids can open it directly. A
permalink bundling several of those same already-not-secret ids
together is not a new kind of exposure, just a convenient way to
restore or share the same list - so no expiry was added either (unlike,
say, a password-reset token): the whole point is that it should still
work long after a normal browser session would've expired, and a
permalink pointing at histories Galaxy has already cleaned up just
renders an empty/partial list once followed (`get_queryset()`'s own
`deleted=False` filter), no separate cleanup needed for the permalink
itself. Encoding the id list directly into a signed token (rather than a
new DB-backed model keyed by a random id) also means zero schema
changes - consistent with this project's general preference to avoid
new migrations when an existing mechanism (here, `SECRET_KEY`-backed
signing, already a Django built-in) does the job.

Regression tests: `workspace.tests.WorkspacePermalinkTest` - the
established `Server` + anonymous `GalaxyUser` DB fixture (`previous_analyses`
is `connection_galaxy`-decorated even though this feature itself never
calls Galaxy), covering: no permalink shown for an empty session; a
permalink is present and points at `/workspace/permalink/`; following it
with a **separate `django.test.Client`** (a genuinely different,
cookie-isolated session - the real scenario this exists for) restores
the right history list and renders that history's name; following it
when the current session already has a *different* history merges
rather than replaces; a garbage token redirects with the error message
instead of raising.

### Inline phylotree.js tree preview on the history detail page

As soon as a run produces a finished newick/nhx dataset - the same file
the dataset table's own `{% if file.extension in "nhx,nwk" %}` branch
shows the "Interactive Tree visualisation"/iTOL buttons for (see the
step-chain section above) - a small tree rendering now appears between
the workflow step chain and the dataset table, without waiting for a
full page reload. Reuses the exact same minimal `d3.layout.phylotree()`
pattern `templates/blast/result.html`'s own "Quick (and inaccurate)
tree" already established (a plain `<svg>`, no toolbar - unlike
`templates/treeviz/tree.html`'s full PRESTO-branded, toolbar-driven
viewer behind the "Interactive Tree visualisation" button, which stays
untouched): `phylotree.js-0.1.9` + `d3.v3.min.js`, both already vendored
under `assets/`, no new dependency added.

`WorkspaceHistoryObjectMixin.get_context_data()` (`workspace/views.py`)
picks the first dataset in `history_content` with `extension` in
`('nhx', 'nwk')` and `'ok'` in `state` - read straight off the
already-fetched `history_content`, no extra Galaxy call - and exposes
`tree_preview_dataset_id`/`tree_preview_url` (the latter a
`reverse('display_raw', ...)` URL for that one dataset's raw content).
Both are empty strings once there's no such dataset yet, same "always
present but empty" convention `dataset_tool_ids`/`tool_names` used
before this.

**The actual requirement - "should not be refreshed at every refresh
event" - takes real care given how this page already works.**
`history_contents_refreshable.html` (chain + table) is entirely
rebuilt/reloaded every 10s (see that section above) - a naive tree
preview living inside that same fragment would refetch the newick and
rebuild the whole d3 layout (losing any pan/zoom) on every single poll,
for a tree that never changes once it exists. Fixed by keeping the
*actual* rendered tree - `#tree-preview-wrapper` (holding `#tree-preview-
svg`) - entirely outside `#history-refreshable-region`/`#history-
refreshable-staging`, in the fixed outer shell
(`history_contents_provenance_ajax.html`, the one part of this page
never rebuilt by a poll): `window.maybeBuildTreePreview(datasetId, url)`
fetches+lays out the tree via `d3.layout.phylotree().svg(...)` **at most
once ever** per history (guarded by `window.treePreviewDatasetId`, only
set after a successful build - a failed fetch leaves it unset so a
later poll can still retry, rather than a transient Galaxy hiccup
permanently giving up on the preview for the rest of the session).

Positioning it correctly - "below the workflow preview, above the
table" - without rebuilding it is the other half of the problem, since
both the step chain and the table (and the space between them) live
inside the auto-rebuilt fragment. `history_contents_refreshable.html`
emits nothing but an empty `#tree-preview-anchor` placeholder div at
exactly that spot, regenerated fresh on every render like the rest of
that fragment; `window.placeTreePreviewIntoLiveRegion()` then
*reparents* the persistent `#tree-preview-wrapper` (a plain DOM move,
which - unlike innerHTML replacement - carries over its already-
rendered SVG and any d3-bound zoom/pan state untouched) into whatever
anchor is currently sitting in the **live** `#history-refreshable-region`
- never into the hidden staging container mid-build, which would
otherwise rip the visible tree out of the page for however long that
poll's staged build takes, then flash it back in on swap. Call sites,
matching this file's existing staged-vs-direct rendering split:
non-staged (the very first render, direct `{% include %}`) places it
immediately, since there's no staging step to race; a staged poll
defers the move to the exact same instant
`history_contents_refreshable.html`'s own swap script already moves the
rest of the freshly-built content into the live region
(`$('#history-refreshable-region').empty().append($root.contents())`) -
one atomic-looking update, not two. The fetch's own `.done()` callback
also calls the same placement function unconditionally (targeting the
live region directly, not `$root`), covering the case where the very
first fetch resolves after that poll's swap already happened.

Regression tests: `workspace.tests.TreePreviewTest` - reuses
`HistoryContentRefreshViewTest`'s own Server/GalaxyUser/local-Tool
fixture pattern, going through the real decorated
`history_content_refresh` view. Asserts `tree_preview_dataset_id`/
`tree_preview_url` are only wired into the rendered fragment once a
dataset is both `nhx`/`nwk` *and* `'ok'` (a still-`running` one doesn't
count), and that `#tree-preview-anchor` sits between
`#workflow-step-chain` and `#myTable` in the rendered HTML - the actual
ordering requirement, not just presence.

**Real bug on the first `deploy-dev` rollout of this feature: the tree
rendered genuinely empty.** Confirmed live (`curl`-fetched the actual
deployed page and its `treePreviewUrl` directly) that detection and
data fetching both worked correctly - `treePreviewDatasetId`/
`treePreviewUrl` were populated with the real dataset id, and that URL
returned a valid, real newick string - so the bug was specifically in
the client-side rendering, not the server-side wiring. Cause:
`#tree-preview-wrapper` started out `style="display:none"`, toggled via
jQuery's `.show()` once built. `maybeBuildTreePreview`'s very first
build (`d3.layout.phylotree()...layout()`) always runs while the
wrapper is still in that initial hidden state (placement into the
visible page only happens *after* building - see this section's own
description of the placement/build split above) - and phylotree.js's
`.layout()` measures the SVG via `getBBox()`/`getComputedTextLength()`
to size and place nodes/branch labels, both of which return zero (or
throw) for anything inside a `display:none` ancestor, since that
removes the element from the render tree entirely, not just from view
- unlike, say, `visibility:hidden` or an off-screen `position:absolute`,
which stay part of the render tree and remain genuinely measurable.
Fixed by hiding the wrapper with Bootstrap's own `.sr-only` class
instead (`position:absolute` + `clip`, never `display:none` - loaded on
every page already, via `base.html`'s Bootstrap 3 CDN link, so no new
dependency) and swapping `.show()` for `.removeClass('sr-only')`.
Regression test: `workspace.tests.TreePreviewTest.
test_wrapper_hidden_via_sr_only_not_display_none`.

**That `.sr-only` fix alone wasn't enough - the tree still rendered
empty on the next real report against the same URL, still with no
thrown exception anywhere.** Debugged directly against the live
deploy-dev page (not guessed): headless Chrome (`--headless=new
--dump-dom`, and separately driven live over the DevTools Protocol -
`websocket-client` in the `ngphylo-py38` env, `Runtime.evaluate`/
`Network.*` events) confirmed the newick fetch itself succeeded (a
direct `curl` to the embedded `treePreviewUrl` also returned a real,
well-formed newick string) and that `#tree-preview-wrapper` really did
get its `sr-only` class removed and get moved into the table's anchor
position - i.e. the fetch/placement machinery all worked - yet
`#tree-preview-svg` still had zero children: `.layout()` was still
being called *before* that removal/move, while the SVG sat inside
`.sr-only`'s 1x1px clip box. Oddly, replaying the exact same
fetch-then-layout sequence a few seconds later in that same live tab's
own console (once the initial page-load moment had passed) built the
tree correctly - and so did two separate offline reproductions of that
same "layout while still sr-only" sequence (one with a minimal
jQuery+d3+phylotree.js page, one loading this app's *entire* real
script/CSS list from the live deploy-dev static host) - neither could
reproduce the empty result outside that one specific real-page moment.
Rather than keep chasing what's likely a narrow, hard-to-isolate
browser layout-timing gap around phylotree.js's own measurement calls
during a real, busier page load, `renderPendingTreePreview()`
(`history_contents_provenance_ajax.html`) now (a) always reparents the
wrapper out of `.sr-only` and into the live anchor *before* ever calling
`.layout()` (previously still done after, despite the comment above
saying otherwise - that gap between comment and code was itself real
and is what let this second bug through), and (b) defensively verifies
`#tree-preview-svg` actually gained real children afterward, retrying
up to 5 times (plain `setTimeout`, not `requestAnimationFrame` - the
latter can be suspended indefinitely for a backgrounded/minimized tab,
which a retry meant to survive a real page's own load timing shouldn't
depend on) before giving up *for that poll* - `window.
treePreviewPendingNewick` is left set on a full give-up, so the next
10s poll's own `placeTreePreviewIntoLiveRegion()` call retries the
whole sequence again from a clean slate regardless. Existing tests
(`workspace.tests.TreePreviewTest`) still cover the server-rendered
side of this (dataset detection, ordering, `.sr-only` markup) unchanged
- the retry/verify logic itself is pure client-side JS with no Python
test harness in this codebase to drive it through a real browser (see
CLAUDE.md's own notes elsewhere on `connection_galaxy`-decorated views
having no such pattern); verified instead via the headless-Chrome
harnesses described above, run directly against this app's real
vendored static assets.

**Still not enough - the very next live report was "no panel at all",
not just an empty one.** Same debugging approach (headless Chrome,
`Runtime.evaluate`/`Network.*` over the DevTools Protocol against the
real deploy-dev URL): the fetch to `treePreviewUrl` had, on the real
page's own very first automatic attempt, resolved with a falsy/empty
body - confirmed live again by manually replaying the *exact same*
request a few seconds later in that same tab's console, which returned
the real 789-byte newick with no trouble at all. `maybeBuildTreePreview`
at the time keyed its "already handled this dataset" guard off having
merely *attempted* the fetch (`window.treePreviewDatasetId = datasetId`
was set unconditionally inside `.done()`, regardless of whether the body
was actually usable) - so that one bad attempt permanently blocked ever
retrying it again. On a *finished* history specifically (this exact
history's tree was already `'ok'`, meaning `object.finished` was almost
certainly true), `history_contents_refreshable.html`'s own `{% if
object.finished %}` block clears `window.historyRefreshTimer` on the
very first render - no poll ever fires again to retry from - so the
panel/heading themselves never even appeared, not just an empty tree.
Editing this comment block itself hit a real, separate gotcha: an
earlier draft literally wrote `{% if object.finished %}` inside a plain
JS `//` comment to explain this - Django's template lexer has no idea
it's inside a JS comment and tried to parse it as a real template tag,
raising `TemplateSyntaxError: Unclosed tag` and breaking the *entire*
page (caught by `manage.py test` before ever being pushed, not live) -
worth remembering for any future comment in this file that needs to
reference a `{% %}`/`{{ }}` construct by name: paraphrase it instead
("the object.finished conditional block"), never spell out the literal
tag syntax.

Fixed by decoupling both the "should we even try" guard and the retry
mechanism from page polling entirely: `maybeBuildTreePreview` now only
ever stops for good once `window.treePreviewRendered` is genuinely
true (set by `renderPendingTreePreview` only on real success - see
above), and on a falsy body or an outright AJAX failure it retries
*itself*, via its own `setTimeout` loop (a handful of attempts a couple
of seconds apart), entirely independently of `window.
historyRefreshTimer`/`refresh()` - so a finished history's one-shot
render still gets several independent chances to recover from a
one-off transient failure, not just the one. A later real poll (a
still-running history) can still short-circuit this early by
succeeding first, same as before. Verified end to end offline (headless
Chrome again, this app's real vendored static assets, `$.ajax`
temporarily replaced with a fake that returns an empty body twice then
the real newick on the third call, mirroring the actual observed
failure) - confirmed the tree renders correctly on that third,
self-triggered retry with no poll involved at all.

### Upgrading both tree-visualization surfaces to a modern phylotree.js

Both remaining consumers of the ancient, vendored `phylotree.js-0.1.9`
(2016-era bundled build) - the inline preview above and the full
"Interactive Tree visualisation" page (`templates/treeviz/tree.html`,
`data.views.tree_visualization`) - were moved to the actively-maintained
npm `phylotree` package (2.6.0, current at the time this was written).
It's the *same* upstream project (`veg/phylotree.js`), just years newer:
D3 v7 core, TypeScript, real tagged releases, none of the ad hoc bundled
concatenation the old vendored file had. `templates/blast/result.html`'s
own "Quick (and inaccurate) tree" was deliberately left on the old
0.1.9 version - out of scope for this change, and not worth touching
opportunistically.

**No build tooling exists in this project (no node, no webpack/esbuild),
so "upgrade the library" meant finding a genuinely CDN-servable browser
bundle, not just an npm package.** Confirmed directly (not assumed):
`phylotree`'s own `package.json` `unpkg`/`jsdelivr` fields point at
`dist/phylotree.js`, a real pre-built UMD bundle - vendored locally
under `assets/phylotree-2.6.0/` (`phylotree.js`, `phylotree.min.js`,
`phylotree.css`, and `phylotree-menus.css` - the last one is a separate
file the npm package ships under `src/render/styles/`, not bundled into
`dist/phylotree.css`, needed for the library's own built-in right-click
context menus to render at all).

**A real bug in that UMD bundle itself, found by actually reading the
unminified output, not assumed:** it declares both `underscore` and
`lodash` as external browser globals, but its own rollup build config
never gave them distinct global names - both defaulted to the
conventional `_`, and rollup only avoided a *local parameter name*
clash inside its own factory function by suffixing the second one to
`_$1`. That's purely an internal variable name; the bundle's actual UMD
header (verified in the unminified `dist/phylotree.js`) reads
`factory(global.phylotree = ..., global._, global._$1)` - i.e. it
genuinely expects a global literally named `_$1` for lodash, which
nothing ever defines on its own. Worked around with
`assets/phylotree-2.6.0/underscore-lodash-shim.js`: load underscore,
capture it (`window.__phylotreeUnderscore = window._`), load lodash
(which overwrites `window._`), then the shim swaps them into the two
slots the bundle actually reads
(`window._$1 = window._; window._ = window.__phylotreeUnderscore;`).
Verified this whole chain works in headless Chrome before writing it
into either template, both against a minimal page and against this
app's *entire* real script/CSS list served from a local instance - it
does.

**API differences from the old vendored version, also verified directly
rather than assumed from the (incomplete) upstream README/API docs:**
- Construction: `new window.phylotree.phylotree(newickString)` (the
  UMD global's own `phylotree` property is the class - confirmed by
  reading the bundle's own `exports.phylotree = Phylotree` line), not
  the old `d3.layout.phylotree()` factory function chain.
- Rendering targets a **container `<div>`**, not a d3 selection of an
  already-existing `<svg>` - `tree.render({container: "#id"|domNode,
  width, height, ...options})` builds the SVG internally. Critically,
  `render()` alone does **not** insert anything into the page - it
  returns a `TreeRender` object; the actual `<svg>` DOM node only comes
  from calling `.show()` on that object, which still has to be
  `appendChild`'d into the container by hand. This isn't documented
  anywhere in the README/API.md fetched during this work - found by
  reading `render()`'s own source and a helper call site
  (`select(rendered_tree.container).node().appendChild(rendered_tree.show())`)
  elsewhere in the bundle, then confirmed empirically (a
  `render()`-only call left the container empty; adding
  `.appendChild(display.show())` produced real content).
- Most option key *names* carried over unchanged from the old version
  (`align-tips`, `show-scale`, `selectable`, `collapsible`,
  `restricted-selectable`, `reroot`, `hide`, `brush`, `zoom`,
  `left-right-spacing`, etc.) - same project, same option vocabulary,
  confirmed by diffing the old vendored file's own internal defaults
  object against the new one's behavior, not by guessing from the name
  similarity alone.

**The inline preview** (`history_contents_provenance_ajax.html`)
follows the exact same build-once/retry/reparent architecture documented
above - only `renderPendingTreePreview`'s actual rendering call changed
(the old `d3.layout.phylotree().svg(d3.select(...))...layout()` become
`new phylotree.phylotree(newick).render({container: $container[0], ...})`
followed by `$container.append(display.show())`), and the retry-on-empty
check now looks for a real `<svg>` inside the container instead of
child nodes of the `<svg>` itself (since the container is a plain `<div>`
now - `#tree-preview-svg` in the markup, kept as-is despite no longer
literally being an `<svg>` tag, to avoid unrelated churn). Nothing about
upgrading the library explains *why* the earlier empty-container timing
bug happened in the first place, so the same defensive retry stayed.

**The full interactive viewer** (`templates/treeviz/tree.html`) is where
the real design decision was: rather than reverse-engineer every one of
the old PRESTO control panel's bespoke controls (phylogram/dendrogram
toggle, linear/radial/slanted layout, ladderizing sort buttons,
bootstrap/branch-length text toggles, tip alignment, manual zoom
buttons, a node-name search/highlight box - all hand-wired in the now-
deleted `assets/presto/load_tree.js` against the old library's API) one
by one against a different, only partially-documented modern API, the
whole bespoke panel was dropped in favor of the new library's own
built-in interaction system (right-click context menus, selection,
zoom/brush) - an explicit choice made with the user rather than
assumed. The page is now just a container `<div>` plus one `render()`
call with `selectable`/`collapsible`/`reroot`/`hide`/`brush`/`zoom`/
`show-scale` all enabled, and a one-line hint below the heading. The
"PRESTO" branding (a name for the old bespoke implementation
specifically) was dropped along with it, in favor of a plain
"Phylogenetic tree viewer" heading - keeping that name on a page that no
longer runs the actual PRESTO code would have been misleading. The
`assets/presto/` directory (d3 v3, `load_tree.js`, the old vendored
`phylotree.js`/`phylotree.css`, bundled font-awesome) is deleted
outright, not left around unreferenced - confirmed nothing else in the
codebase (`grep`'d templates and static references, not assumed) still
pointed at it once this page stopped using it.

**Verified directly, not assumed:** rendering (real SVG with the
correct node/branch counts, no exceptions) for both surfaces, in
headless Chrome, against this app's real vendored static assets served
from a local `docker compose -f docker-compose.standalone.yml -f
docker-compose.dev.yml` instance. This section originally also flagged
the built-in menu as "not verified" - a synthetic `contextmenu` event
never triggered it in headless Chrome - and guessed (wrongly) that this
was just a synthetic-event limitation, writing a "right-click a node"
hint into the page on that assumption. It's actually much simpler: the
menu is wired to a plain `"click"` handler in the library's own source
(`handle_node_click`, called from every node's circle/label click
listener - confirmed by reading it directly, not the incomplete docs),
never `"contextmenu"` at all - so the synthetic `contextmenu` dispatch
was testing an event the library never listens for in the first place,
and the real, live behavior (a real user right-clicking, expecting a
menu per that hint text, and getting the browser's own native menu
instead, since nothing here ever suppresses it) is exactly what got
reported live. Fixed by correcting the hint text to say "click", not
"right-click" - see the toolbar section below for where this surfaced
and how it was actually confirmed working (a real left-click, via
`elementFromPoint` + dispatched mouse events, does open
`.phylotree-context-menu` with real content).
No new Python-level tests were added for either template - both are
pure client-side rendering with no server-side logic change beyond
what `data.tests.TreeVisualizationTest`/`workspace.tests.TreePreviewTest`
already cover (context variables, markup positioning), which still
pass unchanged.

### A real toolbar was added back on top of the built-in UI, on both surfaces

The "lean on built-in UI" decision above still left both tree panels
with no *visible* controls at all - right-click-only discovery, plus
whatever the mouse/scroll gestures happen to be, isn't obviously
findable. A small toolbar was added back to both surfaces (built fresh
against phylotree.js 2.6's own real API, not a port of the old PRESTO
panel's code) rather than reintroducing the old bespoke controls
wholesale - the built-in interaction system (collapse, hide, reroot,
selection) is kept for everything not covered by an explicit button.

Both surfaces share the same shape: a small state object
(`treeOptions`/`window.treePreviewOptions`) plus a `showBranchLengths`
flag, and every control just flips one field and calls a shared rebuild
function - either a full `tree.render(treeOptions)` (layout/align-tips/
support-values - simplest and most reliable given phylotree.js 2.x
*does* expose live accessors for these, see below, but a fresh render is
still what the library's own internal reroot-handling code does
whenever anything changes, verified directly in the vendored bundle)
or, for the branch-length toggle specifically, `display.layout()` on
the *existing* render (keeps the current pan/zoom, since nothing about
that toggle needs a fresh tree). Buttons: a layout (linear/radial)
toggle, align-tips, show-support-values, show-branch-lengths, zoom in/
out, and reset-view (a full rebuild, which conveniently also clears the
zoom transform back to its initial state).

**Verified directly against the vendored bundle, not assumed from the
(incomplete) README/API docs, since the option-name-carried-over
assumption elsewhere in this section turned out right but this needed
checking too:**
- `radial(attr)`/`alignTips(attr)`/`internalNames(attr)` are real, live
  getter/setter methods on the `TreeRender` (`display`) object -
  `display.radial(true).layout()` (matching the *old* API's own
  `tree.radial(!tree.radial()).placenodes().update()` idiom almost
  exactly) re-lays-out in place without needing a fresh `render()` call
  at all. Used for the full-viewer toolbar's own toggles; the inline
  preview's toolbar still does a full rebuild for these three instead,
  for consistency with its own already-established
  `rebuildTreePreviewSvg()`/`treePreviewBranchLengthStyler` split rather
  than mixing two different update strategies in the smaller surface.
- Toggling "show support values" is `internalNames(true)` - internal
  node *names* are what a real newick file's own bootstrap/support
  values usually are, and this option is literally "show internal node
  names as labels" (confirmed by tracing `showInternalName()`'s own
  source, not the option's name alone) - there's no separate
  "bootstrap" boolean that does anything (present in the options
  defaults object but never read anywhere else in the bundle - checked
  directly, a dead/vestigial key in this version).
- There's no built-in option at all for showing branch-length *text*
  labels (only the "branches"/"scaling" options controlling whether
  branch length *proportionally sizes* the drawn geometry, which was
  already the default behavior regardless of this toolbar) - though the
  library does already show a native SVG `<title>` hover tooltip
  ("Length = ...") on every edge, found while tracing this, independent
  of anything added here. The toolbar's own visible text label is a
  small custom `style_edges` callback instead, using
  `this.phylotree.branch_length_accessor(edgeDatum.target)` - the
  library's own accessor (matching how it computes that hover tooltip
  value) - rather than reading the parsed newick's raw `attribute`
  field by hand.
- **A real bug in this library version's `refresh()`**, found by hitting
  it directly (a synthetic reproduction, not a docs warning): its
  edge-styling loop calls the configured `style_edges` callback via
  `edges.each(d => { this.edge_styler(select(this), d); })` - an *arrow*
  function, so `this` inside stays bound to the `TreeRender` instance
  (from `refresh()`'s own enclosing scope) instead of d3's usual
  per-element `this` rebinding for `.each()` callbacks, meaning
  `select(this)` selects the `TreeRender` object, not a DOM node -
  `element.selectAll(...)` inside the callback then throws
  `this.querySelectorAll is not a function`. `layout()`/`update()` call
  the *same* styler through a different, unaffected code path (a plain
  function, not an arrow function, at a different call site - checked
  directly) where `container` is already a proper d3 selection - so the
  branch-length toggle calls `display.layout()`, never `display.refresh()`,
  specifically to avoid this.

Verified end to end in headless Chrome (this app's real vendored static
assets): every toggle on both surfaces, in combination, including that
the inline preview's toolbar survives a simulated poll (the wrapper
gets reparented into a fresh anchor, same as `placeTreePreviewIntoLiveRegion`
already did before this) without losing its already-built tree. No new
Python tests - same reasoning as the section above, this is pure
client-side behavior with no server-side/context-variable change.

### Three real bugs found from actually using the toolbar work above

All three reported live, against a real tree on a real local instance -
none of them were caught by the harness-based verification the toolbar
work above relied on, since none of them are exceptions or missing
DOM nodes; they're all "the thing technically works, but behaves badly
in a way only a real user actually interacting with it would notice."

- **Text far too large relative to the tree in the small inline
  preview.** Confirmed directly (headless Chrome, a real container
  narrower than the SVG's own 800-unit logical width, matching a
  realistic preview panel): phylotree.js 2.x's "responsive" mode does
  scale the tree's *geometry* (nodes, branches) down to fit a smaller
  container via `viewBox` + CSS width, but does **not** scale the text
  down along with it - a node label's `getComputedStyle().fontSize`
  and rendered `getBoundingClientRect().height` both stayed at the
  raw, unscaled `"font-size"` option value (14, the library's own
  default) regardless of how much smaller the SVG was actually
  displayed. In the full-page viewer (plenty of room, no `responsive`
  scaling in play) this never showed up. Fixed by giving the inline
  preview's own `window.treePreviewOptions` a smaller
  `'font-size': 9` - the full viewer's `treeOptions` is untouched,
  since it isn't the "small" one.
- **Click-and-drag getting stuck in "some kind of selection mode."**
  Both `brush` (click-and-drag to rubber-band-select a clade) and
  `zoom` (click-and-drag to pan) were enabled together, on both
  surfaces - the library's own default for `brush` is `true`, which
  this toolbar work had (silently, via not overriding it) inherited.
  The two features bind to the exact same plain-drag gesture, and nothing
  about a brush selection offers a way to back out of it - which is
  exactly the "stuck" feeling reported live. Fixed by setting
  `brush: false` on both `treeOptions` and `window.treePreviewOptions`
  - a drag is now unambiguously "pan," a plain click is unambiguously
  "select/open this node's menu." Verified directly (headless Chrome,
  dispatched a real mousedown/mousemove/mouseup drag sequence): with
  `brush: false`, no `.tree-selection-brush` element exists at all, and
  the same drag correctly updates `display.currentZoomTransform` (a
  clean pan) instead.
- **Left-click on a node did nothing in the small preview, but worked
  in the full viewer - and right-click, which this page's own hint text
  suggested, just opened the browser's native menu.** The second half
  turned out to be a documentation bug in this codebase, not a library
  bug - see this section's own "Verified directly, not assumed"
  paragraph above for the fix (the menu is genuinely wired to plain
  `"click"`, never `"contextmenu"`, confirmed directly in the vendored
  source). The first half was never conclusively separated from the
  `brush`/`zoom` gesture conflict above - a drag-vs-click
  disambiguation problem (which the browser resolves by movement
  distance/timing) is a very plausible explanation for "clicking felt
  like it did nothing" on a small, densely-packed preview where a real
  mouse is more likely to move a pixel or two between mousedown and
  mouseup than on a large, spaced-out tree - and disabling `brush`
  removes the competing gesture entirely. Verified directly (headless
  Chrome, `elementFromPoint` at a node's real on-screen center plus a
  dispatched mousedown/mouseup/click sequence) that a left-click,
  post-fix, correctly opens `.phylotree-context-menu` with real content
  ("Collapse Subtree," a selection-toggle section, etc.) on the full
  viewer - this was, at the time, only actually checked on that one
  surface, not the inline preview too (see immediately below for why
  that gap mattered).

No new Python tests - all three are pure client-side interaction
behavior with no server-side/context-variable change, same reasoning as
the sections above.

### A fourth bug, found immediately after: left-click still did nothing on the inline preview specifically

Reported live right after the three fixes above shipped - the built-in
menu now genuinely opened on the full viewer (confirmed above) but
still silently did nothing on the inline preview, on the same real
tree. Root cause was a real bug in phylotree.js 2.6.0 itself, found by
reading `nodeDropdownMenu`'s own source line by line rather than
guessing: every other DOM lookup in that function correctly uses
`select(container)` (d3's own helper, which accepts a CSS selector
string, a raw DOM node, or an existing d3 selection interchangeably) -
except the one line that actually makes the menu visible, which instead
calls `document.querySelector(container)` directly:

```js
let tree_container = document.querySelector(container); // eslint-disable-line
let rect = tree_container.getBoundingClientRect();
menu_object.style("position", "absolute")....style("display", "block");
```

`document.querySelector()` only ever accepts a selector **string** -
the inline preview's own `rebuildTreePreviewSvg()` was passing
`container: $container[0]` (the actual `<div>` element, not a
selector), which the full viewer's `render({container: "#tree_container", ...})`
never did. Passed a raw DOM element, `document.querySelector()` throws
(the element coerces to the string `"[object HTMLDivElement]"`, not a
valid selector) - an uncaught exception inside a plain d3 `.on("click", ...)`
handler, with nothing higher up the call chain to catch it, silently
aborting `nodeDropdownMenu` on the line *before* `.style("display",
"block")` ever runs. The menu `<div>` itself still gets created earlier
in the same function (a different code path, using the tolerant
`select(container)`), which is exactly why it showed up as "the menu
element exists in the DOM but stays display:none forever" when this
was actually diagnosed (headless Chrome: dispatched a real click at a
node's `elementFromPoint` target, then inspected
`.phylotree-context-menu`'s own computed `display` and
`getBoundingClientRect()` directly, rather than assuming from the
click alone) - not "nothing happens at all" the way it reads from the
outside.

Fixed by passing `container: '#tree-preview-svg'` (a plain selector
string, matching the full viewer's own convention) instead of the raw
element - `#tree-preview-svg` is a unique id in the outer,
never-duplicated shell (see this file's earlier tree-preview sections
for why that's safe), so a string selector is unambiguous here.
Verified directly, the same way the bug itself was found: after the
fix, the same dispatched click makes `.phylotree-context-menu` compute
to `display: block` with a real, correctly-positioned `getBoundingClientRect()`,
on the inline preview specifically. The full toolbar regression (every
toggle, in combination, plus the poll-reparent survival check) was
re-run afterward too, to confirm this one-line change didn't disturb
anything else already verified.

### Copy selected tip names, and export the tree as SVG/PNG

Two more toolbar buttons, on both surfaces: "Copy selected tip names"
and "Export SVG"/"Export PNG". Selection itself is entirely the
library's own existing mechanism (a click on a node opens its menu -
"All terminal nodes"/"Incident branch"/"Path to root" under "Toggle
Selection" already call `display.modifySelection(...)` internally, from
the toolbar work above) - nothing new was built for *making* a
selection, only for *reading it back*.

- **Copy tip names** reads `display.getSelection()` (a real method,
  confirmed directly in the vendored bundle - returns every currently-
  selected node, leaf or internal) and expands any selected *internal*
  node to its own leaf descendants via `display.phylotree.
  selectAllDescendants(node, true, false)`, so selecting a clade's root
  from the menu and selecting its tips one by one both copy the same
  tip name list. Copies via `navigator.clipboard.writeText()`, with a
  fallback (and a fallback-on-*rejection*, not just on the API being
  entirely absent) to the same temporary-textarea-plus-`execCommand`
  approach the Workspace permalink button already uses elsewhere in
  this codebase, for consistency.
- **SVG/PNG export** clones the live `display.svg.node()`, sets explicit
  `width`/`height`/`viewBox` on the clone (needed for the inline
  preview specifically, which renders "responsive" - no fixed
  width/height attributes on the live SVG to begin with), and - the one
  non-obvious part - embeds `phylotree.css`'s own text inside a
  `<style>` element in the exported SVG. Without this, an exported file
  would come out unstyled (plain black shapes, no node colors) - a
  serialized SVG fragment doesn't carry the *page's* external
  stylesheet with it. Fetched once (same-origin, so no CORS concern in
  real use) and cached rather than hand-copied into these templates, so
  it can't drift out of sync with whatever's actually vendored. PNG
  export then rasterizes that exported SVG onto a `<canvas>` (2x scale
  for a sharper export) and calls `canvas.toBlob(..., 'image/png')` -
  the `<img>` used to draw onto the canvas is loaded from a base64
  `data:` URI, not a `URL.createObjectURL()` blob: URL, specifically to
  avoid a known class of browser bug where a blob: URL can leave the
  canvas "tainted" (blocking `toBlob()`) even for same-document,
  same-origin SVG content that a data: URI doesn't hit.

**A second real bug found while building this, more serious than the
branch-length-toggle-specific one already documented above:**
`display.modifySelection()` - the library's *own* selection mechanism,
called internally by the built-in menu's "All terminal nodes"/etc.
actions, i.e. the exact feature "copy tip names" depends on - also
calls the same buggy `refresh()` internally. Since this toolbar's own
`branchLengthStyler`/`treePreviewBranchLengthStyler` is registered via
`style_edges()` unconditionally (not just while the branch-length
toggle is on - the callback runs on every edge regardless, and checks
the toggle itself as its first real line of logic), simply *having*
that toggle installed was enough to make selecting anything via the
built-in menu throw and break, independent of whether branch lengths
were ever shown at all. Caught live, building this exact feature - not
by the earlier toolbar verification, since that never exercised
`modifySelection()` (only this toolbar's own `.layout()`-based toggle,
which was already known to avoid `refresh()`).

Fixing this took two attempts: the first guard checked
`typeof element.selectAll !== 'function'`, which still threw -
`select()` (d3's own selection constructor) always returns a full
Selection object with every method present, *regardless* of what it's
actually wrapping, so `element.selectAll` is a function either way; it
only throws **inside** `selectAll()`, once it tries to call
`querySelectorAll` on the underlying wrapped object. The working guard
instead checks the *wrapped node itself* -
`element.node() && typeof element.node().querySelectorAll ===
'function'` - which correctly distinguishes a real DOM-wrapping
selection from `refresh()`'s broken one (wrapping the `TreeRender`
instance). Applied to both `branchLengthStyler` (full viewer) and
`treePreviewBranchLengthStyler` (inline preview) identically. Verified
directly (headless Chrome): a real `display.modifySelection(tips)`
call, mirroring exactly what the built-in menu's "All terminal nodes"
action does, no longer throws, and `getSelection()` correctly reflects
the selection afterward - on both surfaces.

Verified end to end (headless Chrome, this app's real vendored static
assets, the real vendored `phylotree.css` content - not a stand-in
string - fetched once via `curl` and mocked into `$.get` since a
`file://` test page can't itself cross-origin-fetch it the way a real
same-origin deployment does): a real clade selection, its tip names
correctly extracted (including the internal-node-selection expansion
case), a captured SVG export blob containing both real tree markup and
the embedded CSS, and a captured PNG export blob with a real, non-trivial
size and no canvas-tainting error - on both surfaces, plus the full
existing toggle regression re-run afterward to confirm nothing else
regressed. **Not verified** (a real limitation of headless/synthetic-
click automation, not a known bug): whether `navigator.clipboard.
writeText()` itself actually places text on a real system clipboard -
`navigator.clipboard` is a read-only browser property that can't be
mocked by simple assignment, and a script-triggered synthetic click
doesn't carry the "real user gesture" credit the Clipboard API
requires, so the write silently no-ops in this exact test setup
regardless of whether the underlying code is correct. Worth a real
manual click to confirm.

### A third layout: "unrooted" (distinct from the existing "circular"/radial one)

The toolbar's layout control was a 2-way `is-radial` toggle
("Linear"/"Radial"). Requested: a third, genuinely different "unrooted"
layout - the name is a real source of confusion here, since an
equal-angle unrooted layout is itself sometimes informally called
"radial" too, easy to conflate with phylotree.js's own `is-radial`
option (a circular layout rooted at the center - a different algorithm
entirely). Confirmed directly against the vendored bundle (not assumed
from the name alone) that the library has a real, separate, working
third option for exactly this: `is-unrooted` / the `unrooted()`
accessor, genuinely read in several real layout-computation code paths
(`showInternalName`-adjacent direction logic, the core layout function,
etc.) - not a vestigial/unused option the way `bootstrap` turned out to
be earlier in this file. Verified empirically too: rendering the same
tree with `is-radial` vs `is-unrooted` produces different node
transform coordinates, not the same layout under two names.

The old "Linear"/"Radial" 2-button group became a 3-button
"Linear"/"Circular"/"Unrooted" group on both surfaces (renaming
"Radial" to "Circular" specifically to stop colliding with the new
"Unrooted" button's own informal "radial" nickname) - `is-radial` and
`is-unrooted` are never both true at once (nothing in the library
enforces that itself; `setLayout()`/`setTreePreviewLayout()` always set
both flags together, one true and one false, whichever button was
clicked). Same rebuild-from-scratch pattern as every other layout-
affecting toggle already in this toolbar. Also added proper `active`-
class tracking to the full viewer's layout buttons for the first time
(the inline preview's other toggle buttons already did this; the full
viewer's plain Linear/Radial pair never had it, since with only two
mutually exclusive states the default/non-default one was implied -
worth doing properly now that there are three).

Verified in headless Chrome (both surfaces): clicking each of the three
layout buttons in sequence correctly sets exactly one of
`is-radial`/`is-unrooted` and clears the other, moves the `active` CSS
class to the clicked button and off the other two, and the rest of the
existing toggle regression (align tips, support values, branch lengths,
zoom, reset, the inline preview's poll-reparent survival) still passes
unchanged afterward.

### A negative branch length breaks the new unrooted layout specifically

Reported live against a real tree with `Tant_cDNA:-0.00007` - a small
negative branch length, the kind of numerical-correction artifact some
phylogenetic methods leave near true zero (not a data error to reject,
just not physically meaningful as a length). In the unrooted layout
specifically, that one branch rendered roughly **800 SVG units long**,
against a tree where the next-longest real branch was ~230 units and
most branches were single digits - a real, large visual defect, not a
subtle one. Linear and circular layouts were unaffected - specific to
unrooted's own equal-angle direction-vector math.

Root cause, found by measuring actual rendered edge lengths and
`branch_length_accessor()`'s return value directly against the real
tree (not guessed): `branch_length_accessor` clamps a negative length to
exactly `0` - and it's that **zero**, not the negative sign itself, that
degenerates the unrooted layout's per-edge direction vector (confirmed
by testing a tiny *positive* epsilon in place of `0`, which did **not**
fix it, before finding what does). What actually fixes it: stripping
just the minus sign and keeping the branch's original magnitude
(`-0.00007` -> `0.00007`) - verified directly, both offline and against
the real live page (headless Chrome, clicking the real "Unrooted"
toolbar button and re-measuring: the tree's longest edge afterward was
~313 units, matching its real longest branch, with no ~800-unit
outlier).

Fixed by sanitizing the newick string itself, once, right where each
surface first receives it - `templates/treeviz/tree.html`'s `var
newick = "...".replace(/:-(?=[\d.])/g, ':')`, and the inline preview's
equivalent applied once when `window.treePreviewPendingNewick` is first
set (not repeated on every `rebuildTree()`/`rebuildTreePreviewSvg()`
call) - rather than relying on phylotree.js's own clamping, which by
the time it runs has already fed the broken zero into the layout math.
Applied unconditionally (not just when the unrooted layout is active) -
a negative branch length isn't meaningful in the *other* two layouts
either, even though they happen to tolerate it without visibly
breaking; sanitizing once regardless of which layout is currently
selected is simpler than trying to apply it only when needed.

No new Python tests - this is a pure client-side newick string
transform with no server-side/context-variable change, same reasoning
as this file's other tree-toolbar sections. Verified end to end
(headless Chrome): the full existing toggle regression (layout
switching between all three, align tips, support values, branch
lengths, zoom, reset) still passes unchanged with the sanitizing in
place, on both surfaces.

### Tip label font size +/- controls

Two more buttons, both surfaces: increase/decrease the tip label font
size (2px per click, clamped to [6, 40]). Straightforward-looking, but
the actual mechanism needed real investigation, not just passing a
bigger number to an obvious option:

- **The library's own `font_size()` accessor method is completely
  unreachable.** Confirmed directly (not guessed): the `TreeRender`
  constructor sets `this.font_size = this.options["font-size"]` - a
  plain instance property - which permanently shadows the same-named
  prototype method mixed in elsewhere in the bundle. `display.
  font_size(20)` throws `"display.font_size is not a function"`,
  because property lookup finds the instance's own number (14) long
  before it would ever reach the prototype's method.
- **Setting the `"font-size"` render option alone doesn't make an
  *increase* visible either.** Confirmed empirically: rendering with
  `"font-size": 24` still measured a computed tip-label font-size of
  14px. The actual rendered size
  (`this.shown_font_size`) is computed as `Math.min(this.font_size,
  this.scales[0])`, and under this library's default spacing mode
  ("fixed-step"), `this.scales[0]` is just `this.fixed_width[0]` -
  hardcoded to `14` in the constructor, entirely independent of the
  "font-size" option. A *decrease* (e.g. `"font-size": 8`) still worked
  fine on its own, since 8 is already below that same 14 cap either way
  - which is what made the bug easy to miss testing only in one
  direction.
- **The real lever is the library's own live `spacing_x()` accessor**
  (unaffected by the `font_size()` shadowing bug - nothing sets
  `this.spacing_x` as an instance property) - it maps directly onto
  `fixed_width[0]`, i.e. the same value that caps `shown_font_size`.
  Raising it past the desired font size unlocks that size. Its own
  built-in re-render (calling `this.placenodes()` when set) wasn't
  enough on its own to recompute `shown_font_size` either, confirmed by
  testing - only a full explicit `.layout()` call afterward actually
  did. `applyTipFontSize()` (full viewer) /
  the inline preview's equivalent block in `rebuildTreePreviewSvg()`
  both: set `treeOptions["font-size"]`, do the normal full rebuild, then
  call `display.spacing_x(fontSize, true)` + `display.layout()` only
  when the current `spacing_x()` is smaller than the new font size (a
  no-op skip on a pure decrease, since a bigger spacing_x left over from
  an earlier increase is harmless - `shown_font_size`'s own `Math.min`
  already caps it back down correctly via the smaller "font-size" value
  itself).
- Raising `spacing_x` also means more vertical space between sibling
  tips generally, not just a bigger font cap - an intentional,
  reasonable side effect (bigger text needs more room to avoid
  overlapping labels), not a workaround for something unwanted.

Verified end to end (headless Chrome, both surfaces): clicking increase
three times then decrease four times from the same starting point
produced the expected net-one-decrement font size at every step (14 ->
16 -> 20 -> 12 on the full viewer; 9 -> 11 -> 15 -> 7 on the inline
preview, matching its own smaller default), each time reading the
*actual computed* CSS font-size off a real rendered tip label, not just
trusting the option value was set. The full existing toggle regression
was re-run afterward on both surfaces too. No new Python tests - same
reasoning as this file's other tree-toolbar sections, pure client-side
behavior with no server-side/context-variable change.

### Width/height range sliders - and a real architectural finding they forced

Two `<input type="range">` controls, both surfaces, for the tree's own
`width`/`height` render options. Sounds like it should be as simple as
wiring two sliders to `treeOptions.width`/`.height` and rebuilding - it
wasn't, because of a real, significant finding about how this library
actually sizes a tree:

**Under phylotree.js 2.x's own *default* spacing mode
("left-right-spacing"/"top-bottom-spacing": "fixed-step"), the
"width"/"height" render options are completely ignored for layout
sizing.** Confirmed by reading `do_lr()`'s and the equivalent
top-bottom logic directly: under fixed-step, `this.size[0]` (leaf/
vertical extent) is computed as `this._extents[0][1] *
this.fixed_width[0]` and `this.size[1]` (depth/horizontal extent) as
`this.max_depth * this.fixed_width[1]` - both driven purely by the
tree's own shape (leaf count, depth) and `fixed_width`
(`spacing_x()`/`spacing_y()`, the same accessors the tip-font-size
feature above already uses), **never** by `this.width`/`this.height`.
Verified empirically too: doubling both `width`/`height` options and
rebuilding left the tree's own rendered extent (measured directly off
node transform coordinates) completely unchanged. Only the *other*
branch of that same logic - `"fit-to-size"` mode, which computes
`this.scales[0]`/`[1]` from `(this.size[0/1] - padding) /
extent[0/1]`, i.e. actually starting from the width/height option - ties
rendered size to those options at all.

Both `treeOptions`/`window.treePreviewOptions` were switched to
`"left-right-spacing": "fit-to-size", "top-bottom-spacing":
"fit-to-size"` for this reason - not something the width/height
sliders opt into per-use, a permanent mode change, since there's no
sane way to offer working sliders without it. Verified this doesn't
change the *rendered look* of a normally-sized tree in any of the three
layouts, only how its size responds to width/height.

**Genuinely nice side effect: this made the tip-font-size feature
simpler.** The `spacing_x()`-bumping workaround `applyTipFontSize()`
needed (see the section above) was specifically compensating for
fixed-step mode's `scales[0] = fixed_width[0]` cap - under
`fit-to-size`, `scales[0]` comes from the width/height option instead,
which for any reasonably-sized canvas is comfortably larger than a
normal tip font size. Confirmed empirically: `"font-size": 24` alone,
with no `spacing_x()` involvement at all, now measures as a real 24px
computed font-size. `applyTipFontSize()`/the inline preview's
equivalent were simplified accordingly - the workaround is gone
entirely, not just unused.

For the **inline preview specifically** (the one surface using
`responsive: true`), the width/height sliders behave a bit differently
from a literal "resize the picture" control: since CSS forces the
displayed width to 100% of its container regardless of the `width`
option, changing "Width" doesn't change the pixel width shown - it
changes the tree's internal aspect ratio (`viewBox="0 0 width
height"`), i.e. how tall the preview ends up once that fixed width is
applied via `height: auto`. "Height" has a more direct effect (visibly
taller/shorter, matching the ratio). Confirmed directly (headless
Chrome, an artificially narrowed container matching a realistic preview
column): the *displayed* width stayed capped at the container's own
width in every case, while displayed height scaled exactly with the
viewBox aspect ratio (e.g. 400x900 -> displayed 700x1575, matching
400:900 = 700:1575). Both sliders are still offered - "Width" is a real
tuning knob (how much horizontal vs. vertical room the layout algorithm
gets to work with), just not literally "how many pixels wide it looks."

Verified end to end (headless Chrome, both surfaces): the sliders
genuinely grow/shrink the tree's own rendered extent (measured directly
off node coordinates, not just the SVG's own declared size) - roughly
2x width option -> 2x rendered extent, consistent with `fit-to-size`'s
own scale computation - and the *entire* existing regression suite
(layout switching between all three, align tips, support values,
branch lengths, zoom, font size, selection/copy-tip-names, SVG/PNG
export, the inline preview's poll-reparent survival, and the negative-
branch-length unrooted-layout fix re-checked against the real live
page) was re-run after this mode switch and still passes unchanged -
this touches core layout sizing shared by every other feature in this
toolbar, so all of it needed re-confirming, not just the two new
sliders. No new Python tests - pure client-side behavior, same
reasoning as this file's other tree-toolbar sections.

### Width/height sliders should resize the tree, not the page

Immediate follow-up bug report: the height slider was growing the
*visible panel* on the page along with the tree, not just the tree's
own content - on the full viewer, `#tree_container` (non-responsive
mode) had no CSS height at all, so its box just grew to whatever
`svg.attr("height", ...)` set (see `set_size()`/`initialize_svg()`
above); on the inline preview (`responsive: true`), the svg's own style
is `width: 100%; height: auto`, so a taller `height` option widens the
`viewBox`'s aspect ratio, which - with width pinned to 100% of the
container - renders proportionally taller on the page too. Either way,
the container was sized *by* its content instead of the other way
around, so "scale the tree" and "grow the page" were the same thing.

Fixed by giving both containers a genuinely fixed CSS `height` plus
`overflow: auto` - `#tree_container { height: 800px; }` (matching
`treeOptions`'s own default `height: 800`, so the tree fits with no
scrollbar at its starting size) and `#tree-preview-wrapper .panel-body
{ height: 420px; }` (matching `window.treePreviewOptions`'s default
`height: 420`, same reasoning). Past that starting size, a slider still
genuinely grows the tree's own rendered extent exactly as before (see
the section above) - it just now overflows within a viewport of fixed
visible size, scrollable via the browser's native scrollbars, rather
than pushing the rest of the page down. This is a pure CSS containment
fix - no JS changed, the sliders' own handlers (`rebuildTree()`/
`rebuildTreePreviewSvg()`) are untouched.

Verified directly (headless Chrome, both surfaces, reading
`getBoundingClientRect()` on the container/panel-body before and after
driving each slider to its max): container height/width are bit-for-bit
identical before and after in every case, `scrollHeight`/`scrollWidth`
exceed `clientHeight`/`clientWidth` once a slider is pushed well past
its default (confirming the content is actually overflowing rather than
silently clipping), and shrinking a slider back down still leaves the
container at its original fixed size throughout. No new Python tests -
pure CSS, same reasoning as this file's other tree-toolbar sections.

### Tip font size +/- buttons replaced with a range control

Both surfaces' "Tip font size" +/-2px buttons (`font_size_increase`/
`font_size_decrease`, `tp_font_size_increase`/`tp_font_size_decrease`)
were swapped for a single `<input type="range">`
(`tree_font_size_range`, `tp_font_size_range`, both `min="6" max="40"
step="1"`, matching the buttons' own existing clamp range) - the same
control shape the width/height sliders already established on both
surfaces. `applyTipFontSize()`'s own body (full viewer) and the inline
preview's equivalent inline handler are unchanged - only how
`tipFontSize`/`window.treePreviewOptions['font-size']` gets set
changed, from `Math.min/max(current +/- 2, ...)` on a click to
`parseInt(this.value, 10)` on the range's own `'input'` event (same
"fires continuously while dragging" choice already used for the width/
height sliders, for the same immediate-feedback reason).

No behavior regression from the switch: `shown_font_size`'s own
`Math.min(this.font_size, this.scales[0])` cap (see the tip-font-size
section above) still applies exactly the same way regardless of how
`"font-size"` gets set - confirmed directly (headless Chrome, both
surfaces): dragging the slider to its max (40) on a tree/container
whose current `fit-to-size` `scales[0]` is smaller than 40 still caps
the *computed* rendered size at that smaller value, same as the old
+/- buttons would have if clicked enough times to reach the same
option value. The full existing toggle regression (layout switching,
align tips, support values, branch lengths, zoom, width/height sliders,
the panel-stays-fixed-size CSS above, selection/copy-tip-names, SVG/PNG
export, the inline preview's poll-reparent survival) was re-run
afterward on both surfaces and still passes unchanged. No new Python
tests - pure client-side control swap, same reasoning as this file's
other tree-toolbar sections.

### Width/height/font sliders weren't actually independent on the inline preview

Reported live: dragging the Width slider on the workspace/history
preview visibly changed the tree's height and font size too, not just
its horizontal extent - despite the width handler only ever touching
`window.treePreviewOptions.width`. The full viewer never had this bug
(confirmed directly, not assumed - see below); it was specific to the
inline preview's own `responsive: true` mode.

**Root cause: `responsive: true` renders via an svg `viewBox`, not real
width/height attributes, and viewBox-to-viewport scaling is isotropic
(the same factor on every axis) whenever the *displayed* size is
externally pinned - which it is here** (`style("width", "100%")
.style("height", "auto")`, see `initialize_svg()` in the vendored
bundle). With the container's physical width fixed, widening the
`viewBox` (i.e. raising the "width" option, which sets `viewBox="0 0
width height"`) necessarily shrinks the effective pixels-per-unit
factor - on *both* axes, since a viewBox has one scale, not two
independent ones - so the displayed height shrinks too, and with it
every element sized in svg user units, including tip label text.
Confirmed directly (headless Chrome, measuring real `getBoundingClientRect()`
values, not the svg's own internal viewBox numbers which don't reveal
this): doubling only the "width" option roughly halved the *displayed*
height of the same tree, with "height"/"font-size" never touched.
`getComputedStyle().fontSize` on a text element is *not* a reliable way
to check this, incidentally - it reports the authored/attribute value,
not the effective on-screen size after viewBox scaling, which is
exactly why this file's own earlier verification of the font-size
slider (see that section above) didn't already catch this - it happened
to only ever test font size in isolation, never alongside a width
change.

**Fixed by switching `window.treePreviewOptions.responsive` to `false`**
- matching `templates/treeviz/tree.html`'s own full viewer, which
renders with real `width`/`height` svg attributes and no viewBox at
all, and was never affected by this in the first place (confirmed
directly: an isolated width-only change against the vendored bundle in
non-responsive mode left the rendered vertical extent and `scales[0]`/
`shown_font_size` completely unchanged). With no viewBox, "width" and
"height" size their own axis independently by construction - a tree
wider or taller than the fixed-size panel (see the section above) now
genuinely overflows into that panel's own scrollbars instead of being
silently rescaled to fit, which is also exactly the "the tree can
potentially be larger than the viewport" behavior asked for alongside
this fix. `buildTreePreviewExportableSvg()` (SVG/PNG export) needed no
change - it already set `width`/`height`/`viewBox` explicitly on the
exported clone regardless of the live svg's own mode.

**A second, unrelated bug found while verifying this fix, pre-existing
and not introduced by it:** `tp_height_range`'s declared `value="420"`
was never actually reachable given its own `min="200" step="50"` (valid
steps are 200, 250, 300, ... - 420 isn't one of them) - the HTML range
input spec's value-sanitization algorithm silently snaps an
out-of-step value to the nearest valid one *whenever it's set*,
including the initial parse of the `value` attribute on page load, not
just on drag. So this slider's real starting value was always `400`,
one silent step away from the `420` both its own markup and
`window.treePreviewOptions`'s own default height claimed - invisible in
normal use (dragging from an already-snapped position never re-triggers
the mismatch), but exactly the kind of thing that surfaces as a
confusing "off by one step" once something (here, this fix's own
regression test, setting `.val(420)` programmatically to reset between
checks) tries to treat 420 as the slider's real value. Fixed by
widening its `step` to `10` (200, 210, 220, ... does include 420) -
the width slider and both surfaces' other range inputs were checked
too and don't have this misalignment.

Verified end to end (headless Chrome): width, height, and font-size
changes, applied one at a time, each leave the other two - and the
fixed-size panel from the section above - completely unaffected
(`getBoundingClientRect()` measurements before/after every change), and
a tree pushed well past its default size genuinely overflows
(`scrollHeight`/`scrollWidth` exceeding `clientHeight`/`clientWidth`)
rather than rescaling. The full existing toolbar regression (layout
switching, align tips, support values, branch lengths, zoom, the
poll-reparent survival check) was re-run afterward and still passes
unchanged. No new Python tests - pure client-side rendering-mode
change, same reasoning as this file's other tree-toolbar sections.

### Home page: download banner for the Firefox browser extension

A separate, external project
(`C3BI-pasteur-fr/ngphylogeny-browser-extension`) ships a Firefox
extension that lets a user paste sequences from any web page straight
into an NGPhylogeny One Click or BLAST form - this app has no code of
its own for it, only a place to point people at the release. Added a
promotional banner to the home page (`templates/home.html`, right after
the existing `{% include "include/phylogeny_analysis_list.html" %}`)
linking directly to the versioned GitHub release asset
(`.../releases/download/v0.2.0/....xpi`) - Firefox's own
`.xpi`-extension handling shows its native "Add extension?" install
prompt on navigating straight to that URL, no App-side download view or
storage needed.

Deliberately **not** added to `templates/include/phylogeny_analysis_list.html`
itself, even though the banner sits directly below that include's own
markup - that partial is also reused, unmodified, by
`templates/phylogeny_analysis_choices.html` (a separate "Phylogeny
Analysis" page reachable from elsewhere), and the ask was specifically
for the home page. Keeping the banner in `home.html`'s own `{% block
content %}` instead of inside the shared include means the other page
using that include doesn't pick up a banner nobody asked to put there.

Went through this session's usual "propose the UI change as a mockup
built from the page's real layout/colors before touching templates"
step first - two Artifact-canvas options (a fourth card slotted into
the existing One Click/Advanced/A la Carte row, vs. a separate banner
below that row) were built matching the real page's actual accent
color (`#1f77b4`), card gradient/border treatment, and
`"Helvetica Neue", Helvetica, Arial, sans-serif` font stack (all read
directly from `assets/css/home.css`/`custom.css`, not guessed). The
banner option was chosen - it reads as a distinct callout for a
one-off external tool rather than a fourth, visually-equal workflow
choice competing with the three core submission paths, and doesn't
force `phylo-analysis`'s existing 3-wide grid to become an uneven
4-wide one.

New CSS (`assets/css/home.css`): `.extension-banner` (a
left-accent-bordered callout, same `#1f77b4` accent and light-blue tint
as the mockup, not reusing `.phylo-analysis`'s own card styling since
this isn't one of those cards) plus a generic `.badge-new` pill (amber
`#e8871e`, matching the mockup's own "NEW" tag) - the closest existing
precedent in this codebase for a colored pill is the workspace history
page's `.status-pill-*`/`.itol-badge` set (see the step-chain/dataset-
table sections above), but those are state-specific (finished/running/
error) and not reusable here, so a small new class was added instead
of stretching an unrelated one.

No new Python tests - a static link + CSS addition on an already-
covered page (`data.tests.StaticPagesSmokeTest.test_pages_return_200`
already renders `home.html` end to end and continues to pass unchanged);
nothing here has server-side logic or a context variable to unit test,
same reasoning as this file's other pure-template/CSS sections.
