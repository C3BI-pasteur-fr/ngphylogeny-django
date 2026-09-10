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
