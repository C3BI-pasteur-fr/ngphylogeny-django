# NGPhylogeny.fr


A Django web application to make Phylogeny analysis.

# Installation

NGPhylogeny.fr requires python 2.7.

To install it:

* You may first install and activate a new conda environment:

```
conda create --name ngphylo python=2.7
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

# Docker

You can run NGPhylogeny.fr via its [docker image](https://hub.docker.com/r/evolbioinfo/ngphylogeny/):

## Using a remote galaxy instance:

Note: The remote galaxy instance must contain all necessary tools (listed [here](toolflags.txt)).

```
docker run -p 8080:8000 evolbioinfo/ngphylogeny django_admin_username django_admin_password django_admin_email galaxy_url galaxy_api_key
```

## Using a custom local docker galaxy instance
We provide a [galaxy docker image](https://hub.docker.com/r/evolbioinfo/ngphylogeny-galaxy/) with all tools and workflows already installed (see this [github repo](https://github.com/C3BI-pasteur-fr/ngphylogeny-galaxy)):

```
# Starting Docker image of Galaxy
docker run --privileged=true  \
           -p 8080:80 -p 8121:21 -p 8122:22 \
           evolbioinfo/ngphylogeny-galaxy
# MacOS => Starting NGPhylogeny.fr 
docker run -p 8000:8000 evolbioinfo/ngphylogeny admin admin@admin http://host.docker.internal:8080 admin
# Linux => Starting NGPhylogeny.fr
docker run -p 8000:8000 --net=host evolbioinfo/ngphylogeny admin admin@admin http://localhost:8080 admin
```
