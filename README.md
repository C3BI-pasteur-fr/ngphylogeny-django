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
pip install -r requirements.txt
```

* Create Django databases

```
python manage.py makemigrations
python manage.py migrate
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
```

* Import existing tools from Galaxy to NGPhylogeny.fr

```
python manage.py importtools --galaxyurl=http://url_galaxy:port --query="phylogeny" --flags=toolflags.txt
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
python manage.py importworkflows --galaxyurl=http://url_galaxy:port
```

* Run the django server

```
python manage.py runserver
```

# Docker

You can run NGPhylogeny.fr using the docker image:

## Building the docker image

```
docker build -t ngphylo .
```

## Running NGPhylogeny.fr using a remote galaxy instance:

Note: The remote galaxy instance must contain all necessary tools (listed [here](toolflags.txt)).

```
docker run -p 8080:8000 evolbioinfo/ngphylogeny django_admin_username django_admin_password django_admin_email galaxy_url galaxy_api_key
```

## Running NGPhylogeny.fr using a custom local galaxy instance
```
# Starting Docker image of Galaxy
docker run --privileged=true  \
       -e GALAXY_CONFIG_TOOL_CONFIG_FILE=config/tool_conf.xml.sample,config/shed_tool_conf.xml.sample,/local_tools/tool_conf.xml \
       -e GALAXY_DOCKER_ENABLED=True -p 8080:80 -p 8121:21 -p 8122:22 \
       evolbioinfo/ngphylogeny-galaxy
# Starting Docker image of NGPhylogeny.fr (MacOS)
docker run -p 8000:8000 evolbioinfo/ngphylogeny admin admin@admin http://host.docker.internal:8080 admin
# Starting Docker image of NGPhylogeny.fr (Linux)
docker run -p 8000:8000 --net=host evolbioinfo/ngphylogeny admin admin@admin http://localhost:8080 admin
```
