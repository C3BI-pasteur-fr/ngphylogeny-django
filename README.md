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

* Run the django server

```
python manage.py runserver
```
