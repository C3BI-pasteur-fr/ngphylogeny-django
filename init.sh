#!/bin/bash

GALAXYURL=${1-'https://galaxy.pasteur.fr'}

pip install virtualenv
virtualenv ENV
source ENV/bin/activate
pip install -r requirement.txt
python manage.py makemigrations
python manage.py migrate
python manage.py loaddata tool_flags
python manage.py importtools --galaxyurl $GALAXYURL --query "phylogeny"
python manage.py compute_tools_links