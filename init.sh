#!/bin/bash

pip install -r requirement.txt
python manage.py makemigrations
python manage.py migrate
python manage.py createsuperuser
python manage.py loaddata tool_flags
python manage.py creategalaxyserver --interactive --activate
python manage.py importtools --query "phylogeny"
python manage.py compute_tools_links

