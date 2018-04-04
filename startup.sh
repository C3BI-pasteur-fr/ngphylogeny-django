#!/bin/bash

USERPASS=$1
USEREMAIL=$2
GALAXYSERVER=$3

echo "from django.contrib.auth.models import User; User.objects.filter(username='admin').delete(); User.objects.create_superuser('admin', '$USEREMAIL', '$USERPASS');exit()" | python manage.py shell
python manage.py loaddata tool_flags
python manage.py creategalaxyserver --url=$GALAXYSERVER --activate
python manage.py importtools --galaxyurl=$GALAXYSERVER --query="phylogeny" --flags=toolflags.txt --force
python manage.py import_links --linkfile=toollinks.txt
python manage.py runserver 0.0.0.0:8000
