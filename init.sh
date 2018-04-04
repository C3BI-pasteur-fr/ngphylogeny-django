#!/bin/bash

USERPASS=$1
USEREMAIL=$2
GALAXYSERVER=$3

pip install -r requirement.txt
python manage.py makemigrations
python manage.py migrate
echo "from django.contrib.auth.models import User; User.objects.filter(email='$USEREMAIL').delete(); User.objects.create_superuser('admin', '$USEREMAIL', '$USERPASS')" | python manage.py shell
python manage.py loaddata tool_flags
python manage.py creategalaxyserver --url=$GALAXYSERVER --activate
python manage.py importtools --galaxyurl=$GALAXYSERVER --query="phylogeny" --flags=toolflags.txt
python manage.py import_links --linkfile=toollinks.txt
