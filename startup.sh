#!/bin/bash

USERPASS=$1
USEREMAIL=$2
GALAXYSERVER=$3
GALAXYKEY=$4

# Start required services 
service redis_6379 start

# Initialize databases
python manage.py makemigrations
python manage.py migrate
python manage.py createcachetable
# Create admin user
echo "from django.contrib.auth.models import User; User.objects.filter(username='admin').delete(); User.objects.create_superuser('admin', '$USEREMAIL', '$USERPASS');exit()" | python manage.py shell
python manage.py loaddata tool_flags
python manage.py creategalaxyserver --url=$GALAXYSERVER --activate
# Adding shared galaxy api key
python manage.py addgalaxykey --user admin --galaxyurl $GALAXYSERVER --galaxykey $GALAXYKEY
python manage.py importtools --galaxyurl=$GALAXYSERVER --query="phylogeny" --flags=toolflags.txt --force --inputfields=toolfields.txt
python manage.py import_links --linkfile=toollinks.txt
python manage.py importworkflows --wfnamefile=wfnames.txt --galaxyurl=$GALAXYSERVER

export PYTHONPATH=$PWD:$PYTHONPATH

service celeryd start
service celerybeat start

#celery multi start 3 -l INFO -c:2 1 -c:3 1 -Q:1 default -Q:2 ncbi_blast -Q:3 monitor --app=NGPhylogeny_fr.celery:app
#celery beat --app=NGPhylogeny_fr.celery:app --loglevel=DEBUG --detach
#python manage.py runserver 0.0.0.0:8000

# Run uwsgi
exec uwsgi --ini /home/ngphylo/docker/ngphylogeny_uwsgi.ini &
# Run nginx
nginx -g 'daemon off;'
