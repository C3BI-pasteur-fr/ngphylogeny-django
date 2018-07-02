#!/bin/bash

USERPASS=$1
USEREMAIL=$2
GALAXYSERVER=$3
GALAXYKEY=$4

# Start required services 
#service redis_6379 start
#service celeryd start
#service celerybeat start

# Initialize databases
python manage.py makemigrations
python manage.py migrate
# Create admin user
echo "from django.contrib.auth.models import User; User.objects.filter(username='admin').delete(); User.objects.create_superuser('admin', '$USEREMAIL', '$USERPASS');exit()" | python manage.py shell
python manage.py loaddata tool_flags
python manage.py creategalaxyserver --url=$GALAXYSERVER --activate
# Adding shared galaxy api key
echo "from django.contrib.auth.models import User; from galaxy.models import Server, GalaxyUser; u = User.objects.filter(username='admin'); g = Server.objects.filter(url='$GALAXYSERVER'); print(u);print(g);gu=GalaxyUser(user=u.first(),galaxy_server=g.first(),api_key='$GALAXYKEY',anonymous=True);gu.save()" | python manage.py shell
python manage.py importtools --galaxyurl=$GALAXYSERVER --query="phylogeny" --flags=toolflags.txt --force --inputfields=toolfields.txt
python manage.py import_links --linkfile=toollinks.txt
python manage.py importworkflows --wfnamefile=wfnames.txt --galaxyurl=$GALAXYSERVER

export PYTHONPATH=$PWD:$PYTHONPATH
celery --app=NGPhylogeny_fr.celery:app worker --loglevel=INFO &
celery beat --app=NGPhylogeny_fr.celery:app --loglevel=DEBUG &
python manage.py runserver 0.0.0.0:8000

# Run uwsgi
#exec uwsgi --ini /home/ngphylo/docker/ngphylogeny_uwsgi.ini &
# Run nginx
#nginx -g 'daemon off;'
