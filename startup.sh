#!/bin/bash

USERPASS=$1
USEREMAIL=$2
GALAXYSERVER=$3
GALAXYKEY=$4

# Wait for galaxy
echo ${GALAXYSERVER}
echo $(curl -s ${GALAXYSERVER}/api/version)
until [[ "$(curl -s ${GALAXYSERVER}/api/version | jq -r '.version_major')" == "18.09" ]] ; do
    echo ${GALAXYSERVER}
    echo $(curl -s ${GALAXYSERVER}/api/version)
    sleep 2
done

# Wait for database server if any
if [ ! -z "$NGPHYLO_DATABASE_HOST" ]
then
    until PGPASSWORD=$NGPHYLO_DATABASE_PASSWORD psql -h "$NGPHYLO_DATABASE_HOST" -U "$NGPHYLO_DATABASE_USER" -c '\q'; do
	>&2 echo "Postgres is unavailable - sleeping"
	sleep 1
    done
fi

# Start required services 
service redis_6379 start

# Initialize databases
python manage.py makemigrations
python manage.py migrate
#python manage.py migrate --run-syncdb
python manage.py createcachetable
python manage.py collectstatic --noinput

# Create admin user
echo "from django.contrib.auth.models import User; User.objects.filter(username='admin').exists() or User.objects.create_superuser('admin', '$USEREMAIL', '$USERPASS');exit()" | python manage.py shell
python manage.py loaddata tool_flags
python manage.py creategalaxyserver --url=$GALAXYSERVER --activate
# Adding shared galaxy api key
python manage.py addgalaxykey --user admin --galaxyurl $GALAXYSERVER --galaxykey $GALAXYKEY
python manage.py importtools --galaxyurl=$GALAXYSERVER --query="phylogeny" --flags=toolflags.txt --inputfields=toolfields.txt
python manage.py import_links --linkfile=toollinks.txt
python manage.py importworkflows --wfnamefile=wfnames.txt --galaxyurl=$GALAXYSERVER

export PYTHONPATH=$PWD:$PYTHONPATH

# Re configure environment variables of celery daemons
cp docker/celeryd.default /etc/default/celeryd
if [ ! -z $NGPHYLO_HOST ]
then
    echo "export NGPHYLO_HOST=$NGPHYLO_HOST" >> /etc/default/celeryd
fi

if [ ! -z $NGPHYLO_EMAIL_HOST ]
then
    echo "export NGPHYLO_EMAIL_HOST=$NGPHYLO_EMAIL_HOST" >> /etc/default/celeryd
    echo "export NGPHYLO_EMAIL_PORT=$NGPHYLO_EMAIL_PORT" >>/etc/default/celeryd
    echo "export NGPHYLO_EMAIL_HOST_USER=$NGPHYLO_EMAIL_HOST_USER" >>/etc/default/celeryd
    echo "export NGPHYLO_EMAIL_HOST_PASSWORD=$NGPHYLO_EMAIL_HOST_PASSWORD" >>/etc/default/celeryd
    echo "export NGPHYLO_EMAIL_USE_TLS=$NGPHYLO_EMAIL_USE_TLS" >>/etc/default/celeryd
fi

if [ ! -z $NGPHYLO_DATABASE_ENGINE ]
then
    echo "export NGPHYLO_DATABASE_ENGINE=$NGPHYLO_DATABASE_ENGINE" >>/etc/default/celeryd
    echo "export NGPHYLO_DATABASE_NAME=$NGPHYLO_DATABASE_NAME" >>/etc/default/celeryd
    echo "export NGPHYLO_DATABASE_USER=$NGPHYLO_DATABASE_USER" >>/etc/default/celeryd
    echo "export NGPHYLO_DATABASE_PASSWORD=$NGPHYLO_DATABASE_PASSWORD" >>/etc/default/celeryd
    echo "export NGPHYLO_DATABASE_HOST=$NGPHYLO_DATABASE_HOST" >>/etc/default/celeryd
    echo "export NGPHYLO_DATABASE_PORT=$NGPHYLO_DATABASE_PORT" >>/etc/default/celeryd
fi

service celeryd start
service celerybeat start

#celery multi start 3 -l INFO -c:2 1 -c:3 1 -Q:1 default -Q:2 ncbi_blast -Q:3 monitor --app=NGPhylogeny_fr.celery:app
#celery beat --app=NGPhylogeny_fr.celery:app --loglevel=DEBUG --detach
#python manage.py runserver 0.0.0.0:8000

# Run uwsgi
exec uwsgi --ini /home/ngphylo/docker/ngphylogeny_uwsgi.ini &
# Run nginx
nginx -g 'daemon off;'
