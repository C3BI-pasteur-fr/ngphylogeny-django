#!/bin/bash
# One-shot setup, run by the `init` service in docker-compose.yml before
# web/celery-worker/celery-beat start (they wait on it via
# `depends_on: init: condition: service_completed_successfully`).
set -euo pipefail

if [ -n "${NGPHYLO_DATABASE_HOST:-}" ]; then
    until PGPASSWORD="$NGPHYLO_DATABASE_PASSWORD" psql -h "$NGPHYLO_DATABASE_HOST" -p "${NGPHYLO_DATABASE_PORT:-5432}" -U "$NGPHYLO_DATABASE_USER" -d "$NGPHYLO_DATABASE_NAME" -c '\q' 2>/dev/null; do
        echo "Postgres is unavailable - sleeping"
        sleep 1
    done
fi

# migrations are not committed to git (see CLAUDE.md) - generate them fresh
python manage.py makemigrations
python manage.py migrate
python manage.py createcachetable
python manage.py collectstatic --noinput

# Create the admin user if it doesn't already exist yet (won't clobber its
# password on every restart).
python manage.py shell -c "
from django.contrib.auth.models import User
if not User.objects.filter(username='${NGPHYLO_ADMIN_USER:-admin}').exists():
    User.objects.create_superuser('${NGPHYLO_ADMIN_USER:-admin}', '${NGPHYLO_ADMIN_EMAIL:-admin@example.org}', '${NGPHYLO_ADMIN_PASSWORD:-password}')
"
python manage.py loaddata tool_flags

# Link to a Galaxy server and import its tools/workflows - only if one was
# configured (NGPHYLO_GALAXY_URL/NGPHYLO_GALAXY_KEY - see README.md). Safe to
# skip: the app itself starts fine without a Galaxy server, workflow
# submission just won't work until one is linked (this can also be done
# later - see the one-time management commands in CLAUDE.md).
if [ -n "${NGPHYLO_GALAXY_URL:-}" ] && [ -n "${NGPHYLO_GALAXY_KEY:-}" ]; then
    echo "Linking Galaxy server at $NGPHYLO_GALAXY_URL ..."
    python manage.py creategalaxyserver --url="$NGPHYLO_GALAXY_URL" --activate
    python manage.py addgalaxykey --user "${NGPHYLO_ADMIN_USER:-admin}" --galaxyurl "$NGPHYLO_GALAXY_URL" --galaxykey "$NGPHYLO_GALAXY_KEY"
    python manage.py importtools --galaxyurl="$NGPHYLO_GALAXY_URL" --query="phylogeny" --flags=toolflags.txt --force --inputfields=toolfields.txt
    python manage.py import_links --linkfile=toollinks.txt
    python manage.py importworkflows --galaxyurl="$NGPHYLO_GALAXY_URL" --wfnamefile=wfnames.txt
else
    echo "NGPHYLO_GALAXY_URL/NGPHYLO_GALAXY_KEY not set - skipping Galaxy server setup (see README.md)."
fi

echo "Init done."
