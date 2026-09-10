# NGPhylogeny.fr
# https://ngphylogeny.fr
#
# One image, reused for the web (runserver), celery-worker, and celery-beat
# services in docker-compose.yml (each just overrides the command). Postgres
# and Redis run as their own separate compose services rather than being
# built into this image.

FROM python:3.8-buster

LABEL maintainer="Frederic Lemoine <frederic.lemoine@pasteur.fr>"

WORKDIR /home/ngphylo

# Debian buster is EOL: deb.debian.org no longer mirrors it. Point apt at
# the snapshot.debian.org archive already listed (commented out) in this
# base image's sources.list, and drop the security repo line - its snapshot
# is old enough to report as expired and fail `apt-get update` outright.
RUN sed -i \
        -e 's|^deb http://deb\.debian\.org|# deb http://deb.debian.org|' \
        -e 's|^# deb http://snapshot\.debian\.org|deb http://snapshot.debian.org|' \
        -e '/debian-security/d' \
        /etc/apt/sources.list \
    && apt-get update --fix-missing \
    && apt-get install -y libpq-dev postgresql-client libmagic1 \
    && rm -rf /var/lib/apt/lists/*

COPY requirement.txt .
# numpy must land before biopython, which needs it at build time and doesn't
# declare it - a plain `pip install -r requirement.txt` can fail on a clean
# env depending on pip's resolution order.
RUN pip install numpy==1.24.4 \
    && pip install -r requirement.txt

COPY . .

RUN chmod +x docker/init.sh

CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
