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

# Debian Buster's ca-certificates package (20200601~deb10u2, frozen since
# Buster went EOL - see the apt sources fix above) predates newer root CAs
# such as HARICA TLS RSA Root CA 2021, which e.g. smtp.pasteur.fr's
# certificate chains through - TLS connections to servers using such a CA
# fail with CERTIFICATE_VERIFY_FAILED. certifi ships Mozilla's current CA
# bundle and gets regular updates on PyPI independent of Buster's own
# frozen apt archive; use it as the system bundle.
#
# Just overwriting /etc/ssl/certs/ca-certificates.crt isn't enough on its
# own: this image's Python was built with no working default `cafile`
# (`ssl.get_default_verify_paths()` reports one that doesn't exist on
# disk), so verification falls back to `capath` - a *directory* of
# individual certs plus OpenSSL hash-named symlinks maintained by
# `update-ca-certificates`, not a single concatenated file - and dropping
# one file there doesn't regenerate those symlinks. SSL_CERT_FILE is the
# one override `ssl.get_default_verify_paths()` explicitly documents
# (`openssl_cafile_env`): set it to force every TLS connection in this
# container to use the up-to-date single-file bundle as `cafile` directly,
# sidestepping the capath/hash-symlink path entirely.
RUN CERTIFI_PATH=$(python -c "import certifi; print(certifi.where())") \
    && cp "$CERTIFI_PATH" /etc/ssl/certs/ca-certificates.crt
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt

COPY . .

RUN chmod +x docker/init.sh

# Baked in at build time, not run by docker/init.sh alone: on Kubernetes,
# the one-shot init Job and the web Deployment are separate pods with no
# shared filesystem by default (unlike docker-compose.yml's
# ngphylo-static volume) - collectstatic here means every container from
# this image already has STATIC_ROOT populated, no shared volume needed.
# Doesn't need real secrets/DB connectivity: collectstatic only reads
# each app's own static/ directory and settings.STATIC_ROOT, and
# NGPHYLO_SECRET_KEY/DATABASES both have safe fallback defaults (see
# settings/base.py) that let Django's settings module load at all.
RUN python manage.py collectstatic --noinput

CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
