# NGPhylogeny.fr
# https://ngphylogeny.fr

# base image: python 2.7.14-jessie
FROM python:2.7.14-jessie

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

RUN useradd -m ngphylo
COPY . /home/ngphylo

RUN pip install -r /home/ngphylo/requirement.txt

RUN chown -R ngphylo:ngphylo /home/ngphylo \
    && chmod 755 /home/ngphylo/startup.sh
USER ngphylo
WORKDIR /home/ngphylo/
RUN python manage.py makemigrations
RUN python manage.py migrate --run-syncdb

ENTRYPOINT ["/home/ngphylo/startup.sh"]
