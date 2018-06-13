# NGPhylogeny.fr
# https://ngphylogeny.fr

# base image: python 2.7.14-jessie
FROM python:2.7.14-jessie

# File Author / Maintainer
MAINTAINER Frederic Lemoine <frederic.lemoine@pasteur.fr>

COPY . /home/ngphylo

# Install REDIS-SERVER
RUN wget http://download.redis.io/redis-stable.tar.gz \
    && tar xvzf redis-stable.tar.gz \
    && cd redis-stable \
    && make \
    && cp src/redis-server /usr/local/bin/ \
    && cp src/redis-cli /usr/local/bin/ \
    && mkdir /etc/redis \
    && mkdir /home/ngphylo/redis \
    && cp utils/redis_init_script /etc/init.d/redis_6379 \
    && mv /home/ngphylo/docker/redis.conf /etc/redis/6379.conf \
    && mkdir /home/ngphylo/redis/6379 \
    && cd .. && rm -rf redis-stable*

RUN pip install -r /home/ngphylo/requirement.txt

RUN touch /var/log/redis_6379.log
RUN /etc/init.d/redis_6379 start
WORKDIR /home/ngphylo/
RUN python manage.py makemigrations
RUN python manage.py migrate --run-syncdb

COPY docker/celeryd /etc/init.d/celeryd
COPY docker/celerybeat /etc/init.d/celerybeat
COPY docker/celeryd.default /etc/default/celeryd
RUN chmod 640  /etc/default/celeryd
RUN chmod +x /etc/init.d/celery*
RUN mkdir /var/run/celery

# INSTALL NGINX
COPY docker/nginx /etc/init.d/nginx
COPY docker/nginx.default /etc/default/nginx
COPY docker/ngphylogeny_nginx.conf /etc/nginx/nginx.conf
COPY docker/ngphylogeny_nginx.conf /etc/nginx/nginx.conf
COPY docker/uwsgi.init /etc/init/uwsgi.conf

RUN wget http://nginx.org/download/nginx-1.15.0.tar.gz \
    && tar -xzvf nginx-1.15.0.tar.gz \
    && cd nginx-1.15.0 \
    && ./configure \
    --prefix=/usr/share/nginx \
    --sbin-path=/usr/sbin/nginx \
    --conf-path=/etc/nginx/nginx.conf \
    --pid-path=/var/run/nginx.pid \
    --lock-path=/var/lock/nginx.lock \
    --error-log-path=/var/log/nginx/error.log \
    --http-log-path=/var/log/access.log \
    --user=root \
    --group=root \
    --without-mail_pop3_module \
    --without-mail_imap_module \
    --without-mail_smtp_module \
    --without-http_scgi_module \
    --without-http_memcached_module \
    --with-ipv6 \
    --with-http_ssl_module \
    --with-http_stub_status_module \
    --with-http_gzip_static_module \
    && make && make install \
    && chmod +x /etc/init.d/nginx \
    && chmod 640  /etc/default/nginx

ENTRYPOINT ["/home/ngphylo/startup.sh"]
