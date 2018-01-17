import logging

from django.contrib.sessions.backends.db import SessionStore
from django.contrib.sessions.models import Session
from django.db.models.signals import post_save
from django.dispatch import receiver

from .models import Server

logger = logging.getLogger(__name__)


@receiver(post_save, sender=Server)
def update_active_server_in_session(sender, instance, dispatch_uid='teet', **kwargs):
    sessions = Session.objects.all()
    # sessions.delele()
    # logger.info("Sessions clear")

    logged_in = [s.session_key for s in sessions if s.get_decoded().get('galaxy_server')]

    for session_key in logged_in:
        s = SessionStore(session_key=session_key)
        s['galaxy_server'] = instance.id
        print s
        s.save()

    print [s.get_decoded().get('galaxy_server') for s in sessions]
