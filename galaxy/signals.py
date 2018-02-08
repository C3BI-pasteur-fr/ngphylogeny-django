import logging

from django.contrib.sessions.backends.db import SessionStore
from django.contrib.sessions.models import Session
from django.db.models.signals import post_save
from django.dispatch import receiver

from .models import Server

logger = logging.getLogger(__name__)


@receiver(post_save, sender=Server)
def update_active_server_in_session(sender, instance, dispatch_uid='', **kwargs):
    "Force to update session server storage"
    sessions = Session.objects.all()
    for s in sessions:
        data = SessionStore(session_key=s.session_key).decode(s.session_data)
        data['server_change'] = True
        s.session_data = SessionStore(session_key=s.session_key).encode(data)
        s.save()

