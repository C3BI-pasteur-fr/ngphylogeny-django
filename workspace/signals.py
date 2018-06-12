from django.db.models.signals import pre_delete
from django.dispatch import receiver
from .models import WorkspaceHistory
from tasks import deletegalaxyhistory

@receiver(pre_delete, sender=WorkspaceHistory, dispatch_uid='history_delete_signal')
def send_delete_galaxy_history(sender, instance, using, **kwargs):
    """
    removes history from db and from Galaxy server
    in an asynchronous call to celery
    """
    print "Removing an history signal"
    deletegalaxyhistory.delay(instance.history)
