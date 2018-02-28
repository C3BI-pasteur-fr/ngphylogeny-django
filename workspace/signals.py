from django.db.models.signals import pre_delete
from django.dispatch import receiver
from .models import WorkspaceHistory


@receiver(pre_delete, sender=WorkspaceHistory)
def send_delete_galaxy_history(sender, instance, using, **kwargs):
    """remove history from db and from Galaxy server"""

    gu = instance.get_galaxy_user()
    if gu:
        gi = gu.get_galaxy_instance
        gi.histories.delete_history(instance.history, purge=True)