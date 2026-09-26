from django.conf import settings
from django.db import models


class UserProfile(models.Model):
    """
    Minimal per-account record, one per Django User - not tied to a
    specific Galaxy server the way GalaxyUser is (galaxy/models.py).
    Exists as the natural place to grow account-level attributes without
    overloading auth.User directly - for now just a creation timestamp.

    get_or_create()'d lazily wherever it's needed (see
    galaxy.views.UpdateApiKey's own GalaxyUser precedent for the same
    pattern), not created via a post_save signal - so it needs no
    backfill migration for accounts that already existed before this
    model did.
    """
    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE,
        related_name='profile')
    created_date = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return str(self.user)
