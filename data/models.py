from django.db import models


class ExampleFile(models.Model):
    """
    Store, manage example data file
    """

    name = models.CharField(max_length=32)
    ext = models.CharField(max_length=8)
    type = models.CharField(max_length=12)
    upload = models.FileField(upload_to='example/')
    comment = models.TextField(max_length=200, null=True, blank=True)
