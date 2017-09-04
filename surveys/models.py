# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models


class Feedback(models.Model):
    FEEDBACK_TYPES = (
        ('bug', 'Bug'),
        ('feature', 'Feature'),
        ('support', 'Support request'),
        ('other', 'Other')
    )

    email = models.EmailField(null=False, verbose_name='E-mail')
    title = models.CharField(max_length=120, blank=True, verbose_name='Message subject')
    comment = models.TextField(max_length=1200, blank=True, verbose_name='Your message')
    type = models.CharField(max_length=20, blank=True, choices=(FEEDBACK_TYPES), verbose_name='Type of request')
    timestamp = models.DateTimeField(auto_now_add=True, editable=False)
    resolve = models.BooleanField(default=False)
