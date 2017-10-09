# -*- coding: utf-8 -*-
from django.contrib import admin

from models import Feedback


class FeedbackAdmin(admin.ModelAdmin):
    list_display = ('pk','type','title','email','resolve')


admin.site.register(Feedback,FeedbackAdmin)