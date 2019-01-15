# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.urls import reverse_lazy
from django.views.generic import CreateView, TemplateView

from surveys.forms import FeedbackForm
from surveys.models import Feedback


class FeedbackCreateView(CreateView):
    """

    """
    model = Feedback
    form_class = FeedbackForm
    success_url = reverse_lazy('feedback_success')


class FeedbackSuccessView(TemplateView):

    template_name = "surveys/thankfeedback.html"
