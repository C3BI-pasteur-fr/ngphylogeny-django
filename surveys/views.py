# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
from smtplib import SMTPException

from django.urls import reverse_lazy
from django.views.generic import CreateView, TemplateView

from surveys.emails import send_feedback_notification_email
from surveys.forms import FeedbackForm
from surveys.models import Feedback


class FeedbackCreateView(CreateView):
    """

    """
    model = Feedback
    form_class = FeedbackForm
    success_url = reverse_lazy('feedback_success')

    def form_valid(self, form):
        # Saves the Feedback row first (super().form_valid() does the
        # save) - the email is in addition to that, never a replacement
        # for it. Unlike workspace/blast's own completion emails (sent
        # from a Celery task, where a failure only affects that async
        # job), this runs inline in the submitter's own request - an
        # uncaught SMTP failure here would turn an already-successful
        # submission into a 500 page, so it's swallowed the same way
        # those tasks already swallow their own send failures (see
        # workspace/tasks.py's updateworkspacestatus/blast/tasks.py).
        response = super().form_valid(form)
        try:
            send_feedback_notification_email(self.object)
        except SMTPException as e:
            logging.warning("Problem with smtp server : %s" % (e))
        except Exception as e:
            logging.warning(
                "Unknown problem while sending feedback notification "
                "e-mail: %s" % (e))
        return response


class FeedbackSuccessView(TemplateView):

    template_name = "surveys/thankfeedback.html"
