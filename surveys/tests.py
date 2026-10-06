# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import time
from smtplib import SMTPException
from unittest.mock import Mock, patch

from django.core import signing
from django.test import TestCase, override_settings
from django.urls import reverse

from surveys.emails import (
    build_feedback_notification_email, send_feedback_notification_email)
from surveys.forms import MAX_LINKS, FeedbackForm
from surveys.models import Feedback
from surveys.views import FeedbackCreateView
from utils.antispam import FORM_STARTED_SALT, HONEYPOT_FIELD, MAX_FORM_AGE_SECONDS


def _feedback(**kwargs):
    kwargs.setdefault('email', 'user@example.org')
    kwargs.setdefault('title', 'Something is broken')
    kwargs.setdefault('comment', 'It broke when I clicked the button.')
    kwargs.setdefault('type', 'bug')
    return Feedback(**kwargs)


class FeedbackNotificationEmailTest(TestCase):
    """
    Tests for surveys/emails.py - the contact-form notification sent to
    staff (NGPHYLO_CONTACT_FORM_RECIPIENTS) in addition to (not instead
    of) saving the submitted Feedback row.
    """

    @override_settings(
        NGPHYLO_CONTACT_FORM_RECIPIENTS=['staff@example.org'],
        NGPHYLO_REPORT_FROM_EMAIL='ngphylogeny@pasteur.fr')
    def test_built_email_content_and_recipients(self):
        feedback = _feedback()
        msg = build_feedback_notification_email(feedback)

        self.assertEqual(msg.to, ['staff@example.org'])
        self.assertEqual(msg.from_email, 'ngphylogeny@pasteur.fr')
        self.assertEqual(msg.reply_to, ['user@example.org'])
        self.assertIn('Something is broken', msg.subject)
        self.assertIn('user@example.org', msg.body)
        self.assertIn('It broke when I clicked the button.', msg.body)
        self.assertIn('Bug', msg.body)

    @override_settings(
        NGPHYLO_CONTACT_FORM_RECIPIENTS=[
            'staff1@example.org', 'staff2@example.org'])
    def test_multiple_comma_separated_recipients(self):
        msg = build_feedback_notification_email(_feedback())
        self.assertEqual(
            msg.to, ['staff1@example.org', 'staff2@example.org'])

    @override_settings(NGPHYLO_CONTACT_FORM_RECIPIENTS=[])
    def test_send_no_ops_when_no_recipients_configured(self):
        # Same "safe to leave enabled everywhere" convention as
        # workspace.tasks.send_daily_report/NGPHYLO_REPORT_RECIPIENTS -
        # left unset, sending is skipped entirely rather than raising or
        # sending to nobody.
        with patch('surveys.emails.EmailMessage') as mock_email:
            send_feedback_notification_email(_feedback())
        mock_email.assert_not_called()

    @override_settings(NGPHYLO_CONTACT_FORM_RECIPIENTS=['staff@example.org'])
    def test_send_actually_sends_when_recipients_configured(self):
        with patch(
                'surveys.emails.EmailMessage.send',
                return_value=1) as mock_send:
            send_feedback_notification_email(_feedback())
        mock_send.assert_called_once_with(fail_silently=False)


class FeedbackCreateViewTest(TestCase):
    """
    FeedbackCreateView.form_valid() must save the Feedback row (the
    view's original, pre-existing behavior) *and* send the staff
    notification email - and a failure sending that email (e.g. no SMTP
    configured) must not turn an already-successful submission into a
    500 page, the same way workspace/blast's own completion-email sends
    already swallow SMTP failures from their Celery tasks (see
    workspace/tasks.py).
    """

    def _saving_form(self):
        # A stand-in for a real, captcha-validated FeedbackForm: only
        # form_valid()'s own behavior (save + notify) is under test
        # here, not django-simple-captcha's own validation, which
        # data.tests.StaticPagesSmokeTest already exercises via a real
        # GET render of this same form. ModelFormMixin.form_valid()
        # only ever calls form.save(), so a bare Mock whose save()
        # returns an actually-saved Feedback row is enough.
        saved = Feedback.objects.create(
            type='bug', title='Something is broken',
            comment='It broke when I clicked the button.',
            email='user@example.org')
        form = Mock()
        form.save.return_value = saved
        form.looks_like_spam.return_value = False
        return form, saved

    def test_form_valid_saves_row_and_sends_notification(self):
        view = FeedbackCreateView()
        view.request = self.client.get('/').wsgi_request
        form, saved = self._saving_form()
        with patch(
                'surveys.views.send_feedback_notification_email'
                ) as mock_send:
            view.form_valid(form)

        self.assertEqual(Feedback.objects.count(), 1)
        mock_send.assert_called_once_with(saved)

    def test_form_valid_swallows_smtp_failure(self):
        view = FeedbackCreateView()
        view.request = self.client.get('/').wsgi_request
        form, saved = self._saving_form()
        with patch(
                'surveys.views.send_feedback_notification_email',
                side_effect=SMTPException('boom')):
            # Must not raise - the row is already saved by this point.
            view.form_valid(form)

        self.assertEqual(Feedback.objects.count(), 1)


class FeedbackAntiSpamTest(TestCase):
    """
    Anti-spam layers on FeedbackForm (see surveys/forms.py): a honeypot
    field, a signed form-age check that rejects instant submissions, and
    a cap on links in the message. Form-level checks run the real
    full_clean() - the captcha is intentionally left unanswered here, so
    assertions only look at the non-captcha fields' errors and
    looks_like_spam(), not is_valid() as a whole.
    """

    @staticmethod
    def _data(**overrides):
        data = {
            'type': 'bug',
            'title': 'Something is broken',
            'comment': 'It broke when I clicked the button.',
            'email': 'user@example.org',
            HONEYPOT_FIELD: '',
            'form_started': signing.dumps(
                time.time() - 60, salt=FORM_STARTED_SALT),
        }
        data.update(overrides)
        return data

    def test_normal_submission_is_not_spam(self):
        form = FeedbackForm(data=self._data())
        form.is_valid()
        self.assertNotIn('form_started', form.errors)
        self.assertFalse(form.looks_like_spam())

    def test_filled_honeypot_is_spam(self):
        form = FeedbackForm(data=self._data(**{HONEYPOT_FIELD: 'http://spam.example'}))
        form.is_valid()
        self.assertTrue(form.looks_like_spam())

    def test_too_fast_submission_is_spam(self):
        instant = signing.dumps(time.time(), salt=FORM_STARTED_SALT)
        form = FeedbackForm(data=self._data(form_started=instant))
        form.is_valid()
        self.assertNotIn('form_started', form.errors)
        self.assertTrue(form.looks_like_spam())

    def test_missing_or_tampered_form_started_is_a_validation_error(self):
        for bad in ['', 'not-a-real-token']:
            form = FeedbackForm(data=self._data(form_started=bad))
            form.is_valid()
            self.assertIn('form_started', form.errors)

    def test_expired_form_started_is_a_validation_error(self):
        # signing.dumps() stamps the token with its own creation time, so
        # backdate that clock (not just the payload) to actually expire it.
        with patch('django.core.signing.time.time',
                   return_value=time.time() - MAX_FORM_AGE_SECONDS - 60):
            old = signing.dumps(time.time(), salt=FORM_STARTED_SALT)
        form = FeedbackForm(data=self._data(form_started=old))
        form.is_valid()
        self.assertIn('form_started', form.errors)

    def test_too_many_links_is_a_validation_error(self):
        comment = ' '.join(['https://a.example'] * (MAX_LINKS + 1))
        form = FeedbackForm(data=self._data(comment=comment))
        form.is_valid()
        self.assertIn('comment', form.errors)

    def test_links_up_to_the_cap_are_allowed(self):
        comment = ' '.join(['https://a.example'] * MAX_LINKS)
        form = FeedbackForm(data=self._data(comment=comment))
        form.is_valid()
        self.assertNotIn('comment', form.errors)

    def test_unbound_form_renders_a_fresh_form_started_token(self):
        token = FeedbackForm().initial['form_started']
        started_at = signing.loads(token, salt=FORM_STARTED_SALT)
        self.assertLess(time.time() - started_at, 5)


class FeedbackSpamDropViewTest(TestCase):
    """
    A submission that looks_like_spam() is dropped silently: redirected
    to the thank-you page exactly like a real one, but nothing saved and
    nothing emailed to staff.
    """

    @patch('surveys.views.send_feedback_notification_email')
    def test_spam_is_not_saved_or_emailed_but_still_redirects_as_success(
            self, mock_send):
        form = Mock()
        form.looks_like_spam.return_value = True
        response = FeedbackCreateView().form_valid(form)

        self.assertEqual(Feedback.objects.count(), 0)
        mock_send.assert_not_called()
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, reverse('feedback_success'))
