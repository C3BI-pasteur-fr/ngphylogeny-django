"""
Shared anti-spam checks for public forms (contact form, account sign-up) -
see CLAUDE.md's "Contact form anti-spam layers". Each form mixes in
AntiSpamFormMixin and decides for itself what a flagged submission does:
the contact form drops it silently, the account form shows a visible
error (a real person who got flagged there would otherwise be left
believing their sign-up succeeded).
"""
import time

from crispy_forms.layout import Field
from django import forms
from django.core import signing
from django.core.exceptions import ValidationError

FORM_STARTED_SALT = 'utils.antispam.form_started'
# A real person can't read and fill in a form in under a few seconds;
# bots that submit as soon as the page loads can.
MIN_SUBMIT_SECONDS = 3
# Anything older than this has to be reloaded - stops a token captured
# long ago being replayed indefinitely.
MAX_FORM_AGE_SECONDS = 60 * 60

HONEYPOT_FIELD = 'hp_check'
FORM_STARTED_FIELD = 'form_started'
# Hidden from real visitors with CSS (see .antispam-hp in custom.css), not
# display:none - some bots skip display:none inputs. Deliberately not named
# after a standard autofill token like "website": browsers fill those with
# the visitor's own saved details, which would trip the honeypot for a
# real person.
HONEYPOT_CSS_CLASS = 'antispam-hp'


class AntiSpamFormMixin:
    """
    Adds the honeypot and signed form-age fields to a form, plus the
    checks for them. Must come before the form's own base class in the
    MRO so its __init__ runs first and its clean_* hooks are found.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields[HONEYPOT_FIELD] = forms.CharField(
            required=False, label='',
            widget=forms.TextInput(
                attrs={'autocomplete': 'off', 'tabindex': '-1'}))
        self.fields[FORM_STARTED_FIELD] = forms.CharField(
            required=False, widget=forms.HiddenInput)
        if not self.is_bound:
            self.initial[FORM_STARTED_FIELD] = signing.dumps(
                time.time(), salt=FORM_STARTED_SALT)

    @staticmethod
    def antispam_layout_fields():
        return [Field(HONEYPOT_FIELD, wrapper_class=HONEYPOT_CSS_CLASS),
                Field(FORM_STARTED_FIELD)]

    def clean_form_started(self):
        token = self.cleaned_data.get(FORM_STARTED_FIELD, '')
        try:
            self._started_at = signing.loads(
                token, salt=FORM_STARTED_SALT, max_age=MAX_FORM_AGE_SECONDS)
        except signing.BadSignature:
            raise ValidationError(
                "This form has expired - please reload the page and try again.")
        return token

    def looks_like_spam(self):
        """
        Only meaningful once the form is otherwise valid (captcha included),
        so it only ever runs for a submission that already got past the
        real validation.
        """
        if self.cleaned_data.get(HONEYPOT_FIELD):
            return True
        started_at = getattr(self, '_started_at', None)
        return (started_at is not None and
                time.time() - started_at < MIN_SUBMIT_SECONDS)
