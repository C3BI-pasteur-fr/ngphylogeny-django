from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Submit
from django.forms import ModelForm
from captcha.fields import CaptchaField

from .models import Feedback


class FeedbackForm(ModelForm):
    """Model Feedback form"""
    captcha = CaptchaField()
    class Meta:
        model = Feedback
        fields = ['type', 'title', 'comment', 'email']

    def __init__(self, *args, **kwargs):
        super(FeedbackForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper(self)
        # FormHelper(self) already builds a default layout with every
        # form field (including captcha) via build_default_layout() -
        # appending another Field('captcha', ...) on top of that (the
        # previous code) rendered the captcha widget twice. A typo
        # ('captcha ', trailing space) had been silently masking that
        # second bug: crispy_forms's own FAIL_SILENTLY handling logs a
        # "Could not resolve form field" warning and renders nothing for
        # an unresolvable field name instead of raising, so the page
        # never actually broke - it just silently dropped the intended
        # placeholder customization and, since the captcha field was
        # never resolved, never submitted a real captcha value either,
        # making every contact-form submission fail captcha validation.
        # Building the layout explicitly (rather than appending to the
        # implicit default) fixes both: one captcha widget, with the
        # intended placeholder.
        self.helper.layout = Layout(
            'type', 'title', 'comment', 'email',
            Field('captcha', placeholder="Enter captcha"),
            FormActions(
                Submit('save', 'Send message'),
            ),
        )
