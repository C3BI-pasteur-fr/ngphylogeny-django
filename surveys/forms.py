import re

from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Submit
from django.core.exceptions import ValidationError
from django.forms import ModelForm
from captcha.fields import CaptchaField

from utils.antispam import AntiSpamFormMixin
from .models import Feedback

MAX_LINKS = 2
LINK_RE = re.compile(r'https?://|www\.', re.IGNORECASE)


class FeedbackForm(AntiSpamFormMixin, ModelForm):
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
        # appending another Field('captcha', ...) on top of that rendered
        # the captcha widget twice, so the layout is built explicitly.
        self.helper.layout = Layout(
            'type', 'title', 'comment', 'email',
            *AntiSpamFormMixin.antispam_layout_fields(),
            Field('captcha', placeholder="Enter captcha"),
            FormActions(
                Submit('save', 'Send message'),
            ),
        )

    def clean_comment(self):
        comment = self.cleaned_data.get('comment', '')
        if len(LINK_RE.findall(comment)) > MAX_LINKS:
            raise ValidationError(
                "Please include at most %d links in your message." % MAX_LINKS)
        return comment
