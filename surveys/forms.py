import re

from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import HTML, Layout, Field, Submit
from django.core.exceptions import ValidationError
from django.forms import ModelForm

from utils.antispam import AntiSpamFormMixin
from utils.powcaptcha import PowCaptchaFormMixin, TOKEN_FIELD, NONCE_FIELD
from .models import Feedback

MAX_LINKS = 2
LINK_RE = re.compile(r'https?://|www\.', re.IGNORECASE)


class FeedbackForm(AntiSpamFormMixin, PowCaptchaFormMixin, ModelForm):
    """Model Feedback form"""

    class Meta:
        model = Feedback
        fields = ['type', 'title', 'comment', 'email']

    def __init__(self, *args, **kwargs):
        super(FeedbackForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper(self)
        # FormHelper(self) already builds a default layout with every
        # form field via build_default_layout() - appending another
        # Field() for one of those fields on top of that rendered it
        # twice (see this form's git history for the real bug that hit
        # exactly this with the old captcha field), so the layout is
        # built explicitly instead.
        self.helper.layout = Layout(
            'type', 'title', 'comment', 'email',
            *AntiSpamFormMixin.antispam_layout_fields(),
            Field(TOKEN_FIELD), Field(NONCE_FIELD),
            HTML(self.pow_layout_html()),
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
