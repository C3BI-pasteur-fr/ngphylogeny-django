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
        self.helper.layout.append(
            Field('captcha ', placeholder="Enter captcha")
        )
        self.helper.layout.append(
            FormActions(
                Submit('save', 'Send message'),
            )
        )
