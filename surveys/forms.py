from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit
from django.forms import ModelForm

from models import Feedback


class FeedbackForm(ModelForm):
    """Model Feedback form"""

    class Meta:
        model = Feedback
        fields = ['type', 'title', 'comment', 'email']

    def __init__(self, *args, **kwargs):
        super(FeedbackForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper(self)
        self.helper.layout.append(
        FormActions(
            Submit('save', 'Send message'),
        ))
