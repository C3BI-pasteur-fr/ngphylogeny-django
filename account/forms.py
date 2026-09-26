from captcha.fields import CaptchaField
from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Field, Layout, Submit
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User


class AccountCreationForm(UserCreationForm):
    """
    Public account sign-up. Builds on django.contrib.auth's own
    UserCreationForm (username uniqueness, password confirmation, and
    AUTH_PASSWORD_VALIDATORS - see settings/base.py - all for free)
    rather than hand-rolling that validation, plus an email field (not
    part of the base form) and a captcha - matching surveys.forms.
    FeedbackForm's own established pattern for this codebase's public-
    facing forms, and if anything more warranted here: an account
    sign-up form is a more obvious bot/spam target than a contact form.
    """
    email = forms.EmailField(required=True)
    captcha = CaptchaField()

    class Meta(UserCreationForm.Meta):
        model = User
        fields = ('username', 'email')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            'username', 'email', 'password1', 'password2',
            Field('captcha', placeholder="Enter captcha"),
            FormActions(
                Submit('save', 'Create account'),
            ),
        )

    def clean_email(self):
        email = self.cleaned_data['email']
        if User.objects.filter(email__iexact=email).exists():
            raise forms.ValidationError(
                "An account with this email already exists.")
        return email
