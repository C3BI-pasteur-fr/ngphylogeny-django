# -*- coding: utf-8 -*-
from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout
from django import forms


class UploadForm(forms.Form):
    file = forms.FileField(widget=forms.FileInput(attrs={'multiple': 'true'}))
    helper = FormHelper()
    helper.form_method = 'POST'
    helper.layout = Layout('file',
                           FormActions(Submit('submit', 'Submit')))


class PastedContentForm(forms.Form):
    pasted_text = forms.CharField(widget=forms.Textarea)
    helper = FormHelper()
    helper.form_method = 'POST'
    helper.layout = Layout('pasted_text',
                           FormActions(Submit('submit', 'Submit')))
