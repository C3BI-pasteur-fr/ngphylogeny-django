# -*- coding: utf-8 -*-
from crispy_forms.bootstrap import FormActions
from crispy_forms.bootstrap import StrictButton
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout
from django import forms


class UploadForm(forms.Form):
    input_file = forms.FileField(widget=forms.FileInput(attrs={'multiple': 'true'}))

    def __init__(self, *args, **kwargs):
        super(UploadForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.form_method = 'POST'
        self.helper.layout = Layout('input_file',
                                    FormActions(
                                        Submit('submit', 'Submit')))

    class Media:

        js = ("js/jquery.filer.min.js",)
        css= {
                'all':("css/jquery.filer.css",)
              }


class PastedContentForm(forms.Form):
    pasted_text = forms.CharField(widget=forms.Textarea)
    textarea_id = "id_pasted_text"
    helper = FormHelper()
    helper.form_method = 'POST'
    helper.layout = Layout('pasted_text',
                           FormActions(
                                Submit('submit', 'Submit'),
                                StrictButton( "<span class='glyphicon glyphicon-question-sign'></span> Example",
                                         css_class="btn btn-info",
                                         onclick="populateExample('"+textarea_id+"', FASTA_AA);"),
                                StrictButton("<span class='glyphicon glyphicon-trash'></span>",
                                        css_class="btn btn-default",
                                        onclick="document.getElementById('"+textarea_id+"').value='';" ),
                                ))
    class Media:
        js = ('js/form_examples.js',)