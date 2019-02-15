# -*- coding: utf-8 -*-
from crispy_forms.bootstrap import FormActions
from crispy_forms.bootstrap import StrictButton
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout
from django import forms

from Bio import SeqIO
import StringIO

from blast.models import BlastRun
from utils import biofile

class UploadForm(forms.Form):
    input_file = forms.FileField(widget=forms.FileInput(), required=False)
    pasted_text = forms.CharField(widget=forms.Textarea, required=False)
    blast_run = forms.ChoiceField()
    
    def __init__(self, blastruns=None, *args, **kwargs):
        super(UploadForm, self).__init__(*args, **kwargs)

        # We fill the select options with id of blast runs
        choices= []
        choices.append(
            ('--','--'),
        )
        if blastruns is not None:
            for (bid, bname) in blastruns:
                c = (bid,"%s (%s)" % (str(bname), str(bid)))
                choices.append(c)
        self.fields['blast_run'] = forms.ChoiceField(choices=choices)
        
        textarea_id = "id_pasted_text"

        self.helper = FormHelper()
        self.helper.form_method = 'POST'
        self.helper.layout = Layout('input_file','pasted_text','blast_run',
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
        js = ("js/jquery.filer.min.js",'js/form_examples.js',)
        css= {
                'all':("css/jquery.filer.css",)
              }
