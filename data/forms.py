# -*- coding: utf-8 -*-
from crispy_forms.bootstrap import FormActions
from crispy_forms.bootstrap import StrictButton
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout
from django import forms

from Bio import SeqIO
import StringIO

from blast.models import BlastRun


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
            for b in blastruns:
                obj = BlastRun.objects.get(pk=b)
                c = (b,"%s (%s)" % (obj.query_id, str(obj.id)))
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


            
    def validate_form_inputs(self):
        valid = True
        # Check uploaded file or pasted content
        submitted_file = self.cleaned_data.get('input_file')
        submitted_text = self.cleaned_data.get('pasted_text')
        # Id of a blastrun : We will generate a fasta file from
        # blast results
        blast_run    = self.cleaned_data.get('blast_run')
        if submitted_file:
            nbseq = 0
            #print submitted_file
            for r in SeqIO.parse(submitted_file, "fasta"):
                nbseq += 1
            if nbseq == 0:
                self.add_error(
                    'input_file', "Input file format is not FASTA or file is empty")
                valid = False
            elif nbseq <= 3:
                self.add_error(
                    'input_file',"Input file should contain more than 3 sequences")
                valid = False
        elif submitted_text:
            nbseq = 0
            for r in SeqIO.parse(StringIO.StringIO(submitted_text), "fasta"):
                nbseq += 1
            if nbseq == 0:
                self.add_error(
                    'pasted_text', "Input data format is not FASTA or file is empty")
                valid = False
            elif nbseq <= 3:
                self.add_error(
                    'pasted_text',"Input data should contain more than 3 sequences")
                valid = False
        elif blast_run != '--':
            nbseq = 0
            blastrun = BlastRun.objects.get(pk=blast_run)
            
            for r in SeqIO.parse(StringIO.StringIO(blastrun.to_fasta()), "fasta"):
                nbseq += 1
            if nbseq == 0:
                self.add_error(
                    'pasted_text', "BlastRun data format is not FASTA or file is empty")
                valid = False
            elif nbseq <= 3:
                self.add_error(
                    'pasted_text',"BlastRun data should contain more than 3 sequences")
                valid = False
        else:
            self.add_error(
                None, "No input data has been provided")
            valid = False
        return valid

    class Media:
        js = ("js/jquery.filer.min.js",'js/form_examples.js',)
        css= {
                'all':("css/jquery.filer.css",)
              }
