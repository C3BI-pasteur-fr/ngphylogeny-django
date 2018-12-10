from django import forms
from django.conf import settings
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit

from .models import BlastRun

class BlastForm(forms.Form):
    sequence = forms.CharField(label='Sequence', widget=forms.Textarea)

    init_servers = BlastRun.blast_servers().items()
    init_progs =  BlastRun.blast_progs(init_servers[0][0]).items()
    init_dbs = BlastRun.blast_dbs(init_servers[0][0],init_progs[0][0]).items()

    server = forms.ChoiceField(
        label='Server',
        choices=init_servers,
        initial=init_servers[0],
    )

    program = forms.ChoiceField(
        label='Program',
        choices=init_progs,
        initial=init_progs[0],
    )

    database = forms.ChoiceField(
        label='Database',
        choices=init_dbs,
        initial=init_dbs[0],
    )
    
    evalue = forms.FloatField(label='E.Value Threshold', min_value=0, max_value=1, initial=1.0E-5)
    coverage = forms.FloatField(label='Query coverage Threshold', min_value=0, max_value=1, initial=0.8)
    email = forms.CharField(label='Contact Email (optional)', required=False)

    def __init__(self, *args, **kwargs):
        self.helper = FormHelper()
        self.helper.add_input(Submit('submit', 'Submit'))
        super(BlastForm, self).__init__(*args, **kwargs)
