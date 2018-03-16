from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field, Div, Fieldset
from django import forms
from .models import ToolFieldWhiteList


def map_galaxy_tool_input(attr):
    """Convert Galaxy tool input information to django field attribute dict"""

    field_map = {'initial': attr.get('default_value', attr.get('value', "")),
                 'required': attr.get('optional', ""),
                 'help_text': attr.get('help', ""),
                 'label': attr.get('label', ""),
                 }

    if attr.get("type", "") == "select":
        field_map['choices'] = []
        for opt in attr.get("options", ""):
            field_map['choices'].append((opt[1], opt[0]))

    if attr.get("type", "") == "data":

        opt = attr.get("options", "").get('hda')
        if opt:
            field_map['choices'] = []
        # add dataset of history in a list
        for data in opt:
            field_map['choices'].append((data.get('id'), data.get('name')))

    return field_map


class ToolForm(forms.Form):
    """"""
    tool_id = None
    tool_params = []
    # map field id to galaxy params name
    fields_ids_mapping = {}
    visible_field = []
    n = 0

    def create_field(self, attrfield):
        self.n += 1
        field_id = str(self.n)
        fieldtype = attrfield.get("type", "")
        if fieldtype == "data":
            choices = attrfield.get("options", "").get('hda')
            self.input_file_ids.append(field_id)
            if choices:
                self.fields[field_id] = forms.ChoiceField(**map_galaxy_tool_input(attrfield))
            else:
                self.fields[field_id] = forms.FileField(**map_galaxy_tool_input(attrfield))
                self.fields[field_id].widget.attrs = ({'data-ext': attrfield.get('extensions')})
                self.fields[field_id].widget.attrs = ({'data-name': attrfield.get('name')})

        elif fieldtype == "select":
            if attrfield.get("display", "") == 'radio':
                self.fields[field_id] = forms.ChoiceField(widget=forms.RadioSelect, **map_galaxy_tool_input(attrfield))
            else:
                self.fields[field_id] = forms.ChoiceField(**map_galaxy_tool_input(attrfield))

        elif fieldtype == "boolean":
            data_onvalue = attrfield.get("truevalue", None)
            data_offvalue = attrfield.get("falsevalue", None)

            self.fields[field_id] = forms.CharField(widget=forms.CheckboxInput(
                attrs={'data-toggle': "toggle",
                       'data-on': 'Yes',
                       'data-off': 'No',
                       'data-on-value': data_onvalue,
                       'data-off-value': data_offvalue,
                       }),
                **map_galaxy_tool_input(attrfield))

        elif fieldtype == "text":
            self.fields[field_id] = forms.CharField(**map_galaxy_tool_input(attrfield))

        elif fieldtype == "integer":
            self.fields[field_id] = forms.IntegerField(**map_galaxy_tool_input(attrfield))

        else:
            self.fields[field_id] = forms.CharField(**map_galaxy_tool_input(attrfield))

        return field_id

    def parse_galaxy_input_tool(self, list_inputs, cond_name='', whitelist=list()):

        fields_created = []

        for input_tool in list_inputs:
            if input_tool.get('type') == 'section':
                fields_created.append(Fieldset(input_tool.get('title'),
                                               *self.parse_galaxy_input_tool(input_tool.get('inputs'))))

            elif input_tool.get('name') in whitelist or not whitelist:

                cond_input = input_tool.get('test_param')
                if cond_input:

                    cond_name = input_tool.get('name') + '|'
                    conditional_field = self.create_field(cond_input)
                    self.fields_ids_mapping[conditional_field] = cond_name + cond_input['name']

                    nested_field = []
                    cases = input_tool.get('cases', '')
                    for case in cases:
                        case_inputs = case.get('inputs')
                        if case_inputs:
                            case_value = case.get('value')
                            if case_value is not None:
                                data_test = conditional_field
                                if self.prefix:
                                    data_test = self.prefix + '-' + data_test

                                nested_field.append(Div(data_test=data_test,
                                                        data_case=case_value,
                                                        css_class="well",
                                        *self.parse_galaxy_input_tool(case_inputs, cond_name)))

                    fields_created.append(Div(conditional_field, *nested_field))
                    cond_name = ''

                else:
                    new_field = self.create_field(input_tool)
                    fields_created.append(Field(new_field))
                    self.fields_ids_mapping[new_field] = cond_name + input_tool['name']

        return fields_created

    def __init__(self, tool_params=None, tool_id=None, whitelist=None, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)

        self.tool_params = tool_params or self.tool_params
        self.visible_field = whitelist or self.visible_field
        self.tool_id = tool_id or self.tool_id
        self.input_file_ids = []
        self.helper = FormHelper(self)
        self.helper.form_class = 'blueForms'
        self.helper.form_tag = False
        self.formset = self.parse_galaxy_input_tool(self.tool_params, whitelist=self.visible_field)
        self.helper.layout = Layout(FormActions(Field(*self.formset)))

class ToolFieldWhiteListForm(forms.ModelForm):
    _params = forms.MultipleChoiceField(label='Params')

    model = ToolFieldWhiteList
    fields = ['tool', 'context', '_params', ]

    def clean(self):
        cleaned_data = super(ToolFieldWhiteListForm, self).clean()
        cleaned_data['_params'] = ",".join(self.cleaned_data.get('_params', []))
        return cleaned_data
