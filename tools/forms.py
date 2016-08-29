from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout, Field, Div
from django import forms


def map_galaxy_tool_input(attr):
        """take Galaxy Tool build information return django field dict"""

        field_map = dict()
        field_map['initial'] = attr.get('default_value', attr.get('value', ""))
        field_map['required'] = attr.get('optional', "")
        field_map['help_text'] = attr.get('help', "")
        field_map['label'] = attr.get('label', "")

        if attr.get("type", "") == "select":
            field_map['choices'] = []
            for opt in attr.get("options", ""):
                field_map['choices'].append((opt[1], opt[0]))

        if attr.get("type", "") == "data":

            opt = attr.get("options", "").get('hda')
            if opt:
                field_map['choices'] = []
            for data in opt:
                field_map['choices'].append((data.get('id'), data.get('name')))
                print field_map['choices']

        return field_map



class ToolForm(forms.Form):

    n = 0
    fieds_ids_mapping = {}
    tool_params = None

    def create_field(self, attrfield):

        self.n +=1
        field_id = str(self.n)
        fieldtype = attrfield.get("type", "")

        if fieldtype == "data":
            choices = attrfield.get("options", "").get('hda')
            self.input_file_ids.append(field_id)

            if choices:
                self.fields[field_id] = forms.ChoiceField(**map_galaxy_tool_input(attrfield))
            else:

                self.fields[field_id] = forms.FileField(**map_galaxy_tool_input(attrfield))

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

    def parse_galaxy_input_tool(self, list_inputs, cond_name=''):

        fields_created = []

        for input_tool in list_inputs:
            cond_input = input_tool.get('test_param')

            if cond_input:
                cond_name = input_tool.get('name') + '|'

                conditional_field = self.create_field(cond_input)
                self.fieds_ids_mapping[conditional_field] = cond_name+cond_input['name']

                nested_field = []
                cases = input_tool.get('cases', '')
                for case in cases:
                    case_inputs = case.get('inputs')
                    if case_inputs:
                        case_value = case.get('value')
                        nested_field.append(Div(data_test=conditional_field, data_case=case_value, css_class="well",
                                                *self.parse_galaxy_input_tool(case_inputs, cond_name)))

                fields_created.append(Div(conditional_field, *nested_field))
                cond_name = ''

            else:
                new_field = self.create_field(input_tool)
                fields_created.append(Field(new_field))
                self.fieds_ids_mapping[new_field] = cond_name + input_tool['name']

        return fields_created

    def __init__(self, tool_params=None, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)

        if tool_params:
            self.tool_params = tool_params
        assert self.tool_params, 'tool_params is needed'

        self.input_file_ids = []
        self.helper = FormHelper(self)
        self.helper.form_class = 'blueForms'
        self.helper.form_tag = False
        self.formset = self.parse_galaxy_input_tool(self.tool_params)
        self.helper.layout = Layout(FormActions(Field(*self.formset),
                                                Submit('submit', 'Submit', css_class="pull-right"),
                                                )
                                    )

