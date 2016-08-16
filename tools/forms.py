from crispy_forms.bootstrap import FormActions
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Submit, Layout, Field, Div
from django import forms


def mapGalaxyToolInput(attr):
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

        return field_map



class ToolForm(forms.Form):

    def create_field(self, attrfield):

        fieldtype = attrfield.get("type", "")
        fieldname = attrfield.get("name", "")

        if fieldtype == "data":
            self.fields[fieldname] = forms.FileField(**mapGalaxyToolInput(attrfield))

        elif fieldtype == "select":
            self.fields[fieldname] = forms.ChoiceField(**mapGalaxyToolInput(attrfield))

        elif fieldtype == "boolean":
            data_onvalue = attrfield.pop("truevalue", None)
            data_offvalue = attrfield.pop("falsevalue", None)

            self.fields[fieldname] = forms.CharField(widget=forms.CheckboxInput(
                                                     attrs={'data-toggle': "toggle",
                                                            'data-on': 'Yes',
                                                            'data-off': 'No',
                                                            'data-on-value': data_onvalue,
                                                            'data-off-value': data_offvalue,
                                                            }),
                                                     **mapGalaxyToolInput(attrfield))

        elif fieldtype == "text":
            self.fields[fieldname] = forms.CharField(**mapGalaxyToolInput(attrfield))

        elif fieldtype == "integer":
            self.fields[fieldname] = forms.IntegerField(**mapGalaxyToolInput(attrfield))

        else:
            self.fields[fieldname] = forms.CharField(**mapGalaxyToolInput(attrfield))

        return fieldname

    def parse_galaxy_input_tool(self, list_inputs):

        fields_created = []

        for input_tool in list_inputs:
            cond_input = input_tool.get('test_param')

            if cond_input:
                nested_field = []
                new_name_field = self.create_field(cond_input)

                cases = input_tool.get('cases', '')
                for case in cases:

                    case_inputs = case.get('inputs')
                    if case_inputs:
                        case_value = case.get('value')
                        nested_field.append(Div(id='id_'+new_name_field+"_cond_"+case_value,
                                                *self.parse_galaxy_input_tool(case_inputs)))

                fields_created.append(Div(new_name_field, *nested_field))

            else:
                fields_created.append(Field(self.create_field(input_tool)))

        return fields_created

    def __init__(self, tool_params, *args, **kwargs):
        super(ToolForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper(self)
        self.helper.form_class = 'blueForms'

        self.formset = self.parse_galaxy_input_tool(tool_params)
        FIELDS = self.formset

        self.helper.layout = Layout(FormActions(Field(*FIELDS),
                                                Submit('submit', 'Submit'),
                                            )
                                    )
