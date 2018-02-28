from django import forms
from django.utils.text import slugify

from tools.models import ToolFieldWhiteList
from tools.forms import ToolForm


def tool_form_factory(tool):
    """
    toolform class factory
    :param tool:
    :param prefix:
    :return: Form class
    """

    tool_inputs_details = tool.fetch_tool_json()
    tool_field_white_list, created = ToolFieldWhiteList.objects.get_or_create(tool=tool, context="w")
    class_name = str(slugify(tool.name).title().replace('-', '')) + 'Form'
    prefix = class_name.lower()
    return type(class_name, (ToolForm,),
                {
                    'tool_params': tool_inputs_details.get('inputs'),
                    'tool_id': tool.id_galaxy,
                    'tool_name': tool.name,
                    'prefix': prefix,
                    'visible_field': tool_field_white_list.saved_params,
                    'fields_ids_mapping': {},
                    'n': 0,
                }
                )
