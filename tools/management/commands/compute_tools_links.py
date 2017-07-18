from django.core.management.base import BaseCommand

from tools.models import ToolInputData


class Command(BaseCommand):
    help = 'Compute all Galaxy tools link'

    def handle(self, *args, **options):

        self.handle_noargs()

    def handle_noargs(self):

        for t_input in ToolInputData.objects.all():
            tool_outputs = t_input.search_compatible_outputs()
            for t_output in tool_outputs:
                t_output.compatible_inputs.add(t_input)
                t_output.save()
