from django.core.management.base import BaseCommand

from tools.models import ToolInputData


class Command(BaseCommand):
    help = 'Compute all Galaxy tools link'

    def add_arguments(self, parser):
        # Named (optional) arguments
        parser.add_argument('--full')
        parser.add_argument('--edam')
        parser.add_argument('--ignore', nargs='+', help='ignored file extension')


    def handle(self, *args, **options):

        full = options.get('full')
        ignored_ext = options.get('ignore')

        for t_input in ToolInputData.objects.all():

            tool_outputs = t_input.search_compatible_outputs(ignore=ignored_ext)

            for t_output in tool_outputs:
                print str(t_input) + " -> " + str(t_output)
                t_output.compatible_inputs.add(t_input)
