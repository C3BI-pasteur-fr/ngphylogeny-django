from django.core.management.base import BaseCommand

from tools.models import ToolInputData
from tools.models import ToolOutputData


class Command(BaseCommand):
    help = 'Compute all Galaxy tools link'
    links = []

    def add_arguments(self, parser):
        parser.add_argument('--linkfile')

    def read_link_file(self, linkfile):
        with open(linkfile) as f:
            for line in f:
                link = line.strip().split("\t")
                self.links.append(link)

    def handle(self, *args, **options):
        linkfile = options.get('linkfile')
        print(linkfile)
        self.read_link_file(linkfile)
        for io in self.links:
            print(io)
            first = ToolOutputData.objects.filter(
                tool__name=io[0], name=io[1]).first()
            second = ToolInputData.objects.filter(
                tool__name=io[2], name=io[3]).first()
            first.compatible_inputs.add(second)
