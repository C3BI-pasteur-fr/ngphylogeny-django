import urllib

import requests
from django.conf import settings
from django.core.management.base import BaseCommand

from account.models import GalaxyServer
from tools.models import Tool, ToolDataInput, ToolDataOutput


class Command(BaseCommand):

    help = 'Import Galaxy tools to NGPhylogeny'

    def handle(self, *args, **options):

        self.handle_noargs()

    def handle_noargs(self):

        galaxy_server, created = GalaxyServer.objects.get_or_create(url=settings.GALAXY_SERVER_URL)

        params = urllib.urlencode({'q': "phylogeny"}, True)
        tools_url = '%s/%s/%s/?%s' % (galaxy_server.url, 'api', 'tools', params)
        connection = requests.get(tools_url)
        tools_ids = []
        if connection.status_code == 200:
            tools_ids = connection.json()

            for id_tool in tools_ids:
                params = urllib.urlencode({'io_details': "true"}, True)
                tool_url = '%s/%s/%s/%s/?%s' % (galaxy_server.url, 'api', 'tools', id_tool, params )
                tool_info_request = requests.get(tool_url)
                tool_info = tool_info_request.json()

                t, created = Tool.objects.get_or_create(id_galaxy=id_tool, galaxy_server=galaxy_server)

                if created:
                    t.name = tool_info.get('name')
                    print t.name
                    t.description = tool_info.get('description')
                    t.version = tool_info.get('version')

                    if "toolshed" in t.id_galaxy:
                        t.toolshed = t.id_galaxy.split('/')[0]
                    t.save()

                    inputs_tools = tool_info.get('inputs')
                    for input_d in inputs_tools:

                        if input_d.get('type') == 'data':

                            edam = input_d.get('edam')
                            ed_format = ""
                            if edam:
                                ed_format = edam.get('edam_formats')

                            input_obj = ToolDataInput(name=input_d.get('name'),
                                                      edam_formats=ed_format,
                                                      extensions=input_d.get('extensions'),
                                                      tool=t
                                                      )
                            input_obj.save()
                            print "inputs:", input_d.get('name'), input_d.get('extensions'), input_d.get('edam')

                    outputs_tools = tool_info.get('outputs')

                    for output_d in outputs_tools:
                        output_obj = ToolDataOutput(name=output_d.get('name'),
                                                    edam_formats=output_d.get('edam_format'),
                                                    extensions=output_d.get('format'),
                                                    type="o",
                                                    tool=t
                                                    )
                        output_obj.save()
                        print "output", output_d.get('name'), output_d.get('format'), output_d.get('edam_format')


