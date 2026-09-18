import ast
import logging
import requests
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode

from django.core.exceptions import ValidationError
from django.db import models
from django.utils.functional import cached_property
from django.utils.translation import gettext as _

from data.models import ExampleFile
from galaxy.models import Server

logger = logging.getLogger(__name__)


class Tool(models.Model):
    """

    """
    galaxy_server = models.ForeignKey(Server, on_delete=models.CASCADE)
    id_galaxy = models.CharField(max_length=250)
    toolshed = models.CharField(max_length=100, null=True, blank=True)
    name = models.CharField(max_length=100, blank=True)
    version = models.CharField(max_length=20, blank=True)
    description = models.CharField(max_length=250)
    toolshed_revision = models.CharField(max_length=250, null=True, blank=True)
    visible = models.BooleanField(
        default=True, help_text="Display this tool on the user web interface")
    oneclick = models.BooleanField(default=False)
    rank = models.IntegerField(default=0, help_text="tool order")
    max_nbseq = models.IntegerField(default=-1, help_text="max length")
    max_boot = models.IntegerField(default=-1, help_text="max boot replicates")
    max_lengthxnbseqsquared = models.IntegerField(default=-1, help_text="max length x nbseq^2")
    max_nbseqsquaredxboot = models.IntegerField(default=-1, help_text="max nbseq^2 x boot")
    max_lengthxnbseqsquaredxboot = models.IntegerField(default=-1, help_text="max length x nbseq^2 x boot")
    # If input data is aa : limits are divided by this scaling factor
    aa_scale_factor =  models.IntegerField(default=1, help_text="limit scaling for amino acid")

    @property
    def toolflags(self):
        return ",".join(self.toolflag_set.all().values_list('verbose_name', flat=True))

    @property
    def citations(self):
        f = []
        for c  in self.citation_set.all():
            for d in c.format():
                f.append( d )
        return f

    @property
    def get_params_detail(self):
        return self.tool_json.get('inputs')

    def clean(self):

        if 'err_code' in self.tool_json.keys():
            raise ValidationError({
                'id_galaxy': ValidationError(_("%(err_msg)s"), code=self.tool_json.get('err_code'),
                                             params={'err_msg': self.tool_json.get('err_msg'), })})
        else:
            self.import_tool(self.tool_json)

    def can_run_on_data(self, nseq, length, nboot, seqaa):
        """
        Returns true if the tool can is authorized to run on 
        data of the given size false otherwise
        """
        if (self.max_nbseq > 0 and
            ((nseq > self.max_nbseq) or
             (seqaa and nseq > self.max_nbseq//self.aa_scale_factor))):
            return False
        if self.max_boot > 0 and nboot > self.max_boot :
            return False
        if (self.max_lengthxnbseqsquared > 0 and
            ((length*nseq > self.max_lengthxnbseqsquared) or
             (seqaa and length*nseq > self.max_lengthxnbseqsquared//self.aa_scale_factor))):
            return False
        if self.max_nbseqsquaredxboot > 0 and nseq*nboot > self.max_nbseqsquaredxboot:
            return False
        if (self.max_lengthxnbseqsquaredxboot > 0 and
            ((length*nseq*nboot > self.max_lengthxnbseqsquaredxboot) or
             (seqaa and length*nseq*nboot > self.max_lengthxnbseqsquaredxboot//self.aa_scale_factor))):
            return False
        return True

    def save(self, *args, **kwargs):

        self.clean()
        new = self._state.adding  # _state change after save()
        super(Tool, self).save(*args, **kwargs)

        if new:
            """fetch and created tool io information"""
            self.import_tool_io(self.tool_json)

    class Meta:
        unique_together = (('galaxy_server', 'id_galaxy'),)

    @classmethod
    def import_tools_from_url(cls, galaxy_url, query="phylogeny"):
        try:
            galaxy_server = Server.objects.get(url=galaxy_url)
        except:
            raise Exception("NGPhylogeny server is not properly configured,"
                            "Please ensure that the Galaxy server is correctly set up")
        return cls.import_tools(galaxy_server, query)

    @cached_property
    def tool_json(self):
        return self.fetch_tool_json()

    def search_tools_json(self, galaxy_server='', query="phylogeny"):
        """fetch and return list of id tool from galaxy"""

        galaxy_server = galaxy_server or self.galaxy_server.url

        tools_url = '%s/%s/%s/%s' % (galaxy_server,
                                     'api', 'tools', self.id_galaxy)

        tool_panels = requests.get(tools_url).json()

        for panel in tool_panels:
            if query in panel.get('name').lower():
                return [t.get("id") for t in panel.get('elems')]
        return []

    def fetch_tool_json(self):
        """
            fetch and store tool information
        """
        tool_url = '%s/%s/%s/%s' % (self.galaxy_server.url,
                                    'api', 'tools', self.id_galaxy)
        tool_info_request = requests.get(
            tool_url, params={'io_details': "true"})

        return tool_info_request.json()

    def import_tool(self, tool_info):
        """
        Extract tool information from Galaxy API
        :param dict tool_info: json galaxy tool information
        """
        # tool_info = self.fetch_tool_json()

        self.name = self.name or tool_info.get('name')
        self.description = self.description or tool_info.get('description')
        self.version = self.version or tool_info.get('version')

        toolshed = tool_info.get('tool_shed_repository')

        if toolshed:
            self.toolshed = self.toolshed or toolshed.get('tool_shed')
            self.toolshed_revision = self.toolshed_revision or toolshed.get(
                'changeset_revision')

        return self

    def import_tool_io(self, tool_info):
        """
        Extract tool input/output information from Galaxy API
        :param dict tool_info: json galaxy tool information
        """

        # Inputs
        inputs_tools = tool_info.get('inputs')
        inputs_list = []
        for input_d in inputs_tools:

            if input_d.get('type') == 'data':

                edam = input_d.get('edam')
                ed_format = ""
                if edam:
                    ed_format = edam.get('edam_formats')

                input_obj, created_input = ToolInputData.objects.get_or_create(name=input_d.get('name'),
                                                                               tool=self
                                                                               )
                if created_input:
                    input_obj.edam_formats = ed_format
                    input_obj.extensions = input_d.get('extensions')
                    input_obj.save()

            else:
                # remove input data from whitelist field for workflow context
                # (ie: without data input because on workflow data input come from previous steps)
                inputs_list.append(input_d.get('name'))

        # create whitelist field for workflows
        w, created_w = ToolFieldWhiteList.objects.get_or_create(
            tool_id=self.id, context='w')
        if created_w:
            w._params = ",".join(inputs_list)
            w.save()

        # Outputs
        outputs_tools = tool_info.get('outputs')
        for output_d in outputs_tools:
            output_obj, created_output = ToolOutputData.objects.get_or_create(name=output_d.get('name'),
                                                                              tool=self
                                                                              )
            # TODO issue Galaxy
            if created_output:
                output_obj.edam_format = output_d.get('edam_format')

                if output_d.get('format'):

                    output_obj.extension = output_d.get('format')

                elif output_d.get('model_class') == "ToolOutputCollection":

                    output_obj.extension = "DataCollection"
                else:
                    raise NotImplementedError

                output_obj.save()

    @classmethod
    def import_tools(cls, galaxy_server, tools=None, query="phylogeny", force=False):
        """
        Import a list of tool
        :param galaxy_server: GalaxyServer instance
        :param tools: list of tool id form Galaxy
        :param query: keyword to find specific tools
        :param force: force update
        :return type: dict =   {
                                  'new':[],           #list of tool
                                  'already_exist':[]  #list of tool
                                  'error': [],        #list of id_tool not valid
                                }
        """
        query = query.lower()
        tools_ids = tools or []
        tools_import_report = {
            'new': [],
            'already_exist': [],
            'error': []
        }

        if not tools_ids:
            # fetch tools list
            tools_url = '%s/%s/%s/' % (galaxy_server.url, 'api', 'tools')
            connection = requests.get(tools_url, params={'in_panel': "false"})
            tools_list = connection.json()
            for tool in tools_list:
                _m = [
                    tool.get('id').lower(),
                    tool.get('name').lower(),
                    str(tool.get('panel_section_name')).lower()
                ]
                if query in _m:
                    tools_ids.append(tool.get('id'))

        while tools_ids:
            id_tool = tools_ids.pop()
            try:
                t, created = Tool.objects.get_or_create(
                    id_galaxy=id_tool, galaxy_server=galaxy_server)
                t.save()
                if force:
                    t.import_tool_io(t.tool_json)
                if created or force:
                    # Citations are only (re-)fetched for newly-created or
                    # force-reimported tools - but docker/init.sh always
                    # calls importtools with --force on every container
                    # start, so this branch runs on every single redeploy
                    # for every tool. Replace (clear then re-insert)
                    # rather than blindly append: appending here used to
                    # duplicate every Citation row on every redeploy -
                    # real production tools ended up with 25-75 duplicate
                    # rows of the same 1-3 actual citations.
                    cite_url = '%s/%s/%s/%s/%s' % (galaxy_server.url, 'api', 'tools', id_tool, 'citations')
                    connection = requests.get(cite_url)
                    citations = connection.json()
                    t.citation_set.all().delete()
                    for cite in citations:
                        # No encode/decode roundtrip needed: unlike Python 2's
                        # requests/json stack (which this iso-8859-1->utf8
                        # roundtrip used to correct for), Python 3's
                        # requests.json() already returns correctly-decoded
                        # text - re-encoding it as iso-8859-1 just raises
                        # UnicodeEncodeError on any citation containing a
                        # character outside Latin-1 (en dashes, curly quotes,
                        # ...), which is common in real citation text.
                        c = Citation(reference=cite.get('content',''), tool=t)
                        c.save()
                    tools_import_report['new'].append(t)
                else:
                    tools_import_report['already_exist'].append(t)
            except (ValueError, ValidationError) as e:
                print(e)
                tools_import_report['error'].append(id_tool)
        return tools_import_report

    @property
    def compatible_tool(self):
        """
            return the list of tools that have inputs compatible with the outputs of the tool
            to create workflow.

            ie: tool.output.format == next_tool.input.format
        """
        outputs = ToolOutputData.objects.filter(tool=self)
        outputs_formats = outputs.values_list('edam_formats', flat=True)

        tools_data_input = ToolInputData.objects.filter(tool=self)
        tools_compatible = []

        from ast import literal_eval
        for data_input in tools_data_input:

            list_edam_formats = literal_eval(data_input.edam_formats)
            for ed in list_edam_formats:
                if ed in set(outputs_formats):
                    tools_compatible.append(data_input.tool.pk)

        return Tool.objects.filter(pk__in=tools_compatible)

    def __str__(self):
        return "{} - {}".format(self.name, self.version)


class ToolData(models.Model):
    """
        Generic model for Input/output data
    """

    name = models.CharField(max_length=250)
    edam_data = models.CharField(max_length=250, null=True, blank=True)

    class Meta:
        abstract = True
        unique_together = ['tool', 'name']


class ToolInputData(ToolData):
    """
    Tool input data information extract from Galaxy
    """
    edam_formats = models.CharField(max_length=250, null=True, blank=True)
    extensions = models.CharField(max_length=100)
    examplefile = models.ForeignKey(ExampleFile, null=True, blank=True, on_delete=models.SET_NULL)
    # Wether this field may be linked to the first input data step
    # in the workflow maker.
    # Avoids to link input data to all input file fields in PhyML for example
    # (input alignment + user input tree...)
    galaxy_input_data  = models.BooleanField(default=False)
    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)

    def get_extensions(self):
        return ast.literal_eval(self.extensions)

    def search_compatible_outputs(self, ignore=None):
        """
        Search compatible data produced by others tools
        :param ignore : list of ignored extension
        :return: list of <ToolOutputData>
        """
        ignore = ignore or []

        l_ext = self.get_extensions()
        l_ext_filtered = [ext for ext in l_ext if ext not in ignore]
        print(l_ext_filtered)
        galaxy_server = self.tool.galaxy_server
        return ToolOutputData.objects.filter(extension__in=l_ext_filtered,
                                             tool__galaxy_server=galaxy_server)

    def __str__(self):
        return "%s | %s: %s" % (self.tool, self.name, self.extensions)

    class Meta:
        verbose_name_plural = "Tool input data"


class ToolOutputData(ToolData):
    """
    Tool output data extract from Galaxy
    """
    edam_format = models.CharField(max_length=250, null=True, blank=True)
    extension = models.CharField(max_length=100)

    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)
    compatible_inputs = models.ManyToManyField(ToolInputData, blank=True)

    def search_compatible_inputs(self):
        return ToolInputData.objects.filter(extensions__contains=self.extension)

    def __str__(self):
        return "%s | %s: %s" % (self.tool, self.name, self.extension)

    class Meta:
        verbose_name_plural = "Tool output data"

class Citation(models.Model):
    """
    Tool references
    """
    # TextField, not a length-bounded CharField: this stores raw BibTeX
    # citation text straight from Galaxy's /api/tools/{id}/citations
    # (import_tools() in this module) - real citations (long abstracts,
    # many authors) routinely exceed any fixed length, and Postgres
    # enforces a CharField's max_length as a hard varchar() constraint at
    # the DB level, aborting the whole importtools run with
    # "value too long for type character varying(N)" the moment one
    # citation is too long, rather than failing just that one tool.
    reference = models.TextField(null=True, blank=True)
    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)

    @staticmethod
    def _bibtex_parser():
        """
        bibtexparser.loads() with no parser at all extracts a field's raw
        text with zero LaTeX interpretation - accented author names
        (Guindon et al.'s PhyML citation: "St{\\'{e}}phane"/
        "Jean-Fran{\\c{c}}ois") and brace-protected capitalization
        ("{PhyML}") both showed up completely unresolved and literal on
        the rendered page. convert_to_unicode (via bibtexparser's own
        latex_to_unicode) converts the whole standard set of LaTeX accent
        commands to real Unicode and strips any braces left over
        afterward - covers this class of artifact generally, rather than
        hand-rolling substitutions for every possible accented letter.
        """
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        return parser

    @staticmethod
    def _resolve_latex_escapes(reference):
        """
        Applied to the RAW BibTeX source text before parsing, not to
        already-parsed field values - and that ordering matters.
        convert_to_unicode (see _bibtex_parser()) only knows standard
        LaTeX commands, not $\\less$/$\\greater$ (a real but non-standard
        LaTeX-escaped way some bibliography tools wrap a small-caps tag
        to protect a word's capitalization from BibTeX's automatic
        title-casing - the Newick Utilities citation: "Unix" as
        U$\\less$scp$\\greater$nix$\\less$/scp$\\greater$, meaning
        U<scp>nix</scp>). Worse, run on the *parsed* field value,
        convert_to_unicode's own LaTeX interpretation reads \\l (inside
        \\less) as the real, standard LaTeX command for "ł" and silently
        mangles $\\less$ into $łess$ first - verified directly against
        bibtexparser 1.4.4, not assumed - so by the time a post-parse
        cleanup step would run, the pattern it's looking for is already
        gone. Resolving these escapes on the raw text first avoids the
        collision entirely: a bare "<"/">" character has no LaTeX meaning
        for convert_to_unicode to misinterpret.
        """
        if not reference:
            return reference
        return reference.replace('$\\less$', '<').replace('$\\greater$', '>')

    @staticmethod
    def _strip_scp_tags(text):
        """
        The only thing _resolve_latex_escapes' </>-resolution is ever
        used to wrap here is a <scp>...</scp> small-caps tag (see that
        method's docstring). Strip the tags themselves, applied to
        already-parsed field values, rather than leave literal tag markup
        in txt()'s plain-text output or rely on format()'s HTML output
        silently ignoring an unrecognized tag.
        """
        if not text:
            return text
        return text.replace('<scp>', '').replace('</scp>', '')

    def format(self):
        reference = self._resolve_latex_escapes(self.reference)
        try:
            bib_database = bibtexparser.loads(reference, parser=self._bibtex_parser())
        except Exception:
            # Raw BibTeX straight from Galaxy's /api/tools/{id}/citations
            # (see the reference field's own docstring) - a real one hit
            # live: a bare, non-standard month value ("month = june,"
            # rather than the 3-letter "jun" bibtexparser's common-strings
            # table actually knows, or a quoted/braced string) makes
            # bibtexparser try to resolve it as a @string macro reference
            # and raise bibtexparser.bibdatabase.UndefinedString - which,
            # unhandled, 500'd the entire history detail page over one
            # malformed citation on one tool. Caught broadly (not just
            # UndefinedString) since this is arbitrary external BibTeX
            # text with no schema guarantee - any other parse failure
            # shouldn't crash the page either.
            logger.warning(
                "Could not parse BibTeX citation %s for tool %s",
                self.pk, self.tool_id, exc_info=True)
            return []
        f = []
        for k, v in bib_database.entries_dict.items():
            journal = self._strip_scp_tags(v.get('journal',''))
            title = self._strip_scp_tags(v.get('title',''))
            year = v.get('year','')
            authors = self._strip_scp_tags(v.get('author',''))
            doi = v.get('doi','')
            volume = v.get('volume','')
            pages = v.get('pages','')
            html="""
            <a target="_blank" href="https://doi.org/%s">
            <div class="pub-date">%s</div>
            <div class="pub-content clear">
              <span class="pub-author-list">%s.</span>
              <span class="pub-title">%s.</span>
              <span class="pub-journal-name">%s, %s:%s </span>
            </div>
            </a>
            """ % (doi, year, authors, title, journal, volume, pages)
            f.append(html)                
        return f

    def txt(self):
        reference = self._resolve_latex_escapes(self.reference)
        try:
            bib_database = bibtexparser.loads(reference, parser=self._bibtex_parser())
        except Exception:
            # Same malformed-BibTeX defense as format() above.
            logger.warning(
                "Could not parse BibTeX citation %s for tool %s",
                self.pk, self.tool_id, exc_info=True)
            return ""
        f = ""
        for k, v in bib_database.entries_dict.items():
            journal = self._strip_scp_tags(v.get('journal',''))
            title = self._strip_scp_tags(v.get('title',''))
            year = v.get('year','')
            authors = self._strip_scp_tags(v.get('author',''))
            doi = v.get('doi','')
            volume = v.get('volume','')
            pages = v.get('pages','')
           
            f=f+authors + "("+year+"). "+title+" . "+journal+", "+volume+":"+pages+" - doi:"+doi+"\n"
        return f

    
class ToolFlag(models.Model):
    """
    Flag tools by category
    """

    name = models.CharField(max_length=5, unique=True)
    verbose_name = models.CharField(max_length=250, unique=True)
    tool = models.ManyToManyField(Tool)
    rank = models.IntegerField(default=999, help_text="flags order")

    def __str__(self):
        return self.verbose_name

    class Meta:
        ordering = ['rank', 'verbose_name']


class ToolFieldWhiteList(models.Model):
    """
    Models use to filter Galaxy render form
    """
    CONTEXT_CHOICES = (('w', 'Workflows'),
                       ('aw', 'Advanced Workflows'),
                       ('t', 'Tool Forms'),
                       ('at', 'Advanced Tool Forms'),
                       )

    DICT_CONTEXT_CHOICES = dict(CONTEXT_CHOICES)

    tool = models.ForeignKey(Tool, on_delete=models.CASCADE)
    context = models.CharField(max_length=2, choices=CONTEXT_CHOICES,
                               help_text='In which context use the whitelist')
    _params = models.TextField(max_length=1000, blank=True, default="")

    @property
    def saved_params(self):
        if self._params:
            return self._params.split(',')
        else:
            return []

    @property
    def get_params(self):
        l = []
        for p in self.tool.get_params_detail:
            l.append(p.get('name'))
            cond = p.get('test_param')
            if cond:
                l.append(cond.get('name'))

        return l

    @property
    def get_json_params(self):
        return self.tool.get_params_detail

    def __str__(self):
        return "%s ,%s" % (self.tool.name, self.DICT_CONTEXT_CHOICES.get(self.context))

    class Meta:
        unique_together = ('tool', 'context')
