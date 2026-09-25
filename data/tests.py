import io
from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.http import HttpResponse
from django.test import TestCase, override_settings

from galaxy.models import GalaxyUser, Server
from utils import biofile
from data.views import UploadMixin


class ValidFastaTest(TestCase):
    """
    Regression test for utils.biofile.valid_fasta(): it used to pass an
    uploaded file's binary handle straight to Bio.SeqIO.parse(), which
    reads bytes lines. Biopython's SimpleFastaParser detects end-of-file
    by comparing a line to "" (str) - a check that never matches b""
    under Python 3, so parsing a real uploaded fasta file always crashed
    with IndexError at the true end of the file, every time (this only
    worked in Python 2, where bytes == str). Only caught by actually
    submitting a real workflow through the live UI - none of this
    codebase's other tests exercise an uploaded (as opposed to pasted or
    on-disk) fasta file.
    """

    def test_uploaded_binary_fasta_is_parsed_without_crashing(self):
        fasta = SimpleUploadedFile(
            "test.fa",
            b">seq1\nACGTACGTACGT\n>seq2\nACGTACGTACGT\n",
            content_type="text/plain",
        )
        # seqaa is intentionally not asserted: a short A/C/G/T-only
        # sequence is inherently ambiguous with the (real) amino acid
        # codes Ala/Cys/Gly/Thr, which is a pre-existing property of
        # check_aa()'s heuristic, unrelated to what's being tested here.
        nseq, length, seqaa = biofile.valid_fasta(fasta)
        self.assertEqual(nseq, 2)
        self.assertEqual(length, 12)


class SanitizeFastaContentTest(TestCase):
    """
    Regression test for utils.biofile.sanitize_fasta_id/
    sanitize_fasta_content: a real user's uploaded FASTA had a
    non-breaking space (U+00A0) embedded inside a sequence id
    ("A0A1Q2MHV5\\xa0_1_364" - plausibly a copy-paste artifact from a
    webpage/spreadsheet). Nothing in the regular OneClick/Advanced/Tool
    submission pipeline sanitized sequence ids at all before this
    (cleanseqname() was only ever called from blast/tasks.py, a
    separate code path) - it survived MAFFT/PhyML/PhyML-SMS untouched
    (those tools don't treat a NBSP as a token delimiter) into the
    output tree, where newick_utilities' own Newick parser *does* treat
    it as one: "ERROR: missing ')' at line 0 near '_1_364'" on the
    "Tree image" step - reported live, with the actual .nhx file
    containing 12 such ids. A more lenient parser (gotree) read the
    same tree back with no complaint.
    """

    def test_sanitize_fasta_id_replaces_non_breaking_space(self):
        self.assertEqual(
            biofile.sanitize_fasta_id('A0A1Q2MHV5\xa0_1_364'),
            'A0A1Q2MHV5_1_364')

    def test_sanitize_fasta_id_replaces_newick_special_characters(self):
        self.assertEqual(
            biofile.sanitize_fasta_id('foo(bar);baz,qux:1'),
            'foo_bar_baz_qux_1')

    def test_sanitize_fasta_id_collapses_and_strips_underscores(self):
        self.assertEqual(biofile.sanitize_fasta_id('a  b'), 'a_b')
        self.assertEqual(biofile.sanitize_fasta_id('_a_'), 'a')

    def test_sanitize_fasta_id_never_returns_empty(self):
        self.assertEqual(biofile.sanitize_fasta_id('   '), 'seq')

    def test_sanitize_fasta_content_rewrites_only_the_id_keeps_description(self):
        content = '>A0A1Q2MHV5\xa0_1_364 some description here\nACGT\n'
        result = biofile.sanitize_fasta_content(content)
        self.assertEqual(
            result,
            '>A0A1Q2MHV5_1_364 some description here\nACGT\n')

    def test_sanitize_fasta_content_leaves_sequence_lines_untouched(self):
        content = '>seq1\nACGT\n>seq2\nACGT\n'
        self.assertEqual(biofile.sanitize_fasta_content(content), content)

    def test_sanitize_fasta_content_accepts_and_returns_bytes(self):
        content = '>A\xa0B\nACGT\n'.encode('utf-8')
        result = biofile.sanitize_fasta_content(content)
        self.assertIsInstance(result, bytes)
        self.assertEqual(result, b'>A_B\nACGT\n')

    def test_sanitize_fasta_content_does_not_change_sequence_counts(self):
        # The actual regression: sanitizing ids must never change what
        # valid_fasta() counts - only the id characters change.
        raw = '>a\xa01\nACGT\n>b\xa02\nACGT\n>c\xa03\nACGT\n>d\xa04\nACGT\n'
        nseq_before, length_before, _ = biofile.valid_fasta(io.StringIO(raw))
        sanitized = biofile.sanitize_fasta_content(raw)
        nseq_after, length_after, _ = biofile.valid_fasta(io.StringIO(sanitized))
        self.assertEqual(nseq_before, nseq_after)
        self.assertEqual(length_before, length_after)


class UploadMixinSanitizesFastaTest(TestCase):
    """
    Confirms the sanitizer from SanitizeFastaContentTest above is
    actually wired into UploadMixin.upload_content()/upload_file() -
    the shared choke point both the plain /data/upload page and
    OneClick's own WorkflowFormView.form_valid() (workflows/views/
    generic.py, which inherits these methods) submit through - not
    just that the sanitizer function works in isolation.
    """

    def _mixin(self):
        mixin = UploadMixin()
        mixin.request = Mock()
        return mixin

    def test_upload_content_sanitizes_pasted_text(self):
        mixin = self._mixin()
        mixin.upload_content(
            '>A0A1Q2MHV5\xa0_1_364\nACGT\n', history_id='hist1')
        sent = mixin.request.galaxy.tools.paste_content.call_args.kwargs['content']
        self.assertIn('>A0A1Q2MHV5_1_364\n', sent)
        self.assertNotIn('\xa0', sent)

    def test_upload_file_sanitizes_an_uploaded_file(self):
        mixin = self._mixin()
        fasta = SimpleUploadedFile(
            'test.fa', '>A0A1Q2MHV5\xa0_1_364\nACGT\n'.encode('utf-8'),
            content_type='text/plain')
        captured = {}

        def _capture(path, **kwargs):
            # Read the temp file *while the mocked upload_file() call is
            # still in progress* - it's a NamedTemporaryFile (delete on
            # close), gone by the time upload_file() returns and its
            # own local variable goes out of scope.
            with open(path, 'rb') as f:
                captured['content'] = f.read().decode('utf-8')

        mixin.request.galaxy.tools.upload_file.side_effect = _capture
        mixin.upload_file(fasta, history_id='hist1')
        self.assertIn('>A0A1Q2MHV5_1_364\n', captured['content'])
        self.assertNotIn('\xa0', captured['content'])


class StaticPagesSmokeTest(TestCase):
    """
    Basic smoke tests: these pages are plain TemplateViews with no
    Galaxy/DB dependency, so they must render successfully in any
    environment (including CI, which has no Galaxy server or Redis).
    """

    def test_pages_return_200(self):
        paths = [
            '/',
            '/about',
            '/documentation',
            '/analysis',
            '/status',
            # Exercises real crispy_forms + django-simple-captcha
            # rendering (not just app loading), unlike the other pages.
            '/about/feedback',
        ]
        for path in paths:
            response = self.client.get(path)
            self.assertEqual(
                response.status_code, 200,
                "GET %s returned %s" % (path, response.status_code))

    def test_workspace_nav_link_is_always_clickable(self):
        """
        Regression test: templates/base.html's "Workspace" nav link used
        to only get a real href when request.session['histories'] had
        at least one entry - an empty session (a brand new visitor, or
        a logged-in account whose analyses only exist as
        WorkspaceHistory.user rows, not session ids - see
        PreviousHistoryListView, workspace/views.py) rendered it as
        plain muted, unclickable text instead. The Workspace page itself
        already handles an empty/account-only case gracefully (its own
        {% empty %} block, or the account-owned rows merged in by
        PreviousHistoryListView), so gating the link itself on session
        state was never actually necessary - just made the page
        unreachable from the nav in exactly the cases where a logged-in
        user would most want to reach it from a fresh session.
        """
        response = self.client.get('/')
        self.assertContains(response, 'href="/workspace/histories"')
        self.assertNotContains(response, 'text-muted')

    def test_feedback_form_renders_exactly_one_captcha_widget(self):
        """
        Regression test: surveys.forms.FeedbackForm used to append a
        Field('captcha ', ...) (trailing-space typo) onto
        FormHelper(self)'s own auto-built default layout, which already
        includes every form field - including a correctly-named
        'captcha' entry. crispy_forms's FAIL_SILENTLY handling makes an
        unresolvable field name log a warning and render "" rather than
        raise, so test_pages_return_200 above (a plain 200-status check)
        never caught this: the page always returned 200, it just quietly
        never rendered a usable captcha input at all, so no real visitor
        could ever pass this form's captcha check. Naively fixing just
        the typo surfaced a second bug from the same line: the field
        was then genuinely resolved *twice* (once by the implicit
        default layout, once by the explicit append), rendering two
        captcha widgets on the page - caught by actually reloading the
        page and counting the widgets, not assumed. Both fixed together
        by building the layout explicitly instead of appending onto the
        implicit default.
        """
        response = self.client.get('/about/feedback')
        html = response.content.decode()
        self.assertEqual(html.count('id="id_captcha_1"'), 1)
        self.assertEqual(html.count('name="captcha_0"'), 1)
        self.assertNotIn('Could not resolve form field', html)


class MaintenanceModeTest(TestCase):
    """
    NGPhylogeny_fr.middleware.MaintenanceModeMiddleware: every request
    should get templates/maintenance.html (503) instead of routing
    normally whenever settings.NGPHYLO_MAINTENANCE_MODE is True, and
    normal routing must be completely unaffected when it's False (the
    default).
    """

    @override_settings(NGPHYLO_MAINTENANCE_MODE=True)
    def test_maintenance_mode_on_serves_maintenance_page_for_any_path(self):
        # Including a path that wouldn't otherwise resolve at all - the
        # middleware runs before URL routing, so it overrides even what
        # would normally be a 404, not just real pages.
        for path in ['/', '/about', '/this-path-does-not-exist']:
            response = self.client.get(path)
            self.assertEqual(response.status_code, 503, path)
            self.assertTemplateUsed(response, 'maintenance.html')

    @override_settings(NGPHYLO_MAINTENANCE_MODE=True)
    def test_status_endpoint_stays_up_during_maintenance(self):
        """
        Regression test: manifest.yaml's web Deployment points both its
        readinessProbe and livenessProbe at /status. Before this
        exemption, enabling maintenance mode made those probes see the
        same 503 the middleware sends everywhere else, read that as "the
        container is broken" rather than "intentionally in maintenance",
        and endlessly restart it via failed liveness checks - the
        Service ended up with zero ready endpoints and the site went
        fully down, the opposite of the intended effect. Only caught by
        actually enabling MAINTENANCE=true against a real deployment.
        """
        response = self.client.get('/status')
        self.assertEqual(response.status_code, 200)
        self.assertTemplateNotUsed(response, 'maintenance.html')

    @override_settings(NGPHYLO_MAINTENANCE_MODE=False)
    def test_maintenance_mode_off_routes_normally(self):
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'home.html')


class DisplayViewsAjaxDetectionTest(TestCase):
    """
    Regression test: display_file/display_params/display_msa
    (data/views.py) used to call request.is_ajax(), a method Django
    removed in 3.1 - every real (non-AJAX) page load of a dataset's
    display/params/MSA page crashed with "AttributeError: 'WSGIRequest'
    object has no attribute 'is_ajax'" (reported live against
    /data/display/<id> on the deployed instance; none of these views had
    a test before). Fixed by checking the X-Requested-With header
    directly, which is what is_ajax() itself used to do.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)

    def test_display_file_non_ajax_renders_template_without_crashing(self):
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'history_id': 'hist1'}):
            response = self.client.get('/data/display/ff5476bcf6c921fa')
        self.assertEqual(response.status_code, 200)

    def test_display_file_ajax_dispatches_to_display_raw_without_crashing(self):
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'history_id': 'hist1'}), \
             patch('data.views.display_raw',
                   return_value=HttpResponse('raw content')) as mock_raw:
            response = self.client.get(
                '/data/display/ff5476bcf6c921fa',
                HTTP_X_REQUESTED_WITH='XMLHttpRequest')
        mock_raw.assert_called_once()
        self.assertEqual(response.status_code, 200)

    def test_display_params_non_ajax_renders_template_without_crashing(self):
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'creating_job': 'job1',
                                  'history_id': 'hist1'}):
            response = self.client.get('/data/params/ff5476bcf6c921fa')
        self.assertEqual(response.status_code, 200)

    def test_display_params_ajax_returns_job_json_without_crashing(self):
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'creating_job': 'job1'}), \
             patch('bioblend.galaxy.jobs.JobsClient.show_job',
                   return_value={'tool_id': 'mytool'}):
            response = self.client.get(
                '/data/params/ff5476bcf6c921fa',
                HTTP_X_REQUESTED_WITH='XMLHttpRequest')
        self.assertEqual(response.status_code, 200)
        self.assertIn(b'mytool', response.content)

    def test_display_msa_non_ajax_renders_template_without_crashing(self):
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'history_id': 'hist1'}):
            response = self.client.get('/data/displaymsa/ff5476bcf6c921fa')
        self.assertEqual(response.status_code, 200)


class DownloadFileContentTest(TestCase):
    """
    Regression test: download_file (data/views.py) used to wrap the
    downloaded dataset bytes in a StreamingHttpResponse directly.
    Iterating a bytes object in Python 3 yields one int per byte (unlike
    Python 2, where iterating a str/bytes yielded one character at a
    time) - StreamingHttpResponse then writes out each int as its own
    chunk, so real content came out completely garbled: b'Hello World'
    was served as b'721011081081113287111114108100' (each byte's decimal
    value concatenated) instead of the actual bytes. Reported live via
    /data/display/<id> showing "hex/binary"-looking garbage instead of
    the real dataset content; not covered by any existing test (the
    AJAX-detection tests above mock display_raw/download_file itself, so
    they never exercised this).
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)

    def test_downloaded_content_is_served_intact_not_garbled(self):
        original = b'Hello World\nThis is a real dataset.\n'
        fake_response = Mock()
        fake_response.read.return_value = original
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'download_url': '/api/datasets/x/display',
                                  'name': 'result', 'file_ext': 'txt',
                                  'history_id': 'hist1'}), \
             patch('data.views.urlopen', return_value=fake_response):
            response = self.client.get(
                '/data/display/ff5476bcf6c921fa',
                HTTP_X_REQUESTED_WITH='XMLHttpRequest')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content, original)

    def test_purged_dataset_download_url_404_does_not_crash_the_view(self):
        """
        Regression test: hit live in production - show_dataset() can
        still succeed and hand back a download_url for a dataset whose
        actual file content has already been purged from Galaxy
        (workspace.tasks.deleteoldgalaxyhistory()'s
        delete_history(..., purge=True), see CLAUDE.md) - this endpoint
        has no ownership check, so an old dataset id stays directly
        reachable. urlopen(req) then raised an uncaught
        urllib.error.HTTPError: HTTP Error 404: Not Found, 500ing the
        page instead of degrading gracefully.
        """
        from urllib.error import HTTPError
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'download_url': '/api/datasets/x/display',
                                  'name': 'result', 'file_ext': 'txt',
                                  'history_id': 'hist1'}), \
             patch('data.views.urlopen',
                   side_effect=HTTPError(
                       'http://fake-galaxy.example.org/api/datasets/x/display',
                       404, 'Not Found', None, None)):
            response = self.client.get('/data/download/ff5476bcf6c921fa')
        self.assertEqual(response.status_code, 200)


class TreeVisualizationTest(TestCase):
    """
    Regression test: tree_visualization (data/views.py) used to pass the
    downloaded newick bytes straight into the template context.
    treeviz/tree.html embeds it in a JS string literal via
    {{ newick_tree|escapejs }} - Django's template rendering calls str()
    on a non-string context value, which for bytes produces the Python
    repr (e.g. "b'(A:0.1,B:0.2);\\n'", literal b-quote-backslash-n and
    all) instead of the actual tree text, so every tree visualization
    page loaded with garbage instead of a real Newick string. Same root
    cause class as DownloadFileContentTest above, different code path
    (template context, not a raw HttpResponse) - not covered by any
    existing test.
    """

    def setUp(self):
        with patch('galaxy.models.requests.get',
                   return_value=Mock(status_code=200,
                                      json=lambda: {'version_major': '25.1'})):
            self.server = Server.objects.create(
                url='http://fake-galaxy.example.org', current=True)
        user = User.objects.create_user('admin')
        GalaxyUser.objects.create(
            user=user, galaxy_server=self.server, api_key='fakekey',
            anonymous=True)

    def test_newick_tree_is_rendered_as_real_text_not_a_bytes_repr(self):
        newick = b'(A:0.1,(B:0.2,C:0.3):0.15);'
        fake_response = Mock()
        fake_response.read.return_value = newick
        with patch('bioblend.galaxy.datasets.DatasetClient.show_dataset',
                   return_value={'download_url': '/api/datasets/x/display',
                                  'history_id': 'hist1'}), \
             patch('data.views.urlopen', return_value=fake_response):
            response = self.client.get('/data/displaytree/ff5476bcf6c921fa')

        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        # escapejs legitimately turns the trailing ";" into a unicode
        # escape for JS-string safety - check the actual tree structure
        # survived intact rather than the exact (escaping-dependent)
        # punctuation.
        self.assertIn('(A:0.1,(B:0.2,C:0.3):0.15)', content)
        self.assertNotIn("b'(A:0.1", content)
