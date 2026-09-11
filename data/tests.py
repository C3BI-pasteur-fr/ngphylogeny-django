from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.http import HttpResponse
from django.test import TestCase

from galaxy.models import GalaxyUser, Server
from utils import biofile


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
