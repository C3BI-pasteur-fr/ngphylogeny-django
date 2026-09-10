from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase

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
