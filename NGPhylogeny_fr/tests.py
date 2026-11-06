from unittest.mock import patch

from django import forms
from django.template.loader import render_to_string
from django.test import RequestFactory, TestCase, override_settings
from django.urls import reverse

from NGPhylogeny_fr.context_processors import git_commit


class GitCommitContextProcessorTest(TestCase):
    """
    NGPhylogeny_fr.context_processors.git_commit - exposes
    settings.NGPHYLO_GIT_COMMIT (baked into a plain GIT_COMMIT file at
    Docker build time, see Dockerfile/settings/base.py's own comments
    on this) to every template, for the footer's version link.
    """

    @override_settings(NGPHYLO_GIT_COMMIT='')
    def test_empty_commit_gives_no_short_sha_or_url(self):
        self.assertEqual(git_commit(None), {
            'git_commit_short': '', 'git_commit_url': ''})

    @override_settings(
        NGPHYLO_GIT_COMMIT='9b81d5aabc1234567890abcdef1234567890abcd')
    def test_real_commit_is_truncated_and_linked_to_the_real_github_repo(self):
        result = git_commit(None)
        self.assertEqual(result['git_commit_short'], '9b81d5a')
        self.assertEqual(
            result['git_commit_url'],
            'https://github.com/C3BI-pasteur-fr/ngphylogeny-django/'
            'commit/9b81d5aabc1234567890abcdef1234567890abcd')

    @override_settings(NGPHYLO_GIT_COMMIT='')
    def test_home_page_shows_no_version_line_when_commit_is_empty(self):
        response = self.client.get('/')
        self.assertNotContains(response, 'version <a href=')

    @override_settings(
        NGPHYLO_GIT_COMMIT='9b81d5aabc1234567890abcdef1234567890abcd')
    def test_home_page_shows_the_version_link_when_commit_is_set(self):
        response = self.client.get('/')
        self.assertContains(
            response,
            'https://github.com/C3BI-pasteur-fr/ngphylogeny-django/'
            'commit/9b81d5aabc1234567890abcdef1234567890abcd')
        self.assertContains(response, '9b81d5a')


class UploadSizeLimitTest(TestCase):
    LIMIT = 51 * 1024 * 1024

    def _post(self, content_length):
        return self.client.post(
            reverse('home'), data='', content_type='application/octet-stream',
            CONTENT_LENGTH=content_length)

    def test_oversized_request_is_refused_before_its_body_is_read(self):
        with patch('django.middleware.csrf.CsrfViewMiddleware.process_view') \
                as csrf:
            response = self._post(self.LIMIT + 1)
        csrf.assert_not_called()
        self.assertEqual(response.status_code, 413)
        self.assertTemplateUsed(response, '413.html')
        self.assertContains(response, '50 MB', status_code=413)

    def test_a_file_of_exactly_the_cap_gets_through(self):
        self.assertNotEqual(self._post(self.LIMIT).status_code, 413)

    def test_unparseable_content_length_is_left_alone(self):
        self.assertNotEqual(self._post('abc').status_code, 413)

    def test_browser_check_uses_the_option_the_plugin_actually_reads(self):
        request = RequestFactory().get('/')
        html = render_to_string(
            'workflows/workflows_form.html', {'form': forms.Form()},
            request=request)
        self.assertIn('maxSize: 50', html)
        self.assertNotIn('fileMaxSize', html)

    def test_forms_without_the_plugin_check_the_size_in_the_browser(self):
        request = RequestFactory().get('/')
        for template, context in (
                ('tools/tool_form.html', {'toolform': forms.Form()}),
                ('workflows/workflows_adv_singlepage_form.html', {}),
                ('workflows/workflows_maker_form.html', {})):
            html = render_to_string(template, context, request=request)
            self.assertIn('total > 50 * 1024 * 1024', html, template)
