from django.test import TestCase, override_settings

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
