from django.conf import settings

GITHUB_REPO_URL = 'https://github.com/C3BI-pasteur-fr/ngphylogeny-django'


def git_commit(request):
    """
    Exposes the running image's git commit (settings.NGPHYLO_GIT_COMMIT -
    see that setting's own comment for how it's baked in at Docker build
    time) to every template, for the footer's "Version" link. Empty
    settings.NGPHYLO_GIT_COMMIT (local dev, no CI build) means both
    values come back empty too - the footer's own {% if %} guard then
    just shows no version line, rather than a link to a nonexistent
    commit.
    """
    sha = settings.NGPHYLO_GIT_COMMIT
    return {
        'git_commit_short': sha[:7],
        'git_commit_url': '%s/commit/%s' % (GITHUB_REPO_URL, sha) if sha else '',
    }
