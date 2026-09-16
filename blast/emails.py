"""
HTML BLAST-completion email (see blast.tasks.launch_ncbi_blast/
launch_pasteur_blast/checkblastruns) - reuses workspace.emails' shared
branding (logos, MIME wiring) and its job_completion_email.html template
directly, so a BLAST notification looks the same as the workflow
job-completion email rather than drifting into its own style.
"""
from django.template.loader import render_to_string
from django.urls import reverse

from workspace.emails import (
    CID_FOOTER_LOGO, CID_HEADER_LOGO, build_branded_html_email, site_url)


def build_blast_completion_email(blastrun, recipient):
    """
    Returns an unsent EmailMultiAlternatives for a finished (or errored)
    BlastRun - see send_blast_completion_email(), which actually sends
    it. Split out so tests can inspect the built message without needing
    a working SMTP backend/recipient.
    """
    error = blastrun.status != blastrun.FINISHED
    results_url = site_url(
        reverse('blast_view', kwargs={'pk': blastrun.id}))
    site_home_url = site_url('/')
    html = render_to_string('workspace/job_completion_email.html', {
        'error': error,
        'results_url': results_url,
        'site_home_url': site_home_url,
        'header_logo_cid': CID_HEADER_LOGO,
        'footer_logo_cid': CID_FOOTER_LOGO,
        'success_message': (
            'Your BLAST search has finished successfully. Your results '
            'are ready to view.'),
        'error_message': (
            'Your BLAST search finished with errors. You can inspect '
            'the details from your results page.'),
    })
    if error:
        subject = 'NGPhylogeny.fr - your BLAST search finished with errors'
        status_line = 'finished with errors'
    else:
        subject = 'NGPhylogeny.fr - your BLAST search has finished'
        status_line = 'finished successfully'
    plain_text = (
        'Dear NGPhylogeny.fr user,\n\n'
        'Your BLAST search has %s.\n\n'
        'View your results: %s\n\n'
        'Thank you for using NGPhylogeny.fr.\n'
        'The NGPhylogeny.fr team\n' % (status_line, results_url)
    )
    return build_branded_html_email(subject, plain_text, html, recipient)


def send_blast_completion_email(blastrun, recipient):
    build_blast_completion_email(blastrun, recipient).send(
        fail_silently=False)
