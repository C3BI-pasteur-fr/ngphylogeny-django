"""
HTML job-completion email (see workspace.tasks.updateworkspacestatus) -
NGPhylogeny.fr / Institut Pasteur branded, with inline (Content-ID) logo
images rather than embedded base64 data: URIs - most mail clients,
Outlook chief among them, don't render data: URI images in HTML email at
all (see workspace/reports.py's module docstring - the daily report hit
this exact problem first).
"""
import os
from email.mime.image import MIMEImage

from django.conf import settings
from django.contrib.staticfiles import finders
from django.core.mail import EmailMultiAlternatives
from django.template.loader import render_to_string
from django.urls import reverse

CID_HEADER_LOGO = 'ngphylogeny_logo'
CID_FOOTER_LOGO = 'institut_pasteur_logo'

# Static (STATICFILES_DIRS) paths for each - see templates/base.html's own
# footer for the same Institut Pasteur logo used there.
_LOGO_STATIC_PATHS = {
    CID_HEADER_LOGO: 'images/logo_phylogeny_small.png',
    CID_FOOTER_LOGO: 'images/logo_institut_pasteur.png',
}


def _read_logo(cid):
    path = finders.find(_LOGO_STATIC_PATHS[cid])
    if not path:
        return None
    with open(path, 'rb') as f:
        return f.read()


def _site_url(path):
    """
    Absolute link back to the site for the given path - https if
    NGPHYLO_HTTPS_HOST is configured (the real deployment always sets
    this - see settings/prod.py), http otherwise (local dev, where
    there's no TLS-terminating proxy in front at all).
    """
    https_host = getattr(settings, 'NGPHYLO_HTTPS_HOST', None)
    if https_host:
        return 'https://%s%s' % (https_host, path)
    host = os.environ.get('NGPHYLO_HOST', 'ngphylogeny.fr')
    return 'http://%s%s' % (host, path)


def build_job_completion_email(history_id, recipient, error):
    """
    Returns an unsent EmailMultiAlternatives for a finished job - see
    send_job_completion_email(), which actually sends it. Split out so
    tests can inspect the built message without needing a working SMTP
    backend/recipient.
    """
    results_url = _site_url(
        reverse('history_detail', kwargs={'history_id': history_id}))
    site_home_url = _site_url('/')
    html = render_to_string('workspace/job_completion_email.html', {
        'error': error,
        'results_url': results_url,
        'site_home_url': site_home_url,
        'header_logo_cid': CID_HEADER_LOGO,
        'footer_logo_cid': CID_FOOTER_LOGO,
    })
    if error:
        subject = 'NGPhylogeny.fr - your analysis finished with errors'
        status_line = 'finished with errors'
    else:
        subject = 'NGPhylogeny.fr - your analysis has finished'
        status_line = 'finished successfully'
    plain_text = (
        'Dear NGPhylogeny.fr user,\n\n'
        'Your analysis has %s.\n\n'
        'View your results: %s\n\n'
        'Thank you for using NGPhylogeny.fr.\n'
        'The NGPhylogeny.fr team\n' % (status_line, results_url)
    )

    msg = EmailMultiAlternatives(
        # Reuses the daily report's sender setting rather than a separate
        # hardcoded address: it's the same SMTP account either way, and
        # institutional mail servers (Exchange/O365-style "Send As"
        # restrictions) commonly only let an account send as its own
        # authorized address - see NGPHYLO_REPORT_FROM_EMAIL in
        # settings/base.py.
        subject, plain_text, settings.NGPHYLO_REPORT_FROM_EMAIL, [recipient])
    msg.attach_alternative(html, 'text/html')
    # multipart/related, not multipart/mixed - required for the logos to
    # be treated as inline (cid:-referenced) parts rather than ordinary
    # file attachments alongside the message.
    msg.mixed_subtype = 'related'
    for cid in (CID_HEADER_LOGO, CID_FOOTER_LOGO):
        data = _read_logo(cid)
        if data is None:
            continue
        image = MIMEImage(data, 'png')
        image.add_header('Content-ID', '<%s>' % cid)
        image.add_header('Content-Disposition', 'inline',
                          filename='%s.png' % cid)
        msg.attach(image)
    return msg


def send_job_completion_email(history_id, recipient, error):
    build_job_completion_email(history_id, recipient, error).send(
        fail_silently=False)
