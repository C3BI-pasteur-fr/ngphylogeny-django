"""
Contact-form notification email (surveys.views.FeedbackCreateView) - a
plain internal notice to staff that a Feedback row was submitted, sent
in addition to (not instead of) saving that row - see
NGPHYLO_CONTACT_FORM_RECIPIENTS in settings/base.py. Deliberately plain
text, not workspace/emails.py's branded HTML template: this is an
internal operational notice to staff, not a user-facing notification,
and doesn't need the logo/CID MIME wiring built for those.
"""
from django.conf import settings
from django.core.mail import EmailMessage


def build_feedback_notification_email(feedback):
    """
    Returns an unsent EmailMessage for a just-submitted Feedback row -
    see send_feedback_notification_email(), which actually sends it.
    Split out so tests can inspect the built message without needing a
    working SMTP backend/recipients.
    """
    subject = 'NGPhylogeny.fr contact form: %s' % (
        feedback.title or feedback.get_type_display())
    body = (
        'A new message was submitted through the NGPhylogeny.fr contact '
        'form.\n\n'
        'Type: %s\n'
        'From: %s\n'
        'Subject: %s\n\n'
        '%s\n' % (
            feedback.get_type_display(), feedback.email,
            feedback.title, feedback.comment)
    )
    return EmailMessage(
        subject, body, settings.NGPHYLO_REPORT_FROM_EMAIL,
        settings.NGPHYLO_CONTACT_FORM_RECIPIENTS,
        # So a staff member can just hit "reply" to answer the
        # submitter directly, instead of copy-pasting their address out
        # of the body.
        reply_to=[feedback.email])


def send_feedback_notification_email(feedback):
    """
    No-ops if NGPHYLO_CONTACT_FORM_RECIPIENTS isn't set - same
    "safe to leave enabled everywhere" convention as
    workspace.tasks.send_daily_report/NGPHYLO_REPORT_RECIPIENTS.
    """
    if not settings.NGPHYLO_CONTACT_FORM_RECIPIENTS:
        return
    build_feedback_notification_email(feedback).send(fail_silently=False)
