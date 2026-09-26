import markdown2
from django import template
from django.conf import settings
from django.utils.safestring import mark_safe

register = template.Library()


def _render_markdown(text):
    style = getattr(settings, 'MARKDOWN_STYLES', {}).get('default', {})
    extras = list(style.get('extras', {}).keys())
    safe_mode = style.get('safe_mode', False)
    html = markdown2.markdown(text or '', extras=extras, safe_mode=safe_mode)
    return mark_safe(html)


@register.filter(name='markdown')
def markdown_filter(text):
    return _render_markdown(text)


class MarkdownNode(template.Node):
    def __init__(self, nodelist):
        self.nodelist = nodelist

    def render(self, context):
        return _render_markdown(self.nodelist.render(context))


@register.tag(name='markdown')
def markdown_tag(parser, token):
    nodelist = parser.parse(('endmarkdown',))
    parser.delete_first_token()
    return MarkdownNode(nodelist)
