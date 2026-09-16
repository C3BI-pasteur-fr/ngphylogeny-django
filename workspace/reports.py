"""
Data gathering, chart rendering, and HTML rendering for the daily
workflow-usage report emailed by workspace.tasks.send_daily_report.

Counts every WorkspaceHistory row ever created (a row per job actually
submitted, one per Galaxy history) - including deleted=True ones, since
deleted only means the underlying Galaxy history/workflow was purged for
storage (see workspace.tasks.deleteoldgalaxyhistory), not that the run
didn't happen. This is a usage report, not a "what's still retained"
report.

Charts are rendered as PNG bytes and referenced from the HTML via
cid: URIs (see render_report_email()) rather than embedded as base64
data: URIs - plenty of real-world email clients, Outlook chief among
them, simply don't render data: URI images in HTML email at all, leaving
blank space where a chart should be. Content-ID inline attachments are
the long-standing, universally-supported way to put an image in an email.
"""
import base64
import io
from collections import defaultdict, OrderedDict
from datetime import timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from django.core.cache import cache
from django.db.models import Count
from django.db.models.functions import TruncMonth, TruncWeek
from django.template.loader import render_to_string
from django.utils import timezone

from blast.models import BlastRun
from workspace.models import WorkspaceHistory

# WorkspaceHistory.workflow_category as actually stored (see CLAUDE.md,
# "Workflow duplicates and the Celery cleanup jobs") vs. what it means to a
# reader of this report. 'duplicated' in particular is not a rerun-only
# marker - it's also what the ordinary "Advanced" submission form produces,
# because Workflow.duplicate() (used to give every real run its own Galaxy-
# side workflow copy) unconditionally sets category='duplicated' on the
# copy it returns, and that's what wkadvanced.py's form_valid() then reads
# back to decide the WorkspaceHistory's own workflow_category.
#
# 'blast' is not a real workflow_category value - there's no
# WorkspaceHistory row for a BLAST search at all (blast is a parallel,
# self-contained app - see CLAUDE.md's "App responsibilities" - with its
# own BlastRun model/date/deleted fields, not a Workflow/WorkspaceHistory).
# gather_last_7_days()/gather_all_time() merge BlastRun counts into the
# same by_day_category/by_category dicts under this key so BLAST shows up
# as just one more category in the existing charts/tables, without every
# category-consuming function needing to know it's sourced differently.
CATEGORY_LABELS = {
    'OneClick': 'OneClick',
    'duplicated': 'Advanced',
    'automaker': 'A La Carte',
    'Tool': 'Single Tool',
    'blast': 'BLAST',
}
CATEGORY_COLORS = {
    'OneClick': '#4C72B0',
    'duplicated': '#DD8452',
    'automaker': '#55A868',
    'Tool': '#C44E52',
    'blast': '#937860',
}
DEFAULT_COLOR = '#8172B2'


def _category_label(raw):
    return CATEGORY_LABELS.get(raw, raw or 'Unknown')


def _category_color(raw):
    return CATEGORY_COLORS.get(raw, DEFAULT_COLOR)


def _workflow_label(row):
    """
    Best-effort identification of *which* workflow/tool a WorkspaceHistory
    row represents - there's no single field for this:
    - 'Tool' (single-tool) runs never get a Workflow FK (only
      tools/views.py's create_history(..., wf_steps=tool_obj.name) call),
      so the tool name only lives in workflow_steps.
    - 'automaker' (A La Carte) runs are each a unique, one-off tool
      combination with no reusable identity worth naming individually.
    - everything else (OneClick, Advanced) has a real Workflow FK whose
      name is the actual "<Tool> OneClick" identity.
    """
    cat = row['workflow_category']
    if cat == 'Tool':
        return row['workflow_steps'] or 'Unknown tool'
    if cat == 'automaker':
        return 'A La Carte'
    return row['workflow__name'] or row['workflow_steps'] or 'Unknown'


def _fig_to_png_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=110)
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def _gather_blast_by_day(start):
    """
    {date: count} of BlastRun rows (both servers - NCBI and Pasteur -
    lumped into one 'blast' category, same as how OneClick/Advanced/A La
    Carte are each already a single category regardless of which tool
    ran) with date__date >= start. Counts deleted=True rows too, same
    "usage report, not a what's-still-retained report" reasoning as
    WorkspaceHistory - see this module's docstring.
    """
    rows = (
        BlastRun.objects
        .filter(date__date__gte=start)
        .values('date__date')
        .annotate(count=Count('id'))
    )
    return {row['date__date']: row['count'] for row in rows}


def gather_last_7_days():
    """
    Returns (days, by_day_category, by_day_oneclick_workflow):
    - days: the 7 dates from 6 days ago through today, oldest first.
    - by_day_category: {date: {raw_category: count}} - includes a
      synthetic 'blast' category merged in from BlastRun, not just real
      WorkspaceHistory.workflow_category values (see CATEGORY_LABELS).
    - by_day_oneclick_workflow: {date: {workflow_name: count}}, OneClick
      submissions only.
    """
    today = timezone.localdate()
    start = today - timedelta(days=6)
    days = [start + timedelta(days=i) for i in range(7)]

    rows = (
        WorkspaceHistory.objects
        .filter(created_date__date__gte=start)
        .values('created_date__date', 'workflow_category', 'workflow_steps',
                 'workflow__name')
        .annotate(count=Count('id'))
    )

    by_day_category = OrderedDict((d, defaultdict(int)) for d in days)
    by_day_oneclick_workflow = OrderedDict((d, defaultdict(int)) for d in days)
    for row in rows:
        d = row['created_date__date']
        if d not in by_day_category:
            continue
        cat = row['workflow_category']
        by_day_category[d][cat] += row['count']
        if cat == 'OneClick':
            by_day_oneclick_workflow[d][_workflow_label(row)] += row['count']

    for d, count in _gather_blast_by_day(start).items():
        if d in by_day_category:
            by_day_category[d]['blast'] += count

    return days, by_day_category, by_day_oneclick_workflow


def gather_all_time():
    """
    Returns (by_category, by_workflow): {raw_category: count} and
    {workflow/tool name: count}, over every WorkspaceHistory row ever
    created, plus a synthetic 'blast' entry in by_category merged in from
    every BlastRun row ever created (see CATEGORY_LABELS/
    _gather_blast_by_day's docstring). by_workflow has no BLAST
    breakdown - there's no per-workflow identity to a BLAST search the
    way there is for a OneClick tool.
    """
    rows = (
        WorkspaceHistory.objects
        .values('workflow_category', 'workflow_steps', 'workflow__name')
        .annotate(count=Count('id'))
    )
    by_category = defaultdict(int)
    by_workflow = defaultdict(int)
    for row in rows:
        by_category[row['workflow_category']] += row['count']
        by_workflow[_workflow_label(row)] += row['count']
    by_category['blast'] += BlastRun.objects.count()
    return by_category, by_workflow


def gather_blast_query_lengths():
    """
    Returns a list of every known BlastRun.query_length (all-time, both
    servers, deleted=True included - same reasoning as gather_all_time()).
    query_length is only set from BlastRun.date onward it was added
    (blast/tasks.py's launch_ncbi_blast/launch_pasteur_blast, and
    deleteoldblastruns() derives it retroactively from query_seq before
    clearing that field) - excludes NULL rather than treating them as 0,
    since "unknown" and "empty query" aren't the same thing.
    """
    return list(
        BlastRun.objects.exclude(query_length__isnull=True)
        .values_list('query_length', flat=True))


# Above this total history span, gather_period_totals() switches from
# weekly to monthly buckets - see that function's docstring.
WEEKLY_TO_MONTHLY_SPAN_DAYS = 731  # ~2 years


def _add_month(d):
    return d.replace(year=d.year + 1, month=1) if d.month == 12 \
        else d.replace(month=d.month + 1)


def gather_period_totals():
    """
    Returns (granularity, totals): granularity is 'week' or 'month', and
    totals is a list of (period_start_date, count) tuples, oldest first,
    covering the first ever WorkspaceHistory row through the current
    period - including empty periods, so a sparse stretch of history reads
    as a real gap rather than being silently compressed out of the chart.

    Buckets by ISO week (Monday-starting) while the total history span is
    <= WEEKLY_TO_MONTHLY_SPAN_DAYS, by calendar month otherwise. A real
    production history restored from years of usage (see CLAUDE.md) can
    span enough weeks that a literal one-bar-per-week chart becomes
    hundreds of bars wide - squeezed into a normal page/email width
    (.report-chart's max-width: 100%), that crushes the chart's height
    down to an illegible sliver. Rolling up into monthly buckets once the
    span gets that long keeps "since the beginning" readable indefinitely,
    however much real history eventually accumulates.
    """
    first = (WorkspaceHistory.objects.order_by('created_date')
             .values_list('created_date', flat=True).first())
    if first is None:
        return 'week', []

    span_days = (timezone.now() - first).days
    granularity = 'week' if span_days <= WEEKLY_TO_MONTHLY_SPAN_DAYS else 'month'
    trunc = TruncWeek if granularity == 'week' else TruncMonth

    rows = (
        WorkspaceHistory.objects
        .annotate(period=trunc('created_date'))
        .values('period')
        .annotate(count=Count('id'))
        .order_by('period')
    )
    counts_by_period = {row['period'].date(): row['count'] for row in rows}
    if not counts_by_period:
        return granularity, []

    first_period = min(counts_by_period)
    last_period = max(counts_by_period)
    advance = (lambda d: d + timedelta(weeks=1)) if granularity == 'week' \
        else _add_month
    totals = []
    period = first_period
    while period <= last_period:
        totals.append((period, counts_by_period.get(period, 0)))
        period = advance(period)
    return granularity, totals


def render_daily_category_chart(days, by_day_category):
    categories = sorted(
        {c for day in by_day_category.values() for c in day} |
        set(CATEGORY_LABELS),
        key=lambda c: list(CATEGORY_LABELS).index(c) if c in CATEGORY_LABELS else 99)
    if not any(sum(day.values()) for day in by_day_category.values()):
        return None

    x_labels = [d.strftime('%a %m/%d') for d in days]
    fig, ax = plt.subplots(figsize=(8, 4))
    bottom = [0] * len(days)
    for cat in categories:
        values = [by_day_category[d].get(cat, 0) for d in days]
        if not any(values):
            continue
        ax.bar(x_labels, values, bottom=bottom, label=_category_label(cat),
               color=_category_color(cat))
        bottom = [b + v for b, v in zip(bottom, values)]
    ax.set_ylabel('Workflows run')
    ax.set_title('Workflows per day, by category (last 7 days)')
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    return _fig_to_png_bytes(fig)


def render_daily_oneclick_chart(days, by_day_oneclick_workflow):
    workflow_names = sorted(
        {name for day in by_day_oneclick_workflow.values() for name in day})
    if not workflow_names:
        return None

    x_labels = [d.strftime('%a %m/%d') for d in days]
    palette = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(8, 4))
    bottom = [0] * len(days)
    for i, name in enumerate(workflow_names):
        values = [by_day_oneclick_workflow[d].get(name, 0) for d in days]
        ax.bar(x_labels, values, bottom=bottom, label=name,
               color=palette(i % 10))
        bottom = [b + v for b, v in zip(bottom, values)]
    ax.set_ylabel('OneClick workflows run')
    ax.set_title('OneClick workflows per day, by tool (last 7 days)')
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    return _fig_to_png_bytes(fig)


def render_alltime_category_pie(by_category):
    items = sorted(((c, n) for c, n in by_category.items() if n > 0),
                    key=lambda x: -x[1])
    if not items:
        return None
    labels = [_category_label(c) for c, _ in items]
    values = [n for _, n in items]
    colors = [_category_color(c) for c, _ in items]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(values, labels=labels, autopct='%1.0f%%', colors=colors,
           startangle=90)
    ax.set_title('All-time workflows by category')
    return _fig_to_png_bytes(fig)


def render_alltime_workflow_bar(by_workflow, top_n=15):
    items = sorted(by_workflow.items(), key=lambda x: -x[1])[:top_n]
    if not items:
        return None
    items.reverse()  # horizontal bar: largest at the top
    labels = [name for name, _ in items]
    values = [n for _, n in items]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.4 * len(items))))
    ax.barh(labels, values, color='#4C72B0')
    ax.set_xlabel('Total workflows run')
    ax.set_title('All-time workflows by workflow/tool' +
                  (' (top %d)' % top_n if len(by_workflow) > top_n else ''))
    return _fig_to_png_bytes(fig)


def render_blast_query_length_histogram(lengths, bins=30):
    if not lengths:
        return None
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(lengths, bins=min(bins, len(set(lengths))), color='#937860')
    ax.set_xlabel('Query length (bp/aa)')
    ax.set_ylabel('BLAST searches')
    ax.set_title('All-time BLAST query lengths (n=%d)' % len(lengths))
    return _fig_to_png_bytes(fig)


def render_period_chart(granularity, totals):
    if not totals:
        return None
    date_fmt = '%Y-%m-%d' if granularity == 'week' else '%Y-%m'
    labels = [d.strftime(date_fmt) for d, _ in totals]
    values = [n for _, n in totals]
    # Widen the figure for long histories so bars/labels don't overlap
    # into an unreadable smear - gather_period_totals() itself keeps the
    # bar count bounded by switching to monthly buckets for long spans,
    # this just handles however many of either granularity there are.
    fig, ax = plt.subplots(figsize=(max(9, 0.35 * len(totals)), 4))
    ax.bar(labels, values, color='#4C72B0')
    ax.set_ylabel('Workflows run')
    ax.set_xlabel('Week starting' if granularity == 'week' else 'Month')
    ax.set_title('Workflows per %s, since the beginning' % granularity)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    return _fig_to_png_bytes(fig)


# Content-ID names for each chart, referenced from the template as
# cid:<name> and matched up with the inline MIMEImage attachments built
# by workspace.tasks.send_daily_report.
CID_DAILY_CATEGORY = 'chart_daily_category'
CID_DAILY_ONECLICK = 'chart_daily_oneclick'
CID_WEEKLY = 'chart_weekly'
CID_ALLTIME_CATEGORY = 'chart_alltime_category'
CID_ALLTIME_WORKFLOW = 'chart_alltime_workflow'
CID_BLAST_LENGTH_HISTOGRAM = 'chart_blast_length_histogram'


def build_report_context():
    """
    Returns (context, images): context is the template context (chart
    slots hold a cid: name - see CID_* above - or None if there was no
    data to chart), images is {cid_name: png_bytes} for only the charts
    that actually got rendered.
    """
    days, by_day_category, by_day_oneclick_workflow = gather_last_7_days()
    by_category, by_workflow = gather_all_time()
    blast_query_lengths = gather_blast_query_lengths()

    all_categories = sorted(
        {c for day in by_day_category.values() for c in day} |
        set(CATEGORY_LABELS),
        key=lambda c: list(CATEGORY_LABELS).index(c) if c in CATEGORY_LABELS else 99)

    daily_table = []
    for d in days:
        counts = by_day_category[d]
        daily_table.append({
            'date': d,
            'total': sum(counts.values()),
            'per_category': [(_category_label(c), counts.get(c, 0))
                              for c in all_categories],
        })
    category_totals = [
        sum(by_day_category[d].get(c, 0) for d in days)
        for c in all_categories
    ]

    chart_bytes = {
        CID_DAILY_CATEGORY: render_daily_category_chart(days, by_day_category),
        CID_DAILY_ONECLICK: render_daily_oneclick_chart(
            days, by_day_oneclick_workflow),
        CID_WEEKLY: render_period_chart(*gather_period_totals()),
        CID_ALLTIME_CATEGORY: render_alltime_category_pie(by_category),
        CID_ALLTIME_WORKFLOW: render_alltime_workflow_bar(by_workflow),
        CID_BLAST_LENGTH_HISTOGRAM: render_blast_query_length_histogram(
            blast_query_lengths),
    }
    images = {cid: data for cid, data in chart_bytes.items() if data is not None}

    context = {
        'generated_at': timezone.now(),
        'category_columns': [_category_label(c) for c in all_categories],
        'daily_table': daily_table,
        'category_totals': category_totals,
        'week_total': sum(row['total'] for row in daily_table),
        'daily_category_chart': (
            CID_DAILY_CATEGORY if CID_DAILY_CATEGORY in images else None),
        'daily_oneclick_chart': (
            CID_DAILY_ONECLICK if CID_DAILY_ONECLICK in images else None),
        'weekly_chart': CID_WEEKLY if CID_WEEKLY in images else None,
        'alltime_total': sum(by_category.values()),
        'alltime_category_table': sorted(
            ((_category_label(c), n) for c, n in by_category.items()),
            key=lambda x: -x[1]),
        'alltime_category_chart': (
            CID_ALLTIME_CATEGORY if CID_ALLTIME_CATEGORY in images else None),
        'alltime_workflow_table': sorted(
            by_workflow.items(), key=lambda x: -x[1]),
        'alltime_workflow_chart': (
            CID_ALLTIME_WORKFLOW if CID_ALLTIME_WORKFLOW in images else None),
        'blast_length_count': len(blast_query_lengths),
        'blast_length_chart': (
            CID_BLAST_LENGTH_HISTOGRAM
            if CID_BLAST_LENGTH_HISTOGRAM in images else None),
    }
    return context, images


CHART_CONTEXT_KEYS = [
    'daily_category_chart', 'daily_oneclick_chart', 'weekly_chart',
    'alltime_category_chart', 'alltime_workflow_chart', 'blast_length_chart',
]


def render_report_html(context):
    return render_to_string('workspace/daily_report_email.html', context)


def render_report_email():
    """
    Returns (html, images) ready for workspace.tasks.send_daily_report to
    send: html references each chart as cid:<name>, images is
    {cid_name: png_bytes} to attach as inline (Content-ID) parts.
    """
    context, images = build_report_context()
    for key in CHART_CONTEXT_KEYS:
        if context[key]:
            context[key] = 'cid:' + context[key]
    return render_report_html(context), images


REPORT_WEB_CACHE_KEY = 'workspace_daily_report_web_context'
# Rendering 5 matplotlib charts is the slow part of assembling this
# report - cheap enough for the email (a Celery task nobody's waiting on),
# but noticeably slow as something a person loads in a browser and sits
# waiting for. A usage dashboard doesn't need to be second-fresh, so cache
# the assembled context instead of rebuilding it from scratch on every
# single page view.
REPORT_WEB_CACHE_TTL = 15 * 60


def build_report_web_context(force_refresh=False):
    """
    Same report as render_report_email(), but for viewing directly in a
    browser (workspace.views.daily_report_view): chart context values are
    base64 data: URIs, which - unlike in an actual emailed message (see
    render_report_email() and this module's docstring for why that
    deliberately avoids them) - every browser renders just fine.

    Cached for REPORT_WEB_CACHE_TTL seconds - see that constant. Pass
    force_refresh=True to bypass and rebuild immediately (and re-cache
    the fresh result).
    """
    if not force_refresh:
        cached = cache.get(REPORT_WEB_CACHE_KEY)
        if cached is not None:
            return cached

    context, images = build_report_context()
    for key in CHART_CONTEXT_KEYS:
        cid = context[key]
        if cid:
            context[key] = ('data:image/png;base64,' +
                             base64.b64encode(images[cid]).decode('ascii'))
    cache.set(REPORT_WEB_CACHE_KEY, context, REPORT_WEB_CACHE_TTL)
    return context
