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
import io
import json
from collections import defaultdict, OrderedDict
from datetime import datetime, time, timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from django.contrib.auth.models import User
from django.core.cache import cache
from django.db.models import Count
from django.db.models.functions import TruncMonth, TruncWeek
from django.template.loader import render_to_string
from django.urls import reverse
from django.utils import timezone

from blast.models import BlastRun
from workspace.emails import site_url
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

    # created_date__gte=<a plain datetime>, not created_date__date__gte=
    # start: the latter compiles to django_datetime_cast_date(created_date,
    # UTC, UTC) >= start, which wraps the column in a function - Postgres
    # can't use a plain index on created_date for that, and this table
    # is large enough (701,957 rows) that this alone cost ~4s per report
    # (EXPLAIN ANALYZE, 2026-09-23), full sequential scan every time. A
    # plain range comparison against the raw column is index-friendly
    # and, since start is already a UTC-normalized date
    # (timezone.localdate() under this project's TIME_ZONE='UTC'), gives
    # the exact same boundary.
    start_dt = timezone.make_aware(datetime.combine(start, time.min))
    rows = (
        WorkspaceHistory.objects
        .filter(created_date__gte=start_dt)
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


def gather_user_growth():
    """
    Returns (granularity, totals): same weekly/monthly period-bucketing
    convention as gather_period_totals() above (see
    WEEKLY_TO_MONTHLY_SPAN_DAYS - shares the same threshold, so a long-
    lived deployment's user-growth chart degrades to monthly bars for
    the same "stay readable" reason), but totals is a *cumulative*
    running total of registered accounts (auth.User.date_joined) as of
    the end of each period - not new signups per period. "Evolution of
    the number of accounts" is naturally a growth curve (how many
    accounts exist by this point), not a bar-per-period signup count.
    """
    first = (User.objects.order_by('date_joined')
              .values_list('date_joined', flat=True).first())
    if first is None:
        return 'week', []

    span_days = (timezone.now() - first).days
    granularity = 'week' if span_days <= WEEKLY_TO_MONTHLY_SPAN_DAYS else 'month'
    trunc = TruncWeek if granularity == 'week' else TruncMonth

    rows = (
        User.objects
        .annotate(period=trunc('date_joined'))
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
    running_total = 0
    while period <= last_period:
        running_total += counts_by_period.get(period, 0)
        totals.append((period, running_total))
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


def render_user_growth_chart(granularity, totals):
    if not totals:
        return None
    date_fmt = '%Y-%m-%d' if granularity == 'week' else '%Y-%m'
    labels = [d.strftime(date_fmt) for d, _ in totals]
    values = [n for _, n in totals]
    fig, ax = plt.subplots(figsize=(max(9, 0.35 * len(totals)), 4))
    ax.plot(labels, values, color='#55A868', marker='o', markersize=3)
    ax.set_ylabel('Total registered accounts')
    ax.set_xlabel('Week starting' if granularity == 'week' else 'Month')
    ax.set_title('Registered accounts over time (cumulative)')
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
CID_USER_GROWTH = 'chart_user_growth'


def _gather_all():
    """
    Every gather_* call, once - shared by build_report_context() (email)
    and build_report_web_context() (the interactive web page), each of
    which calls this independently (the web side's own 15-minute cache -
    see REPORT_WEB_CACHE_TTL - already keeps that from doubling real
    request-time DB load; this is just avoiding two copies of the same
    seven-function call list drifting out of sync with each other).
    """
    days, by_day_category, by_day_oneclick_workflow = gather_last_7_days()
    by_category, by_workflow = gather_all_time()
    period_granularity, period_totals = gather_period_totals()
    user_growth_granularity, user_growth_totals = gather_user_growth()
    return {
        'days': days,
        'by_day_category': by_day_category,
        'by_day_oneclick_workflow': by_day_oneclick_workflow,
        'by_category': by_category,
        'by_workflow': by_workflow,
        'blast_query_lengths': gather_blast_query_lengths(),
        'period_granularity': period_granularity,
        'period_totals': period_totals,
        'user_growth_granularity': user_growth_granularity,
        'user_growth_totals': user_growth_totals,
    }


def _build_common_context(g):
    """
    The context both build_report_context() (email) and
    build_report_web_context() (web) need - every table/total
    _report_body.html renders that isn't a chart image/canvas. Each
    caller adds its own chart-specific keys (daily_category_chart etc.)
    on top of this.
    """
    days = g['days']
    by_day_category = g['by_day_category']
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

    return {
        'generated_at': timezone.now(),
        'category_columns': [_category_label(c) for c in all_categories],
        'daily_table': daily_table,
        'category_totals': category_totals,
        'week_total': sum(row['total'] for row in daily_table),
        'alltime_total': sum(g['by_category'].values()),
        'alltime_category_table': sorted(
            ((_category_label(c), n) for c, n in g['by_category'].items()),
            key=lambda x: -x[1]),
        'alltime_workflow_table': sorted(
            g['by_workflow'].items(), key=lambda x: -x[1]),
        'blast_length_count': len(g['blast_query_lengths']),
        'total_users': User.objects.count(),
    }


def build_report_context():
    """
    Returns (context, images): context is the template context (chart
    slots hold a cid: name - see CID_* above - or None if there was no
    data to chart), images is {cid_name: png_bytes} for only the charts
    that actually got rendered. Used only by the email
    (render_report_email() below) - matplotlib PNGs via cid: attachments
    are unaffected by the interactive web charts added alongside this
    (build_report_web_context()/build_report_web_chart_data() further
    down) - see this module's own docstring for why email can't run the
    JS those need.
    """
    g = _gather_all()
    context = _build_common_context(g)

    chart_bytes = {
        CID_DAILY_CATEGORY: render_daily_category_chart(
            g['days'], g['by_day_category']),
        CID_DAILY_ONECLICK: render_daily_oneclick_chart(
            g['days'], g['by_day_oneclick_workflow']),
        CID_WEEKLY: render_period_chart(
            g['period_granularity'], g['period_totals']),
        CID_ALLTIME_CATEGORY: render_alltime_category_pie(g['by_category']),
        CID_ALLTIME_WORKFLOW: render_alltime_workflow_bar(g['by_workflow']),
        CID_BLAST_LENGTH_HISTOGRAM: render_blast_query_length_histogram(
            g['blast_query_lengths']),
        CID_USER_GROWTH: render_user_growth_chart(
            g['user_growth_granularity'], g['user_growth_totals']),
    }
    images = {cid: data for cid, data in chart_bytes.items() if data is not None}

    context.update({
        'chart_mode': 'email',
        # The email can't offer real hover-interactive charts itself (no
        # JS in an email client - see this module's own docstring), so
        # it points at the one place that can: the live web report
        # (workspace.views.daily_report_view), which renders the same
        # data as Chart.js canvases instead of these static PNGs.
        'report_web_url': site_url(reverse('daily_report')),
        'daily_category_chart': (
            CID_DAILY_CATEGORY if CID_DAILY_CATEGORY in images else None),
        'daily_oneclick_chart': (
            CID_DAILY_ONECLICK if CID_DAILY_ONECLICK in images else None),
        'weekly_chart': CID_WEEKLY if CID_WEEKLY in images else None,
        'alltime_category_chart': (
            CID_ALLTIME_CATEGORY if CID_ALLTIME_CATEGORY in images else None),
        'alltime_workflow_chart': (
            CID_ALLTIME_WORKFLOW if CID_ALLTIME_WORKFLOW in images else None),
        'blast_length_chart': (
            CID_BLAST_LENGTH_HISTOGRAM
            if CID_BLAST_LENGTH_HISTOGRAM in images else None),
        'user_growth_chart': (
            CID_USER_GROWTH if CID_USER_GROWTH in images else None),
    })
    return context, images


CHART_CONTEXT_KEYS = [
    'daily_category_chart', 'daily_oneclick_chart', 'weekly_chart',
    'alltime_category_chart', 'alltime_workflow_chart', 'blast_length_chart',
    'user_growth_chart',
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


# Validated via the dataviz skill's color validator (CVD-safe adjacent
# pairs, >=15 normal-vision floor - run `node scripts/validate_palette.js
# "#2a78d6,#eb6834,#1baf7a,#eda100,#e87ba4" --mode light` from that
# skill's own directory to reproduce) - used only by the interactive web
# charts below (build_report_web_chart_data()). The email's own
# matplotlib charts keep their existing CATEGORY_COLORS (a different,
# non-validated legacy palette - confirmed failing the same validator,
# see CLAUDE.md) unchanged; re-coloring the email's charts was out of
# scope for this change. This site has no dark mode (checked - no
# prefers-color-scheme/data-theme anywhere in assets/css or templates),
# so there's no second, dark-surface-validated set to carry here either.
WEB_CATEGORICAL_PALETTE = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
# Single-series charts (a ranked bar list, a time series, a histogram)
# are a magnitude job, not an identity one - one hue, not a rainbow (see
# the dataviz skill's choosing-a-form.md) - the same palette's slot 1.
WEB_SEQUENTIAL_COLOR = WEB_CATEGORICAL_PALETTE[0]


def _web_slot_color(index):
    return WEB_CATEGORICAL_PALETTE[index % len(WEB_CATEGORICAL_PALETTE)]


def _histogram_bins(values, max_bins=30):
    """
    Plain-Python equal-width histogram (min..max split into equal-width
    buckets, matching matplotlib's own ax.hist default behavior) for
    build_report_web_chart_data()'s BLAST query-length distribution -
    Chart.js has no built-in histogram type, and this avoids depending on
    numpy being importable at runtime for one bucket count (it's only
    ever present here as biopython's own build-time dependency - see
    requirement.txt - not a real runtime dependency of this app).
    """
    lo, hi = min(values), max(values)
    if lo == hi:
        return ['%d' % lo], [len(values)]
    nbins = min(max_bins, len(set(values)))
    width = (hi - lo) / nbins
    counts = [0] * nbins
    for v in values:
        counts[min(int((v - lo) / width), nbins - 1)] += 1
    labels = ['%d–%d' % (lo + i * width, lo + (i + 1) * width)
              for i in range(nbins)]
    return labels, counts


def build_report_web_chart_data(g):
    """
    JSON-serializable chart configs for the interactive web report
    (workspace.views.daily_report_view / templates/workspace/
    _report_charts_init.html) - Chart.js renders these client-side with
    real hover tooltips, a crosshair-style "every series at this X"
    readout, and a legend, per the dataviz skill's interaction.md. Not
    used by the email, which keeps rendering the matplotlib PNGs via
    render_report_email()/build_report_context() entirely unchanged -
    see this module's own docstring for why (email clients don't run
    the JS this needs).
    """
    days = g['days']
    by_day_category = g['by_day_category']
    all_categories = sorted(
        {c for day in by_day_category.values() for c in day} |
        set(CATEGORY_LABELS),
        key=lambda c: list(CATEGORY_LABELS).index(c) if c in CATEGORY_LABELS else 99)
    day_labels = [d.strftime('%a %m/%d') for d in days]

    daily_category = {
        'labels': day_labels,
        'datasets': [
            {'label': _category_label(c), 'color': _web_slot_color(i),
             'data': [by_day_category[d].get(c, 0) for d in days]}
            for i, c in enumerate(all_categories)
            if any(by_day_category[d].get(c, 0) for d in days)
        ],
    }

    by_day_oneclick_workflow = g['by_day_oneclick_workflow']
    workflow_names = sorted(
        {name for day in by_day_oneclick_workflow.values() for name in day})
    daily_oneclick = {
        'labels': day_labels,
        'datasets': [
            {'label': name, 'color': _web_slot_color(i),
             'data': [by_day_oneclick_workflow[d].get(name, 0) for d in days]}
            for i, name in enumerate(workflow_names)
        ],
    }

    period_date_fmt = '%Y-%m-%d' if g['period_granularity'] == 'week' else '%Y-%m'
    weekly = {
        'labels': [d.strftime(period_date_fmt) for d, _ in g['period_totals']],
        'data': [n for _, n in g['period_totals']],
        'color': WEB_SEQUENTIAL_COLOR,
    }

    category_items = sorted(
        ((c, n) for c, n in g['by_category'].items() if n > 0),
        key=lambda x: -x[1])
    alltime_category = {
        'labels': [_category_label(c) for c, _ in category_items],
        'data': [n for _, n in category_items],
        'colors': [_web_slot_color(i) for i in range(len(category_items))],
    }

    workflow_items = sorted(g['by_workflow'].items(), key=lambda x: -x[1])[:15]
    alltime_workflow = {
        'labels': [name for name, _ in workflow_items],
        'data': [n for _, n in workflow_items],
        'color': WEB_SEQUENTIAL_COLOR,
    }

    lengths = g['blast_query_lengths']
    hist_labels, hist_counts = _histogram_bins(lengths) if lengths else ([], [])
    blast_histogram = {
        'labels': hist_labels, 'data': hist_counts, 'color': WEB_SEQUENTIAL_COLOR,
    }

    growth_date_fmt = ('%Y-%m-%d' if g['user_growth_granularity'] == 'week'
                        else '%Y-%m')
    user_growth = {
        'labels': [d.strftime(growth_date_fmt)
                   for d, _ in g['user_growth_totals']],
        'data': [n for _, n in g['user_growth_totals']],
        'color': WEB_SEQUENTIAL_COLOR,
    }

    return {
        'dailyCategory': daily_category,
        'dailyOneclick': daily_oneclick,
        'weekly': weekly,
        'alltimeCategory': alltime_category,
        'alltimeWorkflow': alltime_workflow,
        'blastHistogram': blast_histogram,
        'userGrowth': user_growth,
    }


REPORT_WEB_CACHE_KEY = 'workspace_daily_report_web_context'
# Assembling this report (the gather_* queries above) is genuinely slow
# at real production scale - see WorkspaceHistory.created_date's own
# db_index=True comment for measured numbers - cheap enough for the
# email (a Celery task nobody's waiting on), but noticeably slow as
# something a person loads in a browser and sits waiting for. A usage
# dashboard doesn't need to be second-fresh, so cache the assembled
# context instead of rebuilding it from scratch on every single page
# view.
REPORT_WEB_CACHE_TTL = 15 * 60


def build_report_web_context(force_refresh=False):
    """
    Context for workspace.views.daily_report_view - same tables/totals
    as the email (_build_common_context()), but interactive Chart.js
    charts instead of the email's static matplotlib PNGs (see
    build_report_web_chart_data()'s own docstring for why those have to
    differ). Doesn't call build_report_context()/render matplotlib
    figures at all - a plain dict of numbers for the browser's own JS to
    chart is far cheaper to assemble than rendering 7 PNG figures, on
    top of being the whole point of this change.

    Cached for REPORT_WEB_CACHE_TTL seconds - see that constant. Pass
    force_refresh=True to bypass and rebuild immediately (and re-cache
    the fresh result).
    """
    if not force_refresh:
        cached = cache.get(REPORT_WEB_CACHE_KEY)
        if cached is not None:
            return cached

    g = _gather_all()
    context = _build_common_context(g)
    chart_data = build_report_web_chart_data(g)
    context.update({
        'chart_mode': 'web',
        'chart_data_json': json.dumps(chart_data),
        'daily_category_chart': bool(chart_data['dailyCategory']['datasets']),
        'daily_oneclick_chart': bool(chart_data['dailyOneclick']['datasets']),
        'weekly_chart': bool(chart_data['weekly']['data']),
        'alltime_category_chart': bool(chart_data['alltimeCategory']['data']),
        'alltime_workflow_chart': bool(chart_data['alltimeWorkflow']['data']),
        'blast_length_chart': bool(chart_data['blastHistogram']['data']),
        'user_growth_chart': bool(chart_data['userGrowth']['data']),
    })
    cache.set(REPORT_WEB_CACHE_KEY, context, REPORT_WEB_CACHE_TTL)
    return context
