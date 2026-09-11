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
from collections import defaultdict, OrderedDict
from datetime import timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from django.db.models import Count
from django.db.models.functions import TruncWeek
from django.template.loader import render_to_string
from django.utils import timezone

from workspace.models import WorkspaceHistory

# WorkspaceHistory.workflow_category as actually stored (see CLAUDE.md,
# "Workflow duplicates and the Celery cleanup jobs") vs. what it means to a
# reader of this report. 'duplicated' in particular is not a rerun-only
# marker - it's also what the ordinary "Advanced" submission form produces,
# because Workflow.duplicate() (used to give every real run its own Galaxy-
# side workflow copy) unconditionally sets category='duplicated' on the
# copy it returns, and that's what wkadvanced.py's form_valid() then reads
# back to decide the WorkspaceHistory's own workflow_category.
CATEGORY_LABELS = {
    'OneClick': 'OneClick',
    'duplicated': 'Advanced',
    'automaker': 'A La Carte',
    'Tool': 'Single Tool',
}
CATEGORY_COLORS = {
    'OneClick': '#4C72B0',
    'duplicated': '#DD8452',
    'automaker': '#55A868',
    'Tool': '#C44E52',
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


def gather_last_7_days():
    """
    Returns (days, by_day_category, by_day_oneclick_workflow):
    - days: the 7 dates from 6 days ago through today, oldest first.
    - by_day_category: {date: {raw_category: count}}
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

    return days, by_day_category, by_day_oneclick_workflow


def gather_all_time():
    """
    Returns (by_category, by_workflow): {raw_category: count} and
    {workflow/tool name: count}, over every WorkspaceHistory row ever
    created.
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
    return by_category, by_workflow


def gather_weekly_totals():
    """
    Returns a list of (week_start_date, count) tuples, one per ISO week
    (Monday-starting) from the first ever WorkspaceHistory row through the
    current week, oldest first - including weeks with zero runs, so a
    sparse early history reads as a real gap rather than being silently
    compressed out of the chart.
    """
    rows = (
        WorkspaceHistory.objects
        .annotate(week=TruncWeek('created_date'))
        .values('week')
        .annotate(count=Count('id'))
        .order_by('week')
    )
    counts_by_week = {row['week'].date(): row['count'] for row in rows}
    if not counts_by_week:
        return []

    first_week = min(counts_by_week)
    last_week = max(counts_by_week)
    weekly_totals = []
    week = first_week
    while week <= last_week:
        weekly_totals.append((week, counts_by_week.get(week, 0)))
        week += timedelta(weeks=1)
    return weekly_totals


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


def render_weekly_chart(weekly_totals):
    if not weekly_totals:
        return None
    labels = [d.strftime('%Y-%m-%d') for d, _ in weekly_totals]
    values = [n for _, n in weekly_totals]
    # Widen the figure for long histories so weekly bars/labels don't
    # overlap into an unreadable smear.
    fig, ax = plt.subplots(figsize=(max(9, 0.35 * len(weekly_totals)), 4))
    ax.bar(labels, values, color='#4C72B0')
    ax.set_ylabel('Workflows run')
    ax.set_xlabel('Week starting')
    ax.set_title('Workflows per week, since the beginning')
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


def build_report_context():
    """
    Returns (context, images): context is the template context (chart
    slots hold a cid: name - see CID_* above - or None if there was no
    data to chart), images is {cid_name: png_bytes} for only the charts
    that actually got rendered.
    """
    days, by_day_category, by_day_oneclick_workflow = gather_last_7_days()
    by_category, by_workflow = gather_all_time()

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
        CID_WEEKLY: render_weekly_chart(gather_weekly_totals()),
        CID_ALLTIME_CATEGORY: render_alltime_category_pie(by_category),
        CID_ALLTIME_WORKFLOW: render_alltime_workflow_bar(by_workflow),
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
    }
    return context, images


def render_report_html(context):
    return render_to_string('workspace/daily_report_email.html', context)


def render_report_email():
    """
    Returns (html, images) ready for workspace.tasks.send_daily_report to
    send: html references each chart as cid:<name>, images is
    {cid_name: png_bytes} to attach as inline (Content-ID) parts.
    """
    context, images = build_report_context()
    return render_report_html(context), images
