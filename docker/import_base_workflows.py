#!/usr/bin/env python3
"""
One-shot: imports the 4 base "<Tool> OneClick" workflows into a fresh
Galaxy instance, named to match what tools.management.commands.importworkflows
(run by docker/init.sh) picks up, and what real production actually
contains - see docker-compose.standalone.yml, where this runs as the
galaxy-import-workflows service.

Idempotent: skips any workflow whose name is already present on the
target Galaxy, so re-running docker compose up doesn't create duplicates -
see CLAUDE.md's "upgrade vs the old master branch" section for why that
distinction matters (Galaxy accumulates per-run duplicates with the same
name over time; only the count of *base*, never-yet-imported names matters
here).
"""
import json
import os
import sys

from bioblend.galaxy import GalaxyInstance

GALAXY_URL = os.environ['GALAXY_URL']
GALAXY_KEY = os.environ['GALAXY_KEY']
WORKFLOWS_DIR = os.environ.get('WORKFLOWS_DIR', '/galaxytools/workflows')

FILES = {
    'PhyML': 'Galaxy-Workflow-PhyML.ga',
    'PhyML-SMS': 'Galaxy-Workflow-PhyML-SMS.ga',
    'FastTree': 'Galaxy-Workflow-FastTree.ga',
    'FastME': 'Galaxy-Workflow-FastME.ga',
}

gi = GalaxyInstance(url=GALAXY_URL, key=GALAXY_KEY)
existing_names = {w['name'] for w in gi.workflows.get_workflows()}

for tool_name, fname in FILES.items():
    wf_name = tool_name + ' OneClick'
    if wf_name in existing_names:
        print("skipping %s: already imported" % wf_name)
        continue
    path = os.path.join(WORKFLOWS_DIR, fname)
    if not os.path.exists(path):
        print("ERROR: %s not found (WORKFLOWS_DIR=%s)" % (path, WORKFLOWS_DIR),
              file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        wf = json.load(f)
    wf['name'] = wf_name
    imported = gi.workflows.import_workflow_dict(wf)
    print("imported %s: %s" % (wf_name, imported['id']))
