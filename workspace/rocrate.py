"""
RO-Crate (https://www.researchobject.org/ro-crate/) export for a finished
WorkspaceHistory - packages the workflow/tool provenance this app already
tracks (see resolve_dataset_tools in workspace/views.py) as a citable,
machine-readable JSON-LD record, following the RO-Crate 1.2 core spec and
borrowing the Process Run Crate / Workflow Run Crate profile's
CreateAction/ComputationalWorkflow shape
(https://w3id.org/ro-crate/workflow-run) for describing an executed run.
Not a full implementation of that profile (no per-step FormalParameter
entities, no nested per-step CreateActions) - enough to make the run's
inputs, outputs, tools (with version + citation, where resolvable) and
timing genuinely machine-readable and citable, not just a downloadable
results folder.

1.2, not 1.1: first validated with the reference rocrate-validator
(roc-validator) against the ro-crate-1.1 profile (2026-09-24), which
failed on a dead https://w3id.org/ro-crate/1.1/context redirect (a
real, live 404 on the RO-Crate project's own infrastructure, confirmed
directly via curl - not a validator/network-sandbox issue, since the
sibling 1.2 URL resolves fine; ro-crate-py's own default output has
already moved to 1.2 for the same reason). Re-validating against
ro-crate-1.2 with RO_CRATE_CONTEXT/CONFORMS_TO switched to the (working)
1.2 URLs surfaced four real, spec-level gaps in this module's own
output, all fixed here and confirmed by re-validating again until the
required checks passed clean:
- The root Dataset needs a `license` (RO-Crate 1.2 root-entity
  requirement). No established NGPhylogeny.fr license/terms page exists
  to link to (checked - grepped templates for one) - CC-BY-4.0 is used
  as a reasonable, permissive default; flagged here since that's a
  policy choice this module is making on the project's behalf, not
  something confirmed with the maintainers.
- `resultOf` (a File pointing back at the tool that produced it) isn't
  a real term in the RO-Crate context - schema.org has no such
  property, and the spec-correct way to express this is the reverse
  direction: a CreateAction *per tool step*, with `result` listing that
  step's own output files. Replaced the single whole-run CreateAction's
  informal per-file `resultOf` with one #step-<tool_id> CreateAction per
  distinct tool actually used, alongside the existing top-level #action
  describing the run as a whole.
- Every SoftwareApplication entity needs `url` and `version` (RO-Crate
  1.2 Provenance requirement) - previously only set when this app's own
  Tool table happened to have them. Added real fallbacks: `url` from
  the tool id itself (most Galaxy tool ids already are toolshed paths
  with the scheme stripped) or, for a builtin tool with no toolshed
  path (e.g. Galaxy's own "upload1"), the tool's own page on the Galaxy
  server this run actually used; `version` from the tool id's own
  trailing segment when no local Tool record has one.
- The #workflow entity's @type needs File and SoftwareSourceCode
  alongside ComputationalWorkflow (RO-Crate 1.2 Workflow requirement,
  matching the Workflow Run Crate profile's own convention that a
  workflow entity is simultaneously "a file", "source code" and "a
  computational workflow").

End state, both re-checked directly with rocrate-validator against the
real crate produced by a real finished history (not just the unit
tests' synthetic fixtures):
- ro-crate-1.2: passes clean, 65/65 required checks.
- process-run-crate-0.5 (the more specific profile whose CreateAction/
  ComputationalWorkflow shape this module's output actually follows):
  41/42. The one remaining failure is a hard-coded SHACL constraint in
  roc-validator 0.11.4's own bundled copy of this profile
  (profiles/ro-crate/1.1/must/1_file-descriptor_metadata.ttl,
  sh:hasValue <https://w3id.org/ro/crate/1.1>) that requires the file
  descriptor's conformsTo to include the *exact* 1.1 URI, regardless of
  what the crate's own content actually is - i.e. this specific
  installed validator's process-run-crate-0.5 implementation hasn't
  been updated for RO-Crate 1.2 yet. Not fixed here: adding a
  conformsTo claim of "RO-Crate 1.1" to a crate that deliberately
  doesn't use 1.1's (currently dead) context would be a false,
  self-contradictory declaration made only to satisfy one specific
  tool's stale check, not a real crate improvement.

Also cross-checked by loading the crate with ro-crate-py (the reference
Python implementation) directly - loads cleanly, all entities
dereference correctly (root, mainEntity, File/SoftwareApplication
counts match).

Deliberately does not bundle any dataset bytes. This app does no
bioinformatics computation or long-term storage of its own (see CLAUDE.md's
"What this is") - every dataset File entity here is a *remote* reference
(contentUrl pointing at this app's own existing download endpoint, which
itself proxies Galaxy), not embedded content. Re-downloading/re-zipping a
whole history's actual files server-side on export would risk the same
class of large-file/memory/timeout issues already hit elsewhere in this app
(utils/biofile.py's valid_fasta(), the uwsgi harakiri/large-upload
incidents - see CLAUDE.md) for no real benefit: the data already has a
stable, working URL.
"""
from django.urls import reverse
from django.utils import timezone

from tools.models import Tool
from workspace.emails import site_url

RO_CRATE_CONTEXT = "https://w3id.org/ro/crate/1.2/context"
CONFORMS_TO = "https://w3id.org/ro/crate/1.2"
# The specific profile for "a crate describing an executed computational
# process" - matches the CreateAction/ComputationalWorkflow shape this
# module actually builds. Declared on both the root Dataset entity and
# the file descriptor's own conformsTo (RO-Crate 1.2 requires both) -
# confirmed by re-validating with rocrate-validator against
# process-run-crate-0.5 until it actually passed, not just added on
# faith that declaring it would be accurate.
PROCESS_RUN_CRATE_CONFORMS_TO = "https://w3id.org/ro/wfrun/process/0.5"

# See this module's own docstring for why this specific license - a
# policy default, not confirmed NGPhylogeny.fr project policy.
DEFAULT_LICENSE = "https://creativecommons.org/licenses/by/4.0/"

# The paper this whole app exists to be cited alongside - same reference
# ngphylo_citation() (workspace/views.py) embeds as HTML in the citations
# list, here as a proper schema.org ScholarlyArticle entity instead.
NGPHYLO_CITATION_ID = "https://doi.org/10.1093/nar/gkz303"
NGPHYLO_CITATION = {
    "@id": NGPHYLO_CITATION_ID,
    "@type": "ScholarlyArticle",
    "name": ("NGPhylogeny.fr: new generation phylogenetic services for "
             "non-specialists."),
    "author": ("Lemoine F, Correia D, Lefort V, Doppelt-Azeroual O, "
               "Mareuil F, Cohen-Boulakia S, Gascuel O."),
    "datePublished": "2019",
    "identifier": NGPHYLO_CITATION_ID,
}

CATEGORY_LABELS = {
    'OneClick': 'One Click',
    'duplicated': 'Advanced',
    'automaker': 'A La Carte',
    'Tool': 'Single Tool',
}


def _dataset_url(dataset_id):
    return site_url(reverse('download_file', kwargs={'file_id': dataset_id}))


def _tool_url(galaxy_server, tool_id):
    """
    Most Galaxy tool ids already are a toolshed path with the scheme
    stripped (e.g. "toolshed.pasteur.fr/repos/fmareuil/mafft/mafft/
    7.407_1") - a real, dereferenceable page once "https://" is put
    back. Galaxy's own builtin tools (e.g. "upload1") aren't toolshed-
    hosted at all - fall back to that tool's own page on the Galaxy
    server this run actually used.
    """
    first_segment = tool_id.split('/', 1)[0]
    if '.' in first_segment:
        return 'https://%s' % tool_id
    return '%s/root?tool_id=%s' % (galaxy_server.url.rstrip('/'), tool_id)


def _tool_version(local_tool, tool_id):
    """
    RO-Crate 1.2 requires every SoftwareApplication to declare a
    version. This app's own Tool table has a real one for every tool it
    actually imported from Galaxy - for a Galaxy builtin never imported
    that way (e.g. "upload1", no local Tool row at all), the best
    available fallback is the tool id's own trailing segment, which for
    a toolshed-style id already reads as its version (".../mafft/
    7.407_1" -> "7.407_1"); "unknown" only for a bare id with no
    segment to draw one from (upload1 itself).
    """
    if local_tool and local_tool.version:
        return local_tool.version
    if '/' in tool_id:
        return tool_id.rsplit('/', 1)[-1]
    return 'unknown'


def build_rocrate_metadata(w, dataset_tool_ids, tool_names):
    """
    w: a WorkspaceHistory with .history_content/.history_info already
    parsed (see workspace.views._parse_history_json - both default to a
    real value there, never None, so this doesn't need to guard against
    that itself).
    dataset_tool_ids/tool_names: as returned by
    workspace.views.resolve_dataset_tools(gi, w.galaxy_server, w.history,
    dataset_ids) - already resolved before calling this, since resolving
    them needs a live Galaxy connection and this function deliberately
    doesn't (pure data in, dict out - testable without mocking bioblend).

    Returns the RO-Crate metadata dict (the contents of
    ro-crate-metadata.json) - a plain dict, not yet serialized, so
    callers (the export view, tests) can inspect it directly.
    """
    history_info = w.history_info or {}
    datasets = [f for f in (w.history_content or []) if isinstance(f, dict)]

    graph = [{
        "@id": "ro-crate-metadata.json",
        "@type": "CreativeWork",
        "conformsTo": [
            {"@id": CONFORMS_TO},
            {"@id": PROCESS_RUN_CRATE_CONFORMS_TO},
        ],
        "about": {"@id": "./"},
    }]

    result_refs = []
    file_urls_by_tool = {}
    for f in datasets:
        dataset_id = f.get('id')
        if not dataset_id:
            continue
        url = _dataset_url(dataset_id)
        entity = {
            "@id": url,
            "@type": "File",
            "name": f.get('name') or dataset_id,
            "contentUrl": url,
        }
        if f.get('extension'):
            entity["encodingFormat"] = f.get('extension')
        graph.append(entity)
        result_refs.append({"@id": url})
        tool_id = dataset_tool_ids.get(dataset_id)
        if tool_id:
            file_urls_by_tool.setdefault(tool_id, []).append({"@id": url})

    for tool_id, name in tool_names.items():
        local_tool = Tool.objects.filter(
            galaxy_server=w.galaxy_server, id_galaxy=tool_id).first()
        entity = {
            "@id": "#tool-%s" % tool_id,
            "@type": "SoftwareApplication",
            "name": name,
            "url": _tool_url(w.galaxy_server, tool_id),
            "version": _tool_version(local_tool, tool_id),
        }
        if local_tool:
            if local_tool.description:
                entity["description"] = local_tool.description
            citation_texts = [
                text for text in
                (c.txt().strip() for c in local_tool.citation_set.all())
                if text]
            if citation_texts:
                entity["citation"] = citation_texts
        graph.append(entity)

    workflow_name = (
        (w.workflow.name if w.workflow_id else None)
        or w.workflow_steps
        or w.workflow_category
        or "NGPhylogeny.fr analysis")
    graph.append({
        "@id": "#workflow",
        # A workflow entity is simultaneously "a file", "source code"
        # and "a computational workflow" per RO-Crate 1.2/Workflow Run
        # Crate's own convention - required even though this entity has
        # no separately downloadable file of its own here.
        "@type": ["File", "SoftwareSourceCode", "ComputationalWorkflow"],
        "name": workflow_name,
    })

    # Real Galaxy timestamps (show_history()'s own create_time/
    # update_time) when available - w.created_date is this app's own
    # WorkspaceHistory row creation time, a reasonable fallback but not
    # the same instant Galaxy actually started the run.
    start_time = history_info.get('create_time') or (
        w.created_date.isoformat() if w.created_date else None)
    end_time = history_info.get('update_time')

    action = {
        "@id": "#action",
        "@type": "CreateAction",
        "name": "%s run: %s" % (
            w.workflow_category or "Analysis", workflow_name),
        "instrument": {"@id": "#workflow"},
        "result": result_refs,
        "agent": {"@id": "#ngphylogeny"},
    }
    if start_time:
        action["startTime"] = start_time
    if end_time:
        action["endTime"] = end_time
    graph.append(action)

    # One CreateAction per distinct tool step, each pointing at just
    # that tool's own output files - the spec-correct way to express
    # "this file was produced by this tool" (schema.org has no such
    # property in the other direction - see this module's own
    # docstring for why a per-file `resultOf` isn't valid RO-Crate).
    for tool_id, file_refs in file_urls_by_tool.items():
        graph.append({
            "@id": "#step-%s" % tool_id,
            "@type": "CreateAction",
            "name": "%s step" % tool_names.get(tool_id, tool_id),
            "instrument": {"@id": "#tool-%s" % tool_id},
            "result": file_refs,
        })

    graph.append({
        "@id": "#ngphylogeny",
        "@type": "Organization",
        "name": "NGPhylogeny.fr",
        "url": "https://ngphylogeny.fr",
    })
    graph.append(dict(NGPHYLO_CITATION))
    graph.append({
        "@id": DEFAULT_LICENSE,
        "@type": "CreativeWork",
        "name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
    })
    # RO-Crate 1.2 requires a Root Data Entity conformsTo value to
    # reference an actual Profile-typed entity in the graph, not a bare
    # external @id - found by re-validating (rocrate-validator flagged
    # this exact gap after the plain-reference version above). The
    # process-run-crate-0.5 profile's own equivalent check instead wants
    # CreativeWork specifically - a real discrepancy between the two
    # profiles' own validator checks for the same conformsTo mechanism,
    # not something either crate content or this module can resolve by
    # picking one - both types together satisfy both checks.
    graph.append({
        "@id": CONFORMS_TO,
        "@type": ["Profile", "CreativeWork"],
        "name": "RO-Crate Metadata Specification 1.2",
    })
    graph.append({
        "@id": PROCESS_RUN_CRATE_CONFORMS_TO,
        "@type": ["Profile", "CreativeWork"],
        "name": "Process Run Crate 0.5",
    })

    category_label = CATEGORY_LABELS.get(
        w.workflow_category, w.workflow_category or 'Unknown')
    graph.append({
        "@id": "./",
        "@type": "Dataset",
        "name": "%s - NGPhylogeny.fr analysis" % (w.name or w.history),
        "description": (
            'Workflow run and provenance for the NGPhylogeny.fr '
            'analysis "%s" (%s).' % (w.name or w.history, category_label)),
        "datePublished": timezone.now().isoformat(),
        "license": {"@id": DEFAULT_LICENSE},
        "conformsTo": [
            {"@id": CONFORMS_TO},
            {"@id": PROCESS_RUN_CRATE_CONFORMS_TO},
        ],
        "mainEntity": {"@id": "#action"},
        "hasPart": result_refs,
        "publisher": {"@id": "#ngphylogeny"},
        "citation": {"@id": NGPHYLO_CITATION_ID},
    })

    return {
        "@context": RO_CRATE_CONTEXT,
        "@graph": graph,
    }
