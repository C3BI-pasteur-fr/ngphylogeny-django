#!/bin/bash
# Deletes workflows from a Galaxy server that are older than $DAYS days and
# are NOT one of the canonical "<Tool> OneClick" base workflows.
#
# Why this is needed: every real OneClick/A La Carte submission creates its
# own Galaxy-side copy of the workflow it launches (Workflow.duplicate() in
# workflows/models.py) - same name as the base workflow, fresh Galaxy id.
# Nothing has ever bulk-cleaned these up once they've actually been run
# (the existing Celery periodic task, workflows.tasks.deleteoldgalaxyworkflows,
# only removes ones that were created but never run at all), so real usage
# accumulates them forever. Galaxy has no field distinguishing "the base
# workflow" from "a per-run duplicate" - they share the exact same name -
# so the only reliable way to know which ids to KEEP is NGPhylogeny's own
# Django DB (category='base'). Those ids also aren't stable long-term
# (Galaxy re-imports its bundled workflows with a fresh id on every
# restart), so this fetches the current keep-list fresh from Postgres on
# every run rather than hardcoding it.
#
# Runs as a bounded, resumable batch loop instead of listing everything at
# once - repeatedly asks Galaxy for one page of its OLDEST workflows,
# oldest-first, and stops as soon as it sees one that's too recent to
# delete (everything after it, being even newer, is guaranteed too recent
# too). Deleting shrinks the list, so re-querying at the same offset=0
# each time naturally surfaces the next-oldest batch instead of skipping
# items - no index-shifting bugs from paginating over a list you're
# concurrently deleting from.
#
# Requires: curl, jq, and (unless KEEP_IDS is set manually) kubectl access
# to run psql against the app's Postgres pod.
#
# Usage:
#   GALAXY_URL=https://galaxy.pasteur.fr GALAXY_KEY=<admin api key> \
#   NAMESPACE=ngphylogenyfr-dev \
#   ./cleanup_old_galaxy_workflows.sh
#
# Defaults to a dry run (lists what WOULD be deleted, deletes nothing).
# Set DRY_RUN=false to actually delete - this is irreversible.
set -euo pipefail

export GALAXY_URL="${GALAXY_URL:?Set GALAXY_URL, e.g. https://galaxy.pasteur.fr}"
export GALAXY_KEY="${GALAXY_KEY:?Set GALAXY_KEY to a Galaxy admin API key}"
DAYS="${DAYS:-7}"
BATCH_SIZE="${BATCH_SIZE:-200}"
DRY_RUN="${DRY_RUN:-true}"
NAMESPACE="${NAMESPACE:-}"
# How many DELETE requests to fire concurrently per batch (xargs -P). At
# hundreds of thousands of rows, one-at-a-time (the old
# SLEEP_BETWEEN_DELETES-throttled sequential loop) would take unreasonably
# long - Galaxy's REST API has no bulk-delete endpoint for workflows (only
# DELETE /api/workflows/{id}, one at a time - checked directly against
# galaxyproject/galaxy's own api/workflows.py), so this parallelizes plain
# per-id DELETE calls client-side instead.
#
# There's no way to know from here what galaxy.pasteur.fr can actually
# handle (its own worker count, any rate-limiting proxy in front of it,
# how much real concurrent user traffic it's carrying) - so instead of
# guessing a "safe" number, the circuit breaker below (FAILURE_*) makes it
# safe to try raising this: a rising error rate aborts the whole script
# rather than continuing to hammer a server that's visibly struggling.
# Start conservative, watch for WARNINGs, raise gradually.
PARALLELISM="${PARALLELISM:-8}"
# Abort if a batch's failure rate looks like Galaxy is struggling, rather
# than plowing on regardless. Only evaluated once a batch has a
# meaningful sample size (FAILURE_MIN_SAMPLE) - a single failure out of 2
# attempts shouldn't trip this, a single failure out of 50 also shouldn't,
# but 40 failures out of 50 clearly should.
FAILURE_ABORT_PERCENT="${FAILURE_ABORT_PERCENT:-30}"
FAILURE_MIN_SAMPLE="${FAILURE_MIN_SAMPLE:-10}"

for cmd in curl jq; do
    command -v "$cmd" >/dev/null || { echo "Missing dependency: $cmd" >&2; exit 1; }
done

# --- Step 1: the keep-list (base OneClick workflows' current Galaxy ids) ---
if [ -n "${KEEP_IDS:-}" ]; then
    echo "Using manually-provided KEEP_IDS (skipping DB lookup)."
    keep_ids="$KEEP_IDS"
else
    [ -n "$NAMESPACE" ] || { echo "Set NAMESPACE (or provide KEEP_IDS manually) so the keep-list can be fetched from Postgres." >&2; exit 1; }
    command -v kubectl >/dev/null || { echo "Missing dependency: kubectl (or set KEEP_IDS manually)" >&2; exit 1; }
    pg_pod=$(kubectl get pods -n "$NAMESPACE" -l app=postgres -o jsonpath='{.items[0].metadata.name}')
    keep_ids=$(kubectl exec -i "$pg_pod" -n "$NAMESPACE" -- \
        psql -U ngphylo -d ngphylo -t -A -c \
        "SELECT id_galaxy FROM workflows_workflow WHERE category='base';")
fi

if [ -n "$keep_ids" ]; then
    keep_ids_json=$(printf '%s\n' $keep_ids | jq -R . | jq -s .)
else
    keep_ids_json='[]'
fi
echo "Keeping $(echo "$keep_ids_json" | jq 'length') base workflow(s): $(echo "$keep_ids_json" | jq -c .)"

cutoff_epoch=$(( $(date -u +%s) - DAYS * 86400 ))
echo "Deleting workflows older than $DAYS days (before $(date -u -r "$cutoff_epoch" '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || date -u -d "@$cutoff_epoch" '+%Y-%m-%dT%H:%M:%SZ'))."
[ "$DRY_RUN" = "true" ] && echo "DRY_RUN=true - nothing will actually be deleted. Set DRY_RUN=false to really delete."

deleted_count=0
examined_count=0

while true; do
    page=$(curl -sS -H "x-api-key: $GALAXY_KEY" \
        "$GALAXY_URL/api/workflows?limit=$BATCH_SIZE&offset=0&sort_by=create_time&sort_desc=false")

    page_len=$(echo "$page" | jq 'length')
    if [ "$page_len" -eq 0 ]; then
        echo "No more workflows to examine. Done."
        break
    fi

    stop=false
    ids_to_delete=()

    # Phase 1 (sequential): scan this page oldest-first, decide what to
    # keep/delete/stop-at. Cheap (local jq parsing, no network calls) so
    # there's no benefit to parallelizing it, and the early-stop logic
    # depends on processing strictly in oldest-first order.
    for i in $(seq 0 $((page_len - 1))); do
        item=$(echo "$page" | jq ".[$i]")
        id=$(echo "$item" | jq -r '.id')
        name=$(echo "$item" | jq -r '.name')
        create_time=$(echo "$item" | jq -r '.create_time')
        examined_count=$((examined_count + 1))

        # Strip fractional seconds and normalize to a strict ISO8601 "Z"
        # timestamp - jq's fromdateiso8601 needs that exact shape, and
        # Galaxy's create_time doesn't include a trailing Z.
        item_epoch=$(echo "$create_time" | sed -E 's/\.[0-9]+$//' | \
            jq -R '. + "Z" | fromdateiso8601')

        is_kept=$(echo "$keep_ids_json" | jq --arg id "$id" 'any(.[]; . == $id)')

        if [ "$is_kept" = "true" ]; then
            echo "SKIP (base workflow): $id  $name"
            continue
        fi

        if [ "$item_epoch" -gt "$cutoff_epoch" ]; then
            echo "Reached a workflow newer than the cutoff ($name, $create_time) - stopping (oldest-first order)."
            stop=true
            break
        fi

        if [ "$DRY_RUN" = "true" ]; then
            echo "WOULD DELETE: $id  $name  ($create_time)"
        else
            echo "Queued for deletion: $id  $name  ($create_time)"
            ids_to_delete+=("$id")
        fi
    done

    # Phase 2 (parallel, real runs only): fire the actual DELETE calls
    # concurrently instead of one at a time. `{}` isn't interpolated
    # directly into the bash -c string - it's passed as $1 - to avoid any
    # shell-injection risk from an id containing shell metacharacters.
    # Each failure is appended to fail_log (one id per line) so the
    # circuit breaker below can measure this batch's actual failure rate -
    # xargs's own exit status alone can't tell us HOW MANY failed.
    if [ "$DRY_RUN" != "true" ] && [ "${#ids_to_delete[@]}" -gt 0 ]; then
        fail_log=$(mktemp)
        export FAIL_LOG="$fail_log"
        printf '%s\n' "${ids_to_delete[@]}" | xargs -P "$PARALLELISM" -I{} \
            bash -c '
                id="$1"
                http_code=$(curl -sS -o /dev/null -w "%{http_code}" -X DELETE \
                    -H "x-api-key: $GALAXY_KEY" \
                    "$GALAXY_URL/api/workflows/$id")
                if [ "$http_code" != "200" ]; then
                    echo "  WARNING: delete failed (HTTP $http_code) for $id - leaving it, will retry on next run." >&2
                    echo "$id" >> "$FAIL_LOG"
                fi
            ' _ {}

        failed_in_page=$(wc -l < "$fail_log" | tr -d ' ')
        rm -f "$fail_log"

        attempted_in_page="${#ids_to_delete[@]}"
        if [ "$attempted_in_page" -ge "$FAILURE_MIN_SAMPLE" ]; then
            failure_percent=$(( failed_in_page * 100 / attempted_in_page ))
            if [ "$failure_percent" -ge "$FAILURE_ABORT_PERCENT" ]; then
                echo "ABORTING: $failed_in_page/$attempted_in_page deletes failed in this batch (${failure_percent}% >= ${FAILURE_ABORT_PERCENT}% threshold) - Galaxy looks like it's struggling under PARALLELISM=$PARALLELISM. Lower PARALLELISM and re-run (already-deleted workflows won't be re-processed)." >&2
                exit 1
            fi
        fi
    fi

    deleted_in_page="${#ids_to_delete[@]}"
    deleted_count=$((deleted_count + deleted_in_page))

    [ "$stop" = "true" ] && break

    if [ "$deleted_in_page" -eq 0 ]; then
        # Nothing in this page was deleted (all kept or - shouldn't reach
        # here since a too-new item stops the loop above - but as a
        # belt-and-suspenders guard against an infinite loop): treat a
        # fully-kept page as the end, there's nothing left this script can
        # act on.
        echo "Nothing further to delete. Done."
        break
    fi
done

# deleted_count is "attempted", not "confirmed successful" - deletes run
# in parallel subshells, so individual failures (logged as WARNINGs above,
# to stderr) aren't fed back into this counter. Failed ones are simply
# left alone on Galaxy and get retried the next time this script runs.
echo "Examined $examined_count workflow(s), $([ "$DRY_RUN" = "true" ] && echo "would have deleted" || echo "attempted to delete") $deleted_count."
