#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
from collections import Counter

ROOT = Path.cwd()

META = ROOT / "00_metadata" / "Step23_manuscript_update"
RESULT = ROOT / "04_results" / "manuscript_update" / "Step23"

QUEUE = META / "STEP23C_QA2E_REMAINING_REVIEW_QUEUE.csv"
OLD_LOCK = META / "STEP23C_QA2E_DECISION_LOCK.csv"

DOCX = RESULT / "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"

OUT_CSV = META / "STEP23C_QA2G_REMAINING_DECISION_LOCK.csv"
OUT_TXT = RESULT / "STEP23C_QA2G_REMAINING_DECISION_DOSSIER.txt"
OUT_SHA = META / "STEP23C_QA2G_LOCK_CHECKSUMS.sha256"

EXPECTED_DOCX_SHA = (
    "0e91d04f82c3a3353bc6214b94312f72767f7da50778d40d0f33c016bc65f598"
)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ============================================================
# HUMAN-APPROVED DECISIONS FROM QA2F REVIEW
#
# IMPORTANT:
# - This file freezes decisions only.
# - It DOES NOT edit the manuscript.
# - CONTEXTUAL means QA2B must inspect surrounding sentence before edit.
# - KEEP means no manuscript change for this item.
# ============================================================

D = {

"QA2C-001": {
    "decision": "REWRITE",
    "edit_group": "TITLE",
    "approved_text":
    "Single-cell analysis of two ESCC discovery cohorts identifies a reproducible treatment-associated fibroblast CDKN1B–GSN CORE2 response",
},

"QA2C-002": {
    "decision": "REWRITE",
    "edit_group": "RUNNING_TITLE",
    "approved_text":
    "Running title: Fibroblast CORE2 response in two ESCC discovery cohorts",
},

"QA2C-003": {
    "decision": "REWRITE",
    "edit_group": "KEYWORDS",
    "approved_text":
    "Keywords: Esophageal squamous cell carcinoma; neoadjuvant therapy; single-cell RNA-seq; fibroblast; tumor microenvironment; cross-cohort analysis",
},

"QA2C-005": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "ABSTRACT_VALIDATION_TERM",
    "approved_text":
    "Orthogonal regulator activity proxies and a curated ligand–receptor resource were used for focused follow-up evaluation.",
},

"QA2C-007": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP existing inconclusive / association-only boundary language.",
},

"QA2C-008": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP existing no causal inference / no wet-lab validation boundary language.",
},

"QA2C-012": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP: cell-level statistics can overstate reproducibility.",
},

"QA2C-013": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "Defining reproducible treatment-associated cell programs requires",
},

"QA2C-014": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "INTRO_DISCOVERY_SCOPE_DUPLICATE",
    "approved_text":
    "reproducible, sample-level evidence across the two discovery cohorts",
},

"QA2C-015": {
    "decision": "DUPLICATE_LOCK",
    "edit_group": "INTRO_DISCOVERY_SCOPE_DUPLICATE",
    "approved_text":
    "Same manuscript edit as QA2C-014; DO NOT edit twice.",
},

"QA2C-016": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "two discovery cohorts at the sample level to identify a shared treatment-associated fibroblast response",
},

"QA2C-017": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "",
    "approved_text":
    "Our objective was to define a fibroblast response reproducible across the two discovery cohorts",
},

"QA2C-018": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "",
    "approved_text":
    "and to distinguish the stronger discovery-cohort-reproducible effects represented by the CORE2 module from",
},

"QA2C-019": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP prediction / causality / validated-mechanism boundary language.",
},

"QA2C-020": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "same-direction positive treatment-associated effects in both discovery cohorts",
},

"QA2C-021": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "",
    "approved_text":
    "Across the two discovery cohorts, the CORE2 response was the strongest and most reproducible project-level treatment-associated result.",
},

"QA2C-022": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "same-direction in both discovery cohorts",
},

"QA2C-023": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "the most reproducible component of the fibroblast treatment-associated response across the two discovery cohorts.",
},

"QA2C-024": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP: master-regulator validation was overall inconclusive.",
},

"QA2C-025": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP expression-level / coverage-limited / not upgraded to validated boundary.",
},

"QA2C-026": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP: These remain association-supported candidates, not validated drivers.",
},

"QA2C-027": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "",
    "approved_text":
    "showed discordant association directions across the two discovery cohorts and were not conserved between them.",
},

"QA2C-028": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP negative boundary: discordant comparators are not supported as conserved mechanisms; signaling remains association-only.",
},

"QA2C-030": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "DISC_CORE2_OPENING_DUPLICATE",
    "approved_text":
    "Across the two discovery cohorts, the fibroblast CDKN1B–GSN CORE2 response was the most reproducible treatment-associated fibroblast result.",
},

"QA2C-031": {
    "decision": "DUPLICATE_LOCK",
    "edit_group": "DISC_CORE2_OPENING_DUPLICATE",
    "approved_text":
    "Same manuscript edit as QA2C-030; DO NOT edit twice.",
},

"QA2C-032": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "DISC_CORE2_VS_REGULATOR_SENTENCE",
    "approved_text":
    "Within the two discovery cohorts, the CORE2 response was more reproducible than any individual upstream ligand or regulator candidate evaluated in this project.",
},

"QA2C-033": {
    "decision": "MERGED_LOCK",
    "edit_group": "DISC_CORE2_VS_REGULATOR_SENTENCE",
    "approved_text":
    "Handled by the full-sentence rewrite locked under QA2C-032; do not retain 'regulator mechanism identified'.",
},

"QA2C-034": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP orthogonal validation inconclusive boundary language.",
},

"QA2C-035": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP not mechanistically validated / candidate hypotheses boundary language.",
},

"QA2C-036": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP highest focused association scores, not mechanistic confirmation.",
},

"QA2C-037": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP negative boundary: AREG, DLL1, EREG, and TGFB3 should not be reported as conserved axes.",
},

"QA2C-038": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "is essential for evaluating cross-cohort reproducibility",
},

"QA2C-039": {
    "decision": "CONTEXTUAL_REWRITE",
    "edit_group": "",
    "approved_text":
    "Within the two discovery cohorts, the fibroblast CORE2 response was the most stable treatment-associated result of this project.",
},

"QA2C-041": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP full context-dependent / not-uniformly-replicated / no-biological-cause boundary paragraph.",
},

"QA2C-043": {
    "decision": "KEEP",
    "edit_group": "",
    "approved_text":
    "KEEP full discovery-cohort reproducibility versus OMIX secondary directional-discordance paragraph.",
},

"QA2C-044": {
    "decision": "REWRITE_FRAGMENT",
    "edit_group": "",
    "approved_text":
    "5. Direction orientation: treated minus control, applied consistently across datasets",
},

"QA2C-045": {
    "decision": "REWRITE_HEADING",
    "edit_group": "",
    "approved_text":
    "8. Orthogonal regulator evidence overview (Step19N): signed perturbation",
},

"QA2C-046": {
    "decision": "REWRITE_HEADING",
    "edit_group": "",
    "approved_text":
    "11. Focused signaling association analysis (Step19P/19Q): treatment-adjusted,",
},

"QA2C-047": {
    "decision": "REWRITE_HEADING",
    "edit_group": "",
    "approved_text":
    "13. Analysis reproducibility and provenance: scripts in 07_scripts, logs in 08_logs, tables in",
},

}


print()
print("=" * 88)
print("STEP23C-QA2G: FREEZE REMAINING 39 HUMAN DECISIONS")
print("=" * 88)
print()


# ============================================================
# 1. Verify inputs
# ============================================================

for p in [QUEUE, OLD_LOCK, DOCX]:
    if not p.exists():
        raise SystemExit(f"Missing required input: {p}")

docx_sha_before = sha256(DOCX)

print("QA1B DOCX SHA256:")
print(" ", docx_sha_before)
print()

if docx_sha_before != EXPECTED_DOCX_SHA:
    raise SystemExit(
        "STOP: QA1B DOCX SHA256 changed since the frozen review."
    )

queue_rows = read_csv(QUEUE)
old_rows = read_csv(OLD_LOCK)

print("Existing QA2E locked items :", len(old_rows))
print("Remaining QA2F items       :", len(queue_rows))
print()

if len(old_rows) != 15:
    raise SystemExit(
        f"STOP: expected 15 QA2E locked rows, found {len(old_rows)}"
    )

if len(queue_rows) != 39:
    raise SystemExit(
        f"STOP: expected 39 remaining rows, found {len(queue_rows)}"
    )


# ============================================================
# 2. Verify exact ID set
# ============================================================

queue_by_id = {
    (r.get("item_id") or "").strip(): r
    for r in queue_rows
}

queue_ids = set(queue_by_id)
decision_ids = set(D)

missing = sorted(queue_ids - decision_ids)
extra = sorted(decision_ids - queue_ids)

if missing or extra:
    print("MISSING DECISIONS:", missing)
    print("EXTRA DECISIONS  :", extra)
    raise SystemExit("STOP: decision ID set does not equal frozen 39-item queue.")

print("✓ Exact 39/39 QA2F item IDs matched.")
print()


# ============================================================
# 3. Create lock rows
# ============================================================

lock_rows = []

for item_id in sorted(
    queue_ids,
    key=lambda x: int(x.split("-")[-1])
):
    src = queue_by_id[item_id]
    dec = D[item_id]

    lock_rows.append({
        "item_id": item_id,
        "section": (src.get("section") or "").strip(),
        "paragraph": (src.get("paragraph") or "").strip(),
        "severity": (src.get("severity") or "").strip(),
        "risk_family": (src.get("risk_family") or "").strip(),
        "review_bucket": (src.get("review_bucket") or "").strip(),
        "original_text": (src.get("text") or "").strip(),
        "qa2_reason": (src.get("reason") or "").strip(),
        "human_decision": dec["decision"],
        "edit_group": dec["edit_group"],
        "approved_text_or_instruction": dec["approved_text"],
        "manuscript_edit_authorized": "NO_NOT_YET_QA2B",
    })


write_csv(
    OUT_CSV,
    lock_rows,
    [
        "item_id",
        "section",
        "paragraph",
        "severity",
        "risk_family",
        "review_bucket",
        "original_text",
        "qa2_reason",
        "human_decision",
        "edit_group",
        "approved_text_or_instruction",
        "manuscript_edit_authorized",
    ]
)


# ============================================================
# 4. Human-readable dossier
# ============================================================

lines = [
    "STEP23C-QA2G REMAINING 39 HUMAN DECISION DOSSIER",
    "=" * 88,
    "",
    "STATUS:",
    "DECISIONS_LOCKED__MANUSCRIPT_NOT_EDITED",
    "",
    f"QA1B SHA256: {docx_sha_before}",
    "",
    "Previous QA2E decisions: 15",
    "New QA2G decisions: 39",
    "Total QA2C dispositions: 54",
    "",
    "IMPORTANT:",
    "This dossier authorizes NO manuscript edit yet.",
    "QA2B must apply edits contextually and must respect duplicate/merged edit groups.",
    "",
]

for r in lock_rows:
    lines.extend([
        "=" * 88,
        f"{r['item_id']}  |  {r['section']}  |  paragraph {r['paragraph']}",
        "=" * 88,
        "",
        f"Decision   : {r['human_decision']}",
        f"Edit group : {r['edit_group']}",
        "",
        "ORIGINAL:",
        r["original_text"],
        "",
        "LOCKED HUMAN DECISION / TEXT:",
        r["approved_text_or_instruction"],
        "",
    ])

OUT_TXT.write_text("\n".join(lines), encoding="utf-8")


# ============================================================
# 5. Verify DOCX unchanged
# ============================================================

docx_sha_after = sha256(DOCX)

if docx_sha_after != docx_sha_before:
    raise SystemExit("CRITICAL STOP: QA1B DOCX changed during QA2G.")


# ============================================================
# 6. Checksums
# ============================================================

with open(OUT_SHA, "w", encoding="utf-8") as f:
    for p in [OUT_CSV, OUT_TXT]:
        f.write(f"{sha256(p)}  {p.relative_to(ROOT)}\n")


# ============================================================
# 7. Summary
# ============================================================

counts = Counter(r["human_decision"] for r in lock_rows)

print("=" * 88)
print("STEP23C-QA2G COMPLETE")
print("=" * 88)
print()

print("QA1B DOCX MODIFIED          : NO")
print("STEP22A RECOMPUTED          : NO")
print("AUTOMATIC REWRITES          : NO")
print("GLOBAL FIND/REPLACE         : NOT AUTHORIZED")
print()

print("PREVIOUS QA2E LOCK          : 15")
print("NEW QA2G LOCK               : 39 / 39")
print("TOTAL QA2C DISPOSITIONS     : 54 / 54")
print()

print("QA2G DECISION COUNTS:")
for k, v in sorted(counts.items()):
    print(f"  {k:<24s}: {v}")

print()
print("OUTPUTS:")
print(" ", OUT_CSV.relative_to(ROOT))
print(" ", OUT_TXT.relative_to(ROOT))
print(" ", OUT_SHA.relative_to(ROOT))
print()

print("STATUS:")
print("STEP23C_QA2G_REMAINING_39_DECISIONS_LOCKED")
print()

print("NEXT:")
print("PREPARE QA2B CONTEXT-AWARE MANUSCRIPT EDIT ON A NEW COPY.")
print()

