#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP23A
#
# MANUSCRIPT INTEGRATION PREFLIGHT
#
# PURPOSE
#   Prepare a NEW manuscript-update branch after Step22A closure.
#
# THIS STEP DOES NOT:
#   - modify Step19-22 outputs
#   - modify Step21E manuscript
#   - edit DOCX/PDF
#   - regenerate figures
#   - recompute any analysis
#   - change scientific conclusions
#
# It ONLY:
#   1. verifies Step22A is formally archived
#   2. locates manuscript candidates
#   3. locates Step22A reporting materials
#   4. creates a clean Step23 workspace
#   5. writes inventories / integration plan
# ============================================================

from pathlib import Path
import csv
import hashlib
import datetime
import re
import sys


PROJECT = Path.cwd()

META22 = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

RESULT22 = (
    PROJECT
    / "04_results"
    / "external_validation"
    / "Step22A"
)

RESULT_ROOT = PROJECT / "04_results"
TABLE_ROOT = PROJECT / "06_tables"
FIG_ROOT = PROJECT / "05_figures"
SCRIPT_ROOT = PROJECT / "07_scripts"
LOG_ROOT = PROJECT / "08_logs"

# ------------------------------------------------------------
# New Step23 workspace
# ------------------------------------------------------------

STEP23_ROOT = (
    RESULT_ROOT
    / "manuscript_update"
    / "Step23"
)

STEP23_META = (
    PROJECT
    / "00_metadata"
    / "Step23_manuscript_update"
)

STEP23_TABLE = (
    TABLE_ROOT
    / "manuscript_update"
    / "Step23"
)

STEP23_FIG = (
    FIG_ROOT
    / "manuscript_update"
    / "Step23"
)

for d in [
    STEP23_ROOT,
    STEP23_META,
    STEP23_TABLE,
    STEP23_FIG
]:
    d.mkdir(
        parents=True,
        exist_ok=True
    )

# ------------------------------------------------------------
# Required Step22A archival artifacts
# ------------------------------------------------------------

STEP22_CLOSEOUT = (
    META22
    / "STEP22A_FINAL_CLOSEOUT.txt"
)

STEP22_README = (
    META22
    / "STEP22A_FINAL_PROJECT_README.md"
)

STEP22_MANIFEST = (
    META22
    / "STEP22A_FINAL_FILE_MANIFEST.csv"
)

STEP22_CHECKSUMS = (
    META22
    / "STEP22A_FINAL_CHECKSUMS.sha256"
)

STEP22_7D_LOCK = (
    META22
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_LOCK.csv"
)

STEP22_7D_DOSSIER = (
    META22
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_DOSSIER.txt"
)

STEP22_REPORTING = (
    RESULT22
    / "STEP22A_07D_REPORTING_LANGUAGE.txt"
)

STEP22_FINAL_SUMMARY = (
    RESULT22
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_SUMMARY.txt"
)

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

MANUSCRIPT_INVENTORY = (
    STEP23_META
    / "STEP23A_MANUSCRIPT_CANDIDATE_INVENTORY.csv"
)

STEP22_INPUT_INVENTORY = (
    STEP23_META
    / "STEP23A_STEP22A_INTEGRATION_INPUTS.csv"
)

FIGURE_INVENTORY = (
    STEP23_META
    / "STEP23A_EXISTING_FIGURE_TABLE_INVENTORY.csv"
)

PLAN_OUT = (
    STEP23_ROOT
    / "STEP23A_MANUSCRIPT_INTEGRATION_PLAN.txt"
)

STATUS_OUT = (
    STEP23_META
    / "STEP23A_PREFLIGHT_STATUS.csv"
)


print()
print("=" * 80)
print("STEP23A: MANUSCRIPT INTEGRATION PREFLIGHT")
print("=" * 80)
print()


# ============================================================
# Helpers
# ============================================================

def fail(msg):
    print()
    print("ERROR:")
    print(msg)
    print()
    raise SystemExit(1)


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            block = f.read(
                1024 * 1024
            )

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def relative(path):
    try:
        return str(
            path.relative_to(PROJECT)
        )
    except Exception:
        return str(path)


def write_csv(
    path,
    rows,
    fields
):

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def manuscript_score(path):

    text = str(path).lower()

    score = 0

    # Strong historical manuscript markers
    if "step21e" in text:
        score += 100

    if "final_manuscript" in text:
        score += 80

    if "manuscript" in text:
        score += 50

    if "submission" in text:
        score += 40

    if "revised" in text:
        score += 20

    if "final" in path.name.lower():
        score += 15

    # Preferred editable source
    if path.suffix.lower() == ".docx":
        score += 30

    elif path.suffix.lower() in {
        ".md",
        ".txt"
    }:
        score += 10

    # PDF useful for reference, not editing source
    elif path.suffix.lower() == ".pdf":
        score += 5

    # Avoid external validation outputs being mistaken
    if "step22a" in text:
        score -= 200

    if "step23" in text:
        score -= 200

    return score


# ============================================================
# 1. Confirm project root
# ============================================================

print("===== 1. PROJECT ROOT =====")
print()

required_root_dirs = [
    "00_metadata",
    "04_results",
    "06_tables",
    "07_scripts",
    "08_logs"
]

missing_dirs = [
    x
    for x in required_root_dirs
    if not (
        PROJECT
        / x
    ).exists()
]

print(
    "Project root:",
    PROJECT
)

print(
    "Missing expected root directories:",
    len(
        missing_dirs
    )
)

if missing_dirs:

    print(
        ", ".join(
            missing_dirs
        )
    )

    fail(
        "PROJECT_ROOT_NOT_CONFIRMED"
    )

print(
    "✓ Project root confirmed."
)

print()


# ============================================================
# 2. Verify Step22A is CLOSED
# ============================================================

print("===== 2. VERIFY STEP22A ARCHIVE =====")
print()

required_22 = [
    STEP22_CLOSEOUT,
    STEP22_README,
    STEP22_MANIFEST,
    STEP22_CHECKSUMS,
    STEP22_7D_LOCK,
    STEP22_7D_DOSSIER,
    STEP22_REPORTING,
    STEP22_FINAL_SUMMARY
]

step22_rows = []

for path in required_22:

    exists = path.exists()

    print(
        f"{path.name:<65s}",
        "FOUND"
        if exists
        else "MISSING"
    )

    step22_rows.append({
        "file":
            path.name,

        "path":
            relative(path),

        "exists":
            "YES"
            if exists
            else "NO",

        "size_bytes":
            (
                path.stat().st_size
                if exists
                else ""
            ),

        "sha256":
            (
                sha256_file(path)
                if exists
                else ""
            )
    })


if not all(
    x["exists"] == "YES"
    for x in step22_rows
):

    fail(
        "STEP22A_ARCHIVE_INCOMPLETE"
    )


closeout_text = STEP22_CLOSEOUT.read_text(
    encoding="utf-8",
    errors="replace"
)


required_closeout_phrase = (
    "STEP22A_EXTERNAL_VALIDATION_ARCHIVED_AND_CLOSED"
)


if required_closeout_phrase not in closeout_text:

    fail(
        "STEP22A_CLOSEOUT_STATUS_NOT_CONFIRMED"
    )


print()
print(
    "✓ Step22A is formally archived and CLOSED."
)
print()


write_csv(
    STEP22_INPUT_INVENTORY,
    step22_rows,
    [
        "file",
        "path",
        "exists",
        "size_bytes",
        "sha256"
    ]
)


# ============================================================
# 3. Locate manuscript candidates
#
# Search project only.
# Do not edit or copy anything.
# ============================================================

print("===== 3. LOCATE EXISTING MANUSCRIPT =====")
print()

allowed_suffixes = {
    ".docx",
    ".md",
    ".txt",
    ".pdf"
}

candidate_paths = []


for path in PROJECT.rglob("*"):

    if not path.is_file():
        continue

    if path.suffix.lower() not in allowed_suffixes:
        continue

    rel_lower = relative(
        path
    ).lower()

    # Ignore Step22A / Step23 / logs / raw downloads.
    if (
        "step22a" in rel_lower
        or
        "step23" in rel_lower
        or
        "/08_logs/" in
            (
                "/"
                + rel_lower
            )
        or
        "/01_raw/" in
            (
                "/"
                + rel_lower
            )
    ):
        continue

    # Keep files likely related to manuscript/submission.
    manuscript_terms = [
        "manuscript",
        "submission",
        "step21e",
        "final",
        "revised"
    ]

    if any(
        term in rel_lower
        for term in manuscript_terms
    ):

        candidate_paths.append(
            path
        )


candidate_rows = []

for path in candidate_paths:

    score = manuscript_score(
        path
    )

    candidate_rows.append({
        "rank_score":
            score,

        "path":
            relative(path),

        "filename":
            path.name,

        "suffix":
            path.suffix.lower(),

        "size_bytes":
            path.stat().st_size,

        "modified_time":
            datetime.datetime.fromtimestamp(
                path.stat().st_mtime
            ).isoformat(
                timespec="seconds"
            ),

        "sha256":
            sha256_file(
                path
            )
    })


candidate_rows.sort(
    key=lambda x: (
        -x[
            "rank_score"
        ],
        x[
            "path"
        ]
    )
)


for i, row in enumerate(
    candidate_rows,
    start=1
):

    row["rank"] = i


write_csv(
    MANUSCRIPT_INVENTORY,
    candidate_rows,
    [
        "rank",
        "rank_score",
        "path",
        "filename",
        "suffix",
        "size_bytes",
        "modified_time",
        "sha256"
    ]
)


print(
    "Manuscript candidates found:",
    len(
        candidate_rows
    )
)

print()


for row in candidate_rows[:15]:

    print(
        f"[{row['rank']:02d}] "
        f"score={row['rank_score']:3d}  "
        f"{row['path']}"
    )


print()


editable_candidates = [
    x
    for x in candidate_rows
    if x[
        "suffix"
    ]
    in {
        ".docx",
        ".md",
        ".txt"
    }
]


if not editable_candidates:

    fail(
        "NO_EDITABLE_MANUSCRIPT_CANDIDATE_FOUND"
    )


top_candidate = editable_candidates[0]


print(
    "Top editable manuscript candidate:"
)

print(
    "  ",
    top_candidate[
        "path"
    ]
)

print()


# ============================================================
# 4. Inventory existing figures/tables relevant to manuscript
# ============================================================

print("===== 4. INVENTORY EXISTING FIGURES / TABLES =====")
print()

figure_table_rows = []


search_roots = [
    FIG_ROOT,
    TABLE_ROOT,
    RESULT_ROOT
]


for root in search_roots:

    if not root.exists():
        continue

    for path in root.rglob("*"):

        if not path.is_file():
            continue

        rel = relative(
            path
        )

        rel_lower = rel.lower()

        # Skip Step23 outputs being created now.
        if "step23" in rel_lower:
            continue

        interesting = (
            "step19" in rel_lower
            or
            "step20" in rel_lower
            or
            "step21" in rel_lower
            or
            "step22a" in rel_lower
            or
            "core2" in rel_lower
            or
            "core4" in rel_lower
            or
            "caf" in rel_lower
            or
            "external_validation" in rel_lower
        )

        if not interesting:
            continue

        figure_table_rows.append({
            "path":
                rel,

            "filename":
                path.name,

            "suffix":
                path.suffix.lower(),

            "size_bytes":
                path.stat().st_size,

            "modified_time":
                datetime.datetime.fromtimestamp(
                    path.stat().st_mtime
                ).isoformat(
                    timespec="seconds"
                )
        })


figure_table_rows.sort(
    key=lambda x:
        x["path"]
)


write_csv(
    FIGURE_INVENTORY,
    figure_table_rows,
    [
        "path",
        "filename",
        "suffix",
        "size_bytes",
        "modified_time"
    ]
)


print(
    "Relevant existing result/figure/table files:",
    len(
        figure_table_rows
    )
)

print()


# ============================================================
# 5. Freeze Step23 scientific hierarchy
#
# Important:
# This is a manuscript integration hierarchy, not new analysis.
# ============================================================

print("===== 5. STEP23 INTEGRATION HIERARCHY =====")
print()


integration_hierarchy = [
    (
        "DISCOVERY",
        "GSE197677 / GSE221561",
        "Frozen discovery cohorts and previously established CORE framework"
    ),

    (
        "PRIMARY_EXTERNAL_OMIX",
        "OMIX005710",
        "UNAVAILABLE because individual EASC identity remains unresolved"
    ),

    (
        "SECONDARY_EXTERNAL_OMIX",
        "OMIX005710",
        "22/22 pre-specified histology scenarios show robustly opposite CORE direction"
    ),

    (
        "TECHNICAL_ROBUSTNESS",
        "OMIX005710",
        "Independent recalculation plus 105/105 paired LOPO support non-artifactual directional discordance"
    ),

    (
        "BIOLOGICAL_CAUSE",
        "OMIX005710",
        "NOT ESTABLISHED; no post-hoc mechanism analysis performed"
    )
]


for level, dataset, conclusion in integration_hierarchy:

    print(
        f"{level:<26s} | "
        f"{dataset:<24s} | "
        f"{conclusion}"
    )


print()


# ============================================================
# 6. Write integration plan
# ============================================================

plan = f"""STEP23A MANUSCRIPT INTEGRATION PREFLIGHT
======================================================================

STATUS:
READY_FOR_MANUSCRIPT_INTEGRATION

STEP22A:
ARCHIVED_AND_CLOSED

STEP22A MAY BE RECOMPUTED:
NO

STEP21E ORIGINAL MANUSCRIPT MAY BE OVERWRITTEN:
NO

TOP EDITABLE MANUSCRIPT CANDIDATE:
{top_candidate["path"]}

TOP MANUSCRIPT SHA256:
{top_candidate["sha256"]}

MANUSCRIPT CANDIDATES FOUND:
{len(candidate_rows)}

RELEVANT EXISTING RESULT / FIGURE / TABLE FILES:
{len(figure_table_rows)}

----------------------------------------------------------------------
SCIENTIFIC EVIDENCE HIERARCHY
----------------------------------------------------------------------

1. DISCOVERY / ORIGINAL FROZEN EVIDENCE

   GSE197677 and GSE221561 remain the frozen discovery evidence
   supporting the previously defined CORE framework.

2. OMIX005710 PRIMARY ESCC-ONLY EXTERNAL VALIDATION

   UNAVAILABLE.

   Reason:
   the single EASC patient's individual identity cannot be resolved
   from the available public patient-level metadata.

   This must not be rewritten as a negative primary validation.

3. OMIX005710 PRE-SPECIFIED SECONDARY HISTOLOGY SENSITIVITY

   All 22 hypothetical EASC assignments were evaluated.

   CORE2_STRONG treatment-related direction:
   opposite to frozen discovery direction in 22/22 scenarios.

   CORE4:
   opposite direction in 22/22 scenarios.

   No histology scenario was selected or promoted to primary.

4. TECHNICAL ROBUSTNESS

   Independent score reconstruction reproduced Step22A within
   floating-point precision.

   Paired LOPO:
   105/105 runs retained the negative CORE2 direction.

   Therefore the secondary directional discordance is not explained
   by a simple sign, pairing, reference-membership, P15, or
   single-patient artifact.

5. BIOLOGICAL EXPLANATION

   NOT ESTABLISHED.

   No post-hoc CAF subtype or mechanistic analysis was performed.

   Any future mechanism work must be explicitly labeled exploratory
   and hypothesis-generating.

----------------------------------------------------------------------
STEP23 MANUSCRIPT PLAN
----------------------------------------------------------------------

STEP23A
  Preflight / manuscript source identification
  STATUS: COMPLETE after this script.

STEP23B
  Freeze manuscript-level three-cohort evidence statements.
  Produce a section-by-section change map:
    - Abstract
    - Methods
    - Results
    - Discussion
    - Limitations
    - Supplementary material
  No DOCX editing yet.

STEP23C
  Create NEW manuscript working copy.
  Never overwrite Step21E source.
  Integrate:
    - OMIX005710 Methods
    - secondary sensitivity Results
    - negative/non-replication interpretation
    - histology limitation
    - LOPO robustness
  Keep biological mechanism claims restrained.

STEP23D
  Build supplementary external-validation table(s) and provenance
  summary from frozen Step22A files only.

STEP23E
  Assemble updated manuscript package.

STEP23F
  Final consistency audit:
    - cohort counts
    - patient counts
    - sample counts
    - CORE definitions
    - direction wording
    - primary vs secondary hierarchy
    - figure/table references
    - no unsupported mechanistic claims

----------------------------------------------------------------------
STRICT REPORTING LANGUAGE
----------------------------------------------------------------------

DO SAY:

  "Primary ESCC-only validation was unavailable because the single
   EASC patient could not be resolved at the individual-patient level
   from the available public metadata."

  "In the pre-specified secondary histology-sensitivity analysis,
   OMIX005710 showed a treatment-associated CORE direction opposite
   to that observed in the discovery cohorts across all 22
   hypothetical EASC assignments."

  "The directional discordance remained robust in independent score
   reconstruction and paired leave-one-patient-out analyses."

DO NOT SAY:

  "OMIX005710 definitively validated the opposite mechanism."

  "The third cohort proves that CORE2 decreases after treatment."

  "The ESCC-only external validation failed."

  "One CAF subtype explains the discordance."

  "CORE4 rescued the validation."

----------------------------------------------------------------------
NEXT STEP
----------------------------------------------------------------------

STEP23B:
MANUSCRIPT EVIDENCE-MAP FREEZE

No manuscript text will be modified until Step23B is complete.
"""


PLAN_OUT.write_text(
    plan,
    encoding="utf-8"
)


# ============================================================
# 7. Final preflight status
# ============================================================

status_rows = [
    {
        "item":
            "project_root_confirmed",
        "status":
            "PASS"
    },

    {
        "item":
            "Step22A_archived_and_closed",
        "status":
            "PASS"
    },

    {
        "item":
            "Step22A_required_archival_files",
        "status":
            "8_OF_8_FOUND"
    },

    {
        "item":
            "editable_manuscript_candidates",
        "status":
            str(
                len(
                    editable_candidates
                )
            )
    },

    {
        "item":
            "top_editable_manuscript_candidate",
        "status":
            top_candidate[
                "path"
            ]
    },

    {
        "item":
            "original_manuscript_modified",
        "status":
            "NO"
    },

    {
        "item":
            "Step22A_recomputed",
        "status":
            "NO"
    },

    {
        "item":
            "Step23_workspace_created",
        "status":
            "YES"
    },

    {
        "item":
            "next_step",
        "status":
            "STEP23B_MANUSCRIPT_EVIDENCE_MAP_FREEZE"
    }
]


write_csv(
    STATUS_OUT,
    status_rows,
    [
        "item",
        "status"
    ]
)


# ============================================================
# FINAL
# ============================================================

print("=" * 80)
print("STEP23A COMPLETE")
print("=" * 80)
print()

print(
    "STEP22A STATUS                : ARCHIVED / CLOSED"
)

print(
    "STEP22A RECOMPUTED            : NO"
)

print(
    "ORIGINAL MANUSCRIPT MODIFIED  : NO"
)

print()

print(
    "MANUSCRIPT CANDIDATES FOUND   :",
    len(
        candidate_rows
    )
)

print(
    "EDITABLE CANDIDATES           :",
    len(
        editable_candidates
    )
)

print()

print(
    "TOP EDITABLE MANUSCRIPT:"
)

print(
    " ",
    top_candidate[
        "path"
    ]
)

print()

print(
    "STEP23 WORKSPACE:"
)

print(
    " ",
    relative(
        STEP23_ROOT
    )
)

print()

print(
    "STEP23A PLAN:"
)

print(
    " ",
    relative(
        PLAN_OUT
    )
)

print()

print(
    "SCIENTIFIC HIERARCHY:"
)

print(
    "  DISCOVERY COHORTS            : FROZEN"
)

print(
    "  OMIX PRIMARY ESCC-ONLY       : UNAVAILABLE"
)

print(
    "  OMIX SECONDARY SENSITIVITY   : ROBUSTLY OPPOSITE"
)

print(
    "  TECHNICAL ARTIFACT           : NOT DETECTED"
)

print(
    "  BIOLOGICAL CAUSE             : NOT ESTABLISHED"
)

print()

print(
    "NEXT:"
)

print(
    "STEP23B_MANUSCRIPT_EVIDENCE_MAP_FREEZE"
)

print()

print(
    "DO NOT EDIT THE MANUSCRIPT YET."
)

print()
