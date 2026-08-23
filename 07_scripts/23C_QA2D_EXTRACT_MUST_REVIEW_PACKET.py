#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
from collections import Counter, defaultdict

PROJECT = Path.cwd()

META = (
    PROJECT
    / "00_metadata"
    / "Step23_manuscript_update"
)

RESULT = (
    PROJECT
    / "04_results"
    / "manuscript_update"
    / "Step23"
)

DOCX = (
    RESULT
    / "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"
)

QUEUE = (
    META
    / "STEP23C_QA2C_54_ITEM_REVIEW_QUEUE.csv"
)

OUT_TXT = (
    RESULT
    / "STEP23C_QA2D_MUST_REVIEW_PACKET.txt"
)

OUT_CSV = (
    META
    / "STEP23C_QA2D_MUST_REVIEW_PACKET.csv"
)

STATUS = (
    META
    / "STEP23C_QA2D_STATUS.csv"
)


def fail(msg):
    print()
    print("ERROR:")
    print(msg)
    print()
    raise SystemExit(1)


def clean(x):
    return str(x or "").replace("\r", "").strip()


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)

    return h.hexdigest()


def read_csv(path):
    with open(
        path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with open(
        path,
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        w = csv.DictWriter(
            f,
            fieldnames=fields,
            extrasaction="ignore"
        )

        w.writeheader()
        w.writerows(rows)


if not DOCX.exists():
    fail(f"QA1B DOCX missing:\n{DOCX}")

if not QUEUE.exists():
    fail(f"QA2C review queue missing:\n{QUEUE}")


docx_sha_before = sha256_file(DOCX)

rows = read_csv(QUEUE)

if len(rows) != 54:
    fail(
        f"Expected 54 frozen QA2C rows; found {len(rows)}"
    )


# ============================================================
# Select the first human-review packet
# ============================================================

must = [
    x for x in rows
    if clean(x.get("review_bucket"))
    == "MUST_REVIEW_AND_REVISE"
]

packaging = [
    x for x in rows
    if clean(x.get("risk_family"))
    == "PACKAGING_COMPLETENESS"
]

# union without duplicates
selected = []
seen = set()

for x in must + packaging:

    key = (
        clean(x.get("item_id")),
        clean(x.get("text"))
    )

    if key in seen:
        continue

    seen.add(key)
    selected.append(x)


# preserve existing review order
selected.sort(
    key=lambda x:
        int(clean(x.get("review_order")) or 9999)
)


# ============================================================
# Add deterministic editorial ACTION TYPE only.
#
# IMPORTANT:
# These are NOT rewritten sentences.
# They tell us what kind of human revision is needed.
# ============================================================

def action_type(row):

    section = clean(
        row.get("section")
    )

    family = clean(
        row.get("risk_family")
    )

    text = clean(
        row.get("text")
    ).lower()


    if family == "PACKAGING_COMPLETENESS":

        if section == "Data Availability":
            return (
                "ADD_OMIX_DATASET_SOURCE_AND_ACCESSION_OR_PUBLIC_IDENTIFIER"
            )

        if section == "Supplementary":
            return (
                "ADD_STEP22A_EXTERNAL_VALIDATION_SUPPLEMENTARY_ENTRIES"
            )

        if section == "References":
            return (
                "ADD_REAL_SOURCE_CITATION_FOR_OMIX005710_AND_RELEVANT_SOURCE_STUDY"
            )

        return "COMPLETE_SUBMISSION_PACKAGE_INFORMATION"


    if section == "Abstract":

        return (
            "RECALIBRATE_CROSS_COHORT_CLAIM_AND_PRESERVE_PRIMARY_VS_SECONDARY_BOUNDARY"
        )


    if section == "Results":

        return (
            "DISTINGUISH_DISCOVERY_REPLICATION_FROM_OMIX_SECONDARY_NONREPLICATION"
        )


    if section == "Discussion":

        return (
            "LIMIT_GENERALIZATION_AND_STATE_CONTEXT_DEPENDENT_EXTERNAL_DISCORDANCE"
        )


    if section == "Methods":

        if (
            "validation"
            in text
            or
            family == "VALIDATION_STRENGTH"
        ):
            return (
                "CLARIFY_PRIMARY_UNAVAILABLE_VS_SECONDARY_SENSITIVITY_DESIGN"
            )

        return (
            "CLARIFY_SCOPE_WITHOUT_CHANGING_ANALYTIC_METHOD"
        )


    return "HUMAN_SCIENTIFIC_REWRITE_REQUIRED"


for row in selected:
    row["action_type"] = action_type(row)


# ============================================================
# Summary
# ============================================================

by_section = Counter(
    clean(x.get("section"))
    for x in selected
)

by_family = Counter(
    clean(x.get("risk_family"))
    for x in selected
)

print()
print("=" * 82)
print("STEP23C-QA2D: MUST-REVIEW PACKET")
print("=" * 82)
print()

print(
    "Frozen QA2C items total        :",
    len(rows)
)

print(
    "MUST_REVIEW_AND_REVISE         :",
    len(must)
)

print(
    "Packaging-completeness flags   :",
    len(packaging)
)

print(
    "Unique items in review packet  :",
    len(selected)
)

print()

print("BY SECTION:")

for section, n in by_section.most_common():
    print(
        f"  {section:<24s}: {n}"
    )

print()

print("BY RISK FAMILY:")

for family, n in by_family.most_common():
    print(
        f"  {family:<32s}: {n}"
    )

print()


# ============================================================
# Save CSV
# ============================================================

write_csv(
    OUT_CSV,
    selected,
    [
        "review_order",
        "item_id",
        "severity",
        "section",
        "paragraph",
        "risk_family",
        "review_bucket",
        "action_type",
        "trigger",
        "rule",
        "reason",
        "qa2_suggestion",
        "text"
    ]
)


# ============================================================
# Human-readable review packet
# ============================================================

lines = [
    "STEP23C-QA2D MUST-REVIEW SCIENTIFIC WORDING PACKET",
    "=" * 82,
    "",
    "STATUS:",
    "READ_ONLY_NO_MANUSCRIPT_EDIT",
    "",
    f"QA1B manuscript SHA256: {docx_sha_before}",
    "",
    f"Frozen QA2C total items: {len(rows)}",
    f"MUST_REVIEW_AND_REVISE: {len(must)}",
    f"Packaging flags: {len(packaging)}",
    f"Unique packet items: {len(selected)}",
    "",
    "IMPORTANT:",
    "No rewritten language has been generated automatically.",
    "The ACTION TYPE below only indicates what kind of human",
    "scientific/editorial decision is required.",
    "",
]


current_section = None

for i, row in enumerate(selected, start=1):

    section = clean(
        row.get("section")
    )

    if section != current_section:

        current_section = section

        lines.extend([
            "",
            "=" * 82,
            section.upper(),
            "=" * 82,
            ""
        ])


    lines.extend([
        f"{i:02d}. QA2C item {clean(row.get('item_id'))}",
        f"Severity     : {clean(row.get('severity'))}",
        f"Risk family  : {clean(row.get('risk_family'))}",
        f"Paragraph    : {clean(row.get('paragraph'))}",
        f"Trigger      : {clean(row.get('trigger'))}",
        f"Action type  : {clean(row.get('action_type'))}",
        "",
        "ORIGINAL TEXT:",
        clean(row.get("text")),
        ""
    ])


    reason = clean(
        row.get("reason")
    )

    if reason:

        lines.extend([
            "QA2 REASON:",
            reason,
            ""
        ])


    suggestion = clean(
        row.get("qa2_suggestion")
    )

    if suggestion:

        lines.extend([
            "EXISTING QA2 SUGGESTION:",
            suggestion,
            ""
        ])


    lines.extend([
        "HUMAN REVISION:",
        "[PENDING — DO NOT AUTO-REWRITE]",
        "",
        "-" * 82,
        ""
    ])


lines.extend([
    "",
    "=" * 82,
    "NEXT DECISION",
    "=" * 82,
    "",
    "This packet should now be reviewed item by item.",
    "",
    "Recommended order:",
    "",
    "1. Abstract direct-conflict items.",
    "2. Results / Discussion direct generalization conflicts.",
    "3. Methods primary-vs-secondary validation terminology.",
    "4. Data Availability.",
    "5. Supplementary external-validation index.",
    "6. References / real-source citation completion.",
    "",
    "After these are approved, process discovery-scope qualification",
    "items separately; do not use global find/replace.",
    "",
    "STATUS:",
    "STEP23C_QA2D_PACKET_READY_FOR_HUMAN_REWRITE",
    ""
])


OUT_TXT.write_text(
    "\n".join(lines),
    encoding="utf-8"
)


# ============================================================
# Confirm DOCX unchanged
# ============================================================

docx_sha_after = sha256_file(DOCX)

if docx_sha_after != docx_sha_before:
    fail(
        "CRITICAL: QA1B DOCX changed during read-only QA2D"
    )


# ============================================================
# Status
# ============================================================

status_rows = [
    {
        "item":
            "QA1B_DOCX_modified",

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
            "automatic_rewrites",

        "status":
            "NO"
    },

    {
        "item":
            "MUST_REVIEW_count",

        "status":
            str(len(must))
    },

    {
        "item":
            "packaging_count",

        "status":
            str(len(packaging))
    },

    {
        "item":
            "review_packet_unique_items",

        "status":
            str(len(selected))
    },

    {
        "item":
            "status",

        "status":
            "STEP23C_QA2D_PACKET_READY_FOR_HUMAN_REWRITE"
    }
]


write_csv(
    STATUS,
    status_rows,
    [
        "item",
        "status"
    ]
)


print("=" * 82)
print("STEP23C-QA2D COMPLETE")
print("=" * 82)
print()

print(
    "QA1B DOCX MODIFIED            : NO"
)

print(
    "STEP22A RECOMPUTED            : NO"
)

print(
    "AUTOMATIC REWRITES            : NO"
)

print()

print(
    "MUST REVIEW ITEMS             :",
    len(must)
)

print(
    "PACKAGING FLAGS               :",
    len(packaging)
)

print(
    "UNIQUE REVIEW PACKET ITEMS    :",
    len(selected)
)

print()

print("OUTPUTS:")

print(
    " ",
    OUT_TXT.relative_to(PROJECT)
)

print(
    " ",
    OUT_CSV.relative_to(PROJECT)
)

print()

print("STATUS:")

print(
    "STEP23C_QA2D_PACKET_READY_FOR_HUMAN_REWRITE"
)

print()

print(
    "NEXT: REVIEW THESE ITEMS BEFORE QA2B EDITING."
)

print()
