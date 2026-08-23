#!/usr/bin/env python3

# ============================================================
# STEP23C-QA2F
#
# GROUP THE 39 REMAINING QA2 ITEMS FOR HUMAN REVIEW
#
# READ-ONLY:
#   - NO DOCX edit
#   - NO Step22A recomputation
#   - NO automatic scientific rewrite
#   - NO global find/replace
#
# Outputs:
#   1. discovery-scope packet
#   2. context-review packet
#   3. likely-keep-boundary packet
#   4. section/bucket summary
#
# ============================================================

from pathlib import Path
from collections import Counter
import csv
import hashlib
import re


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

REMAINING = (
    META
    / "STEP23C_QA2E_REMAINING_REVIEW_QUEUE.csv"
)

LOCK = (
    META
    / "STEP23C_QA2E_DECISION_LOCK.csv"
)


SCOPE_TXT = (
    RESULT
    / "STEP23C_QA2F_DISCOVERY_SCOPE_PACKET.txt"
)

CONTEXT_TXT = (
    RESULT
    / "STEP23C_QA2F_CONTEXT_REVIEW_PACKET.txt"
)

KEEP_TXT = (
    RESULT
    / "STEP23C_QA2F_KEEP_BOUNDARY_PACKET.txt"
)

SUMMARY_CSV = (
    META
    / "STEP23C_QA2F_REMAINING_SUMMARY.csv"
)

STATUS_CSV = (
    META
    / "STEP23C_QA2F_STATUS.csv"
)


def fail(msg):
    print()
    print("ERROR:")
    print(msg)
    print()
    raise SystemExit(1)


def clean(x):
    return str(x or "").replace("\r", "").strip()


def norm(x):
    return re.sub(r"\s+", " ", clean(x))


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


def relative(path):
    return str(path.relative_to(PROJECT))


# ============================================================
# 1. Verify inputs
# ============================================================

print()
print("=" * 84)
print("STEP23C-QA2F: GROUP REMAINING 39 ITEMS")
print("=" * 84)
print()


for path in [DOCX, REMAINING, LOCK]:

    if not path.exists():

        fail(
            "Missing required input:\n"
            + str(path)
        )


docx_sha_before = sha256_file(DOCX)

rows = read_csv(REMAINING)
lock_rows = read_csv(LOCK)


if len(rows) != 39:

    fail(
        f"Expected 39 remaining QA2 items; found {len(rows)}."
    )


if len(lock_rows) != 15:

    fail(
        f"Expected 15 locked QA2E decisions; found {len(lock_rows)}."
    )


print("QA1B SHA256:")
print(" ", docx_sha_before)
print()

print("Locked QA2E items :", len(lock_rows))
print("Remaining items   :", len(rows))
print()


# ============================================================
# 2. Verify bucket counts
# ============================================================

bucket_counts = Counter(
    clean(x.get("review_bucket"))
    for x in rows
)


EXPECTED = {
    "QUALIFY_TO_DISCOVERY_SCOPE": 22,
    "CONTEXT_REVIEW": 10,
    "LIKELY_KEEP_BOUNDARY_LANGUAGE": 5,
    "QUALIFY_OR_KEEP_DISCOVERY_SCOPE": 2,
}


print("===== FROZEN BUCKET COUNTS =====")
print()


for bucket, expected in EXPECTED.items():

    observed = bucket_counts.get(
        bucket,
        0
    )

    print(
        f"{bucket:<42s}: "
        f"{observed} "
        f"(expected {expected})"
    )

    if observed != expected:

        fail(
            f"Frozen bucket count changed for {bucket}."
        )


print()
print("✓ Remaining 39-item count lock confirmed.")
print()


# ============================================================
# 3. Add human-review priority labels
#
# No editorial decision is made here.
# ============================================================

def priority(row):

    section = clean(
        row.get("section")
    )

    bucket = clean(
        row.get("review_bucket")
    )

    text = norm(
        row.get("text")
    ).lower()


    if section in {
        "Title",
        "Running Title"
    }:
        return "P0_TITLE_RUNNING_TITLE"


    if section == "Abstract":
        return "P1_ABSTRACT"


    if bucket in {
        "QUALIFY_TO_DISCOVERY_SCOPE",
        "QUALIFY_OR_KEEP_DISCOVERY_SCOPE"
    }:

        if section in {
            "Introduction",
            "Results",
            "Discussion"
        }:
            return "P2_DISCOVERY_SCOPE_MAIN_TEXT"

        return "P3_DISCOVERY_SCOPE_OTHER"


    if bucket == "CONTEXT_REVIEW":

        return "P4_CONTEXT_REVIEW"


    if bucket == "LIKELY_KEEP_BOUNDARY_LANGUAGE":

        return "P5_KEEP_BOUNDARY_CHECK"


    return "P9_OTHER"


def suggested_review_question(row):

    bucket = clean(
        row.get("review_bucket")
    )

    section = clean(
        row.get("section")
    )

    family = clean(
        row.get("risk_family")
    )


    if section in {
        "Title",
        "Running Title"
    }:

        return (
            "Does this title wording over-generalize beyond "
            "GSE197677/GSE221561 after OMIX005710 directional "
            "non-replication?"
        )


    if bucket == "QUALIFY_TO_DISCOVERY_SCOPE":

        return (
            "Is the scientific statement still valid if explicitly "
            "limited to the two discovery cohorts?"
        )


    if bucket == "QUALIFY_OR_KEEP_DISCOVERY_SCOPE":

        return (
            "Does the existing sentence already make the two-discovery-"
            "cohort scope sufficiently explicit, or should one short "
            "scope qualifier be added?"
        )


    if bucket == "LIKELY_KEEP_BOUNDARY_LANGUAGE":

        return (
            "Does the existing not-validated / inconclusive / "
            "association-only / hypothesis-generating boundary already "
            "protect against overclaiming?"
        )


    if bucket == "CONTEXT_REVIEW":

        if family == "MECHANISM_CAUSALITY":

            return (
                "Does this wording imply mechanism, causality, driver, "
                "or validated regulator beyond the available evidence?"
            )

        if family == "VALIDATION_STRENGTH":

            return (
                "Does this wording correctly distinguish discovery "
                "replication, primary validation unavailable, and "
                "secondary sensitivity analysis?"
            )

        return (
            "Does the statement need scientific scope calibration after "
            "the OMIX005710 result?"
        )


    return "Human scientific review required."


for row in rows:

    row["qa2f_priority"] = priority(row)

    row[
        "qa2f_review_question"
    ] = suggested_review_question(row)


# ============================================================
# 4. Split into three packets
# ============================================================

scope_rows = [
    x for x in rows
    if clean(x.get("review_bucket"))
    in {
        "QUALIFY_TO_DISCOVERY_SCOPE",
        "QUALIFY_OR_KEEP_DISCOVERY_SCOPE"
    }
    or clean(x.get("section"))
    in {
        "Title",
        "Running Title"
    }
]


context_rows = [
    x for x in rows
    if clean(x.get("review_bucket"))
    == "CONTEXT_REVIEW"
]


keep_rows = [
    x for x in rows
    if clean(x.get("review_bucket"))
    == "LIKELY_KEEP_BOUNDARY_LANGUAGE"
]


if (
    len(scope_rows)
    + len(context_rows)
    + len(keep_rows)
    != 39
):

    fail(
        "QA2F grouping did not account for all 39 items."
    )


# ============================================================
# 5. Sort manuscript order
# ============================================================

SECTION_ORDER = {
    "Title": 1,
    "Running Title": 2,
    "Abstract": 3,
    "Introduction": 4,
    "Results": 5,
    "Discussion": 6,
    "Methods": 7,
    "Data Availability": 8,
    "Code Availability": 9,
    "References": 10,
    "Figure Legends": 11,
    "Tables": 12,
    "Supplementary": 13,
    "Unspecified": 99,
}


def para_num(x):

    m = re.search(
        r"\d+",
        clean(x.get("paragraph"))
    )

    return int(m.group()) if m else 999999


def sort_key(x):

    return (
        SECTION_ORDER.get(
            clean(x.get("section")),
            50
        ),
        para_num(x),
        int(
            clean(x.get("review_order"))
            or 9999
        )
    )


scope_rows.sort(key=sort_key)
context_rows.sort(key=sort_key)
keep_rows.sort(key=sort_key)


# ============================================================
# 6. Packet writer
# ============================================================

def write_packet(
    path,
    title,
    packet_rows,
    instructions
):

    lines = [
        title,
        "=" * 84,
        "",
        "STATUS:",
        "READ_ONLY_HUMAN_REVIEW_PACKET",
        "",
        f"Items: {len(packet_rows)}",
        "",
        instructions,
        "",
        "No manuscript text has been changed.",
        "",
    ]


    current_section = None


    for i, row in enumerate(
        packet_rows,
        start=1
    ):

        section = clean(
            row.get("section")
        ) or "Unspecified"


        if section != current_section:

            current_section = section

            lines.extend([
                "",
                "=" * 84,
                section.upper(),
                "=" * 84,
                ""
            ])


        lines.extend([
            f"{i:02d}. {clean(row.get('item_id'))}",
            f"Priority      : {clean(row.get('qa2f_priority'))}",
            f"Severity      : {clean(row.get('severity'))}",
            f"Risk family   : {clean(row.get('risk_family'))}",
            f"Review bucket : {clean(row.get('review_bucket'))}",
            f"Paragraph     : {clean(row.get('paragraph'))}",
            "",
            "ORIGINAL TEXT:",
            clean(row.get("text")),
            "",
            "QA2 REASON:",
            clean(row.get("reason")),
            "",
            "HUMAN REVIEW QUESTION:",
            clean(row.get("qa2f_review_question")),
            "",
            "HUMAN DECISION:",
            "[PENDING — NO AUTO-REWRITE]",
            "",
            "-" * 84,
            ""
        ])


    path.write_text(
        "\n".join(lines),
        encoding="utf-8"
    )


write_packet(
    SCOPE_TXT,
    "STEP23C-QA2F DISCOVERY-SCOPE REVIEW PACKET",
    scope_rows,
    (
        "Primary question: decide whether each conserved / reproducible / "
        "replicated statement should be retained, explicitly restricted "
        "to the two discovery cohorts, or rewritten."
    )
)


write_packet(
    CONTEXT_TXT,
    "STEP23C-QA2F CONTEXT-REVIEW PACKET",
    context_rows,
    (
        "Primary question: inspect validation/mechanism/causality wording "
        "that requires scientific judgment rather than simple scope "
        "qualification."
    )
)


write_packet(
    KEEP_TXT,
    "STEP23C-QA2F LIKELY-KEEP-BOUNDARY PACKET",
    keep_rows,
    (
        "Primary question: confirm that existing cautionary language "
        "(not validated / inconclusive / association-only / "
        "hypothesis-generating) should be preserved rather than weakened."
    )
)


# ============================================================
# 7. Summary CSV
# ============================================================

section_counts = Counter(
    clean(x.get("section"))
    or "Unspecified"
    for x in rows
)


priority_counts = Counter(
    clean(x.get("qa2f_priority"))
    for x in rows
)


summary_rows = []


for key, value in priority_counts.most_common():

    summary_rows.append({
        "summary_type":
            "priority",

        "group":
            key,

        "count":
            value
    })


for key, value in section_counts.most_common():

    summary_rows.append({
        "summary_type":
            "section",

        "group":
            key,

        "count":
            value
    })


for key, value in bucket_counts.most_common():

    summary_rows.append({
        "summary_type":
            "review_bucket",

        "group":
            key,

        "count":
            value
    })


write_csv(
    SUMMARY_CSV,
    summary_rows,
    [
        "summary_type",
        "group",
        "count"
    ]
)


# ============================================================
# 8. Confirm DOCX unchanged
# ============================================================

docx_sha_after = sha256_file(DOCX)


if docx_sha_after != docx_sha_before:

    fail(
        "CRITICAL: QA1B DOCX changed during QA2F."
    )


# ============================================================
# 9. Status
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
            "remaining_items",
        "status":
            "39"
    },

    {
        "item":
            "discovery_scope_packet",
        "status":
            str(len(scope_rows))
    },

    {
        "item":
            "context_review_packet",
        "status":
            str(len(context_rows))
    },

    {
        "item":
            "keep_boundary_packet",
        "status":
            str(len(keep_rows))
    },

    {
        "item":
            "status",
        "status":
            "STEP23C_QA2F_REVIEW_PACKETS_READY"
    },

    {
        "item":
            "next_step",
        "status":
            "HUMAN_REVIEW_SCOPE_PACKET_FIRST"
    },
]


write_csv(
    STATUS_CSV,
    status_rows,
    [
        "item",
        "status"
    ]
)


# ============================================================
# FINAL
# ============================================================

print("=" * 84)
print("STEP23C-QA2F COMPLETE")
print("=" * 84)
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
    "REMAINING ITEMS               : 39"
)

print(
    "DISCOVERY-SCOPE PACKET        :",
    len(scope_rows)
)

print(
    "CONTEXT-REVIEW PACKET         :",
    len(context_rows)
)

print(
    "KEEP-BOUNDARY PACKET          :",
    len(keep_rows)
)

print()

print("PRIORITY COUNTS:")

for key, value in priority_counts.most_common():

    print(
        f"  {key:<34s}: {value}"
    )


print()

print("SECTION COUNTS:")

for key, value in section_counts.most_common():

    print(
        f"  {key:<24s}: {value}"
    )


print()

print("OUTPUTS:")

print(
    " ",
    relative(SCOPE_TXT)
)

print(
    " ",
    relative(CONTEXT_TXT)
)

print(
    " ",
    relative(KEEP_TXT)
)

print(
    " ",
    relative(SUMMARY_CSV)
)

print()

print("STATUS:")

print(
    "STEP23C_QA2F_REVIEW_PACKETS_READY"
)

print()

print("NEXT:")

print(
    "REVIEW DISCOVERY-SCOPE PACKET FIRST."
)

print()

