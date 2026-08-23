#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP23C-QA2C
#
# SECTION-LEVEL SUMMARY OF QA2 SCIENTIFIC WORDING RISKS
#
# READ-ONLY AUDIT ONLY.
#
# DOES NOT:
#   - modify any DOCX
#   - rewrite any sentence
#   - recompute Step22A
#   - change Step22A conclusions
#   - make final editorial decisions
#
# PURPOSE:
#   Take the frozen QA2 scientific-wording audit
#   (expected: CRITICAL 0 / HIGH 34 / MEDIUM 20)
#   and organize all 54 flagged items by:
#
#     - manuscript section
#     - severity
#     - risk family
#     - trigger
#     - paragraph
#     - original wording
#     - provisional human-review bucket
#
# The "review_bucket" is triage guidance only.
# It is NOT an automatic revision decision.
# ============================================================

from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import re
import sys


PROJECT = Path.cwd()

STEP23_META = (
    PROJECT
    / "00_metadata"
    / "Step23_manuscript_update"
)

STEP23_RESULT = (
    PROJECT
    / "04_results"
    / "manuscript_update"
    / "Step23"
)

STEP23_META.mkdir(
    parents=True,
    exist_ok=True
)

STEP23_RESULT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Fixed QA1B manuscript
# ============================================================

DOCX = (
    STEP23_RESULT
    / "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"
)


# ============================================================
# Expected QA2 files
#
# We allow either meta/ or result/ location because the prior
# QA2 script may have written the audit into either directory.
# ============================================================

QA2_AUDIT_CANDIDATES = [
    STEP23_META
    / "STEP23C_QA2_SCIENTIFIC_WORDING_AUDIT.csv",

    STEP23_RESULT
    / "STEP23C_QA2_SCIENTIFIC_WORDING_AUDIT.csv",
]

QA2_REVISION_MAP_CANDIDATES = [
    STEP23_META
    / "STEP23C_QA2_REVISION_MAP.txt",

    STEP23_RESULT
    / "STEP23C_QA2_REVISION_MAP.txt",
]

QA2_STATUS_CANDIDATES = [
    STEP23_META
    / "STEP23C_QA2_STATUS.csv",

    STEP23_RESULT
    / "STEP23C_QA2_STATUS.csv",
]


# ============================================================
# Outputs
# ============================================================

SECTION_SUMMARY_OUT = (
    STEP23_META
    / "STEP23C_QA2C_SECTION_RISK_SUMMARY.csv"
)

REVIEW_QUEUE_OUT = (
    STEP23_META
    / "STEP23C_QA2C_54_ITEM_REVIEW_QUEUE.csv"
)

STATUS_OUT = (
    STEP23_META
    / "STEP23C_QA2C_STATUS.csv"
)

REPORT_OUT = (
    STEP23_RESULT
    / "STEP23C_QA2C_SECTION_REPORT.txt"
)


print()
print("=" * 82)
print("STEP23C-QA2C: SECTION-LEVEL SCIENTIFIC WORDING RISK REPORT")
print("=" * 82)
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


def clean(x):

    return (
        str(x or "")
        .replace("\r", "")
        .strip()
    )


def norm(x):

    return re.sub(
        r"\s+",
        " ",
        clean(x)
    )


def sha256_file(path):

    h = hashlib.sha256()

    with open(
        path,
        "rb"
    ) as f:

        while True:

            b = f.read(
                1024 * 1024
            )

            if not b:
                break

            h.update(
                b
            )

    return h.hexdigest()


def read_csv(path):

    with open(
        path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        return list(
            csv.DictReader(f)
        )


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
        writer.writerows(
            rows
        )


def first_existing(paths):

    for path in paths:

        if path.exists():

            return path

    return None


def relative(path):

    try:

        return str(
            path.relative_to(
                PROJECT
            )
        )

    except Exception:

        return str(path)


def find_column(
    fieldnames,
    aliases
):

    lower_map = {
        str(x).lower():
            x
        for x in fieldnames
    }

    for alias in aliases:

        if alias.lower() in lower_map:

            return lower_map[
                alias.lower()
            ]

    # loose normalized match
    def canonical(s):

        return re.sub(
            r"[^a-z0-9]",
            "",
            str(s).lower()
        )

    canon_map = {
        canonical(x):
            x
        for x in fieldnames
    }

    for alias in aliases:

        c = canonical(
            alias
        )

        if c in canon_map:

            return canon_map[
                c
            ]

    return None


# ============================================================
# 1. Verify inputs
# ============================================================

print("===== 1. VERIFY INPUTS =====")
print()


if not DOCX.exists():

    fail(
        "QA1B manuscript not found:\n"
        + str(DOCX)
    )


DOCX_SHA_BEFORE = sha256_file(
    DOCX
)


QA2_AUDIT = first_existing(
    QA2_AUDIT_CANDIDATES
)


if QA2_AUDIT is None:

    fail(
        "STEP23C_QA2_SCIENTIFIC_WORDING_AUDIT.csv not found "
        "in Step23 metadata/results directories."
    )


QA2_REVISION_MAP = first_existing(
    QA2_REVISION_MAP_CANDIDATES
)

QA2_STATUS = first_existing(
    QA2_STATUS_CANDIDATES
)


print(
    "QA1B manuscript:"
)

print(
    " ",
    relative(
        DOCX
    )
)

print(
    "QA1B SHA256:"
)

print(
    " ",
    DOCX_SHA_BEFORE
)

print()

print(
    "QA2 audit:"
)

print(
    " ",
    relative(
        QA2_AUDIT
    )
)

if QA2_REVISION_MAP:

    print(
        "QA2 revision map:"
    )

    print(
        " ",
        relative(
            QA2_REVISION_MAP
        )
    )

if QA2_STATUS:

    print(
        "QA2 status:"
    )

    print(
        " ",
        relative(
            QA2_STATUS
        )
    )

print()


# ============================================================
# 2. Load QA2 audit and detect schema
# ============================================================

print("===== 2. LOAD QA2 AUDIT =====")
print()


rows = read_csv(
    QA2_AUDIT
)


if not rows:

    fail(
        "QA2 audit CSV is empty."
    )


fieldnames = list(
    rows[0].keys()
)


severity_col = find_column(
    fieldnames,
    [
        "severity",
        "risk_level",
        "priority",
        "level"
    ]
)

section_col = find_column(
    fieldnames,
    [
        "section",
        "manuscript_section",
        "chapter"
    ]
)

paragraph_col = find_column(
    fieldnames,
    [
        "paragraph_index",
        "paragraph",
        "paragraph_no",
        "para_index",
        "para"
    ]
)

rule_col = find_column(
    fieldnames,
    [
        "rule",
        "rule_id",
        "trigger_rule",
        "audit_rule",
        "category"
    ]
)

text_col = find_column(
    fieldnames,
    [
        "original_text",
        "paragraph_text",
        "original",
        "text",
        "sentence",
        "context"
    ]
)

trigger_col = find_column(
    fieldnames,
    [
        "trigger",
        "matched_term",
        "keyword",
        "match",
        "term"
    ]
)

reason_col = find_column(
    fieldnames,
    [
        "reason",
        "rationale",
        "explanation",
        "issue"
    ]
)

suggestion_col = find_column(
    fieldnames,
    [
        "suggestion",
        "recommended_action",
        "recommended_revision",
        "revision",
        "proposed_action"
    ]
)


if severity_col is None:

    fail(
        "Could not identify severity column.\n"
        "Columns found:\n"
        + ", ".join(
            fieldnames
        )
    )


if text_col is None:

    fail(
        "Could not identify original-text column.\n"
        "Columns found:\n"
        + ", ".join(
            fieldnames
        )
    )


print(
    "Rows loaded:",
    len(
        rows
    )
)

print()

print(
    "Detected schema:"
)

for label, col in [
    (
        "severity",
        severity_col
    ),
    (
        "section",
        section_col
    ),
    (
        "paragraph",
        paragraph_col
    ),
    (
        "rule",
        rule_col
    ),
    (
        "text",
        text_col
    ),
    (
        "trigger",
        trigger_col
    ),
    (
        "reason",
        reason_col
    ),
    (
        "suggestion",
        suggestion_col
    ),
]:

    print(
        f"  {label:<12s}: "
        f"{col or '<not present>'}"
    )


print()


# ============================================================
# 3. Normalize audit rows
# ============================================================

def normalize_severity(value):

    s = clean(
        value
    ).upper()

    mapping = {
        "CRIT":
            "CRITICAL",

        "CRITICAL":
            "CRITICAL",

        "HIGH":
            "HIGH",

        "H":
            "HIGH",

        "MED":
            "MEDIUM",

        "MEDIUM":
            "MEDIUM",

        "M":
            "MEDIUM",

        "LOW":
            "LOW",

        "L":
            "LOW",
    }

    return mapping.get(
        s,
        s or "UNKNOWN"
    )


def normalize_section(value):

    s = norm(
        value
    )

    if not s:

        return "Unspecified"

    sl = s.lower()

    if "title" in sl or "running title" in sl:

        return "Title"

    if "abstract" in sl:

        return "Abstract"

    if "intro" in sl:

        return "Introduction"

    if "result" in sl:

        return "Results"

    if "discussion" in sl:

        return "Discussion"

    if "limitation" in sl:

        return "Limitations"

    if "method" in sl:

        return "Methods"

    if "data availability" in sl:

        return "Data Availability"

    if "code availability" in sl:

        return "Code Availability"

    if "supplement" in sl:

        return "Supplementary"

    if "reference" in sl:

        return "References"

    if "figure legend" in sl:

        return "Figure Legends"

    if "table" in sl:

        return "Tables"

    return s


# ============================================================
# 4. Risk-family classification
#
# This is triage only, not a scientific decision.
# ============================================================

GENERALIZATION_TERMS = [
    "conserved",
    "reproducible",
    "reproduced",
    "replicated",
    "replication",
    "uniform",
    "across independent cohorts",
    "across cohorts",
    "cross-cohort",
    "consistent across",
    "same-direction"
]

VALIDATION_TERMS = [
    "validation",
    "validated",
    "validates",
    "external validation",
    "confirm",
    "confirmed",
    "confirmation"
]

MECHANISM_TERMS = [
    "mechanism",
    "mechanistic",
    "driver",
    "causal",
    "causality",
    "master regulator",
    "upstream regulator"
]

BOUNDARY_TERMS = [
    "not validated",
    "not mechanistically validated",
    "inconclusive",
    "hypothesis-generating",
    "association-only",
    "not causal",
    "not available",
    "not evaluable",
    "cannot establish",
    "does not establish",
    "not confirmed"
]

PACKAGING_TERMS = [
    "data availability",
    "supplementary",
    "references",
    "citation needed",
    "accession",
    "omix005710"
]


def combined_lower(item):

    return " ".join([
        clean(
            item.get(
                "text"
            )
        ),
        clean(
            item.get(
                "rule"
            )
        ),
        clean(
            item.get(
                "trigger"
            )
        ),
        clean(
            item.get(
                "reason"
            )
        )
    ]).lower()


def risk_family(item):

    s = combined_lower(
        item
    )

    section = item[
        "section"
    ]

    if (
        section
        in {
            "Data Availability",
            "Supplementary",
            "References"
        }
        or
        any(
            x in s
            for x in PACKAGING_TERMS
        )
        and
        (
            "availability" in s
            or
            "supplement" in s
            or
            "reference" in s
            or
            "citation" in s
        )
    ):

        return "PACKAGING_COMPLETENESS"

    if any(
        x in s
        for x in MECHANISM_TERMS
    ):

        return "MECHANISM_CAUSALITY"

    if any(
        x in s
        for x in VALIDATION_TERMS
    ):

        return "VALIDATION_STRENGTH"

    if any(
        x in s
        for x in GENERALIZATION_TERMS
    ):

        return "CROSS_COHORT_GENERALIZATION"

    return "OTHER_SCIENTIFIC_WORDING"


def contains_explicit_boundary(text):

    t = text.lower()

    return any(
        x in t
        for x in BOUNDARY_TERMS
    )


def scope_explicit_to_discovery(text):

    t = text.lower()

    explicit_markers = [
        "two cohorts",
        "two escc",
        "gse197677",
        "gse221561",
        "discovery cohorts",
        "original analytical context",
        "within the original"
    ]

    return any(
        x in t
        for x in explicit_markers
    )


# ============================================================
# 5. Provisional review bucket
#
# IMPORTANT:
# These are NOT automatic revision decisions.
#
# Buckets:
#
# MUST_REVIEW_AND_REVISE
#   likely direct conflict / missing packaging item
#
# QUALIFY_TO_DISCOVERY_SCOPE
#   probably valid if explicitly limited to the two discovery
#   cohorts / original analytical context
#
# LIKELY_KEEP_BOUNDARY_LANGUAGE
#   already says not validated / inconclusive / association-only
#
# CONTEXT_REVIEW
#   human judgment still required
# ============================================================

def review_bucket(item):

    family = item[
        "risk_family"
    ]

    section = item[
        "section"
    ]

    text = item[
        "text"
    ]


    if family == "PACKAGING_COMPLETENESS":

        return "MUST_REVIEW_AND_REVISE"


    if family == "MECHANISM_CAUSALITY":

        if contains_explicit_boundary(
            text
        ):

            return "LIKELY_KEEP_BOUNDARY_LANGUAGE"

        return "CONTEXT_REVIEW"


    if family == "VALIDATION_STRENGTH":

        if contains_explicit_boundary(
            text
        ):

            return "LIKELY_KEEP_BOUNDARY_LANGUAGE"

        if (
            "omix005710"
            in text.lower()
        ):

            return "MUST_REVIEW_AND_REVISE"

        return "CONTEXT_REVIEW"


    if family == "CROSS_COHORT_GENERALIZATION":

        if section in {
            "Title",
            "Abstract"
        }:

            return "MUST_REVIEW_AND_REVISE"

        if scope_explicit_to_discovery(
            text
        ):

            return "QUALIFY_OR_KEEP_DISCOVERY_SCOPE"

        return "QUALIFY_TO_DISCOVERY_SCOPE"


    return "CONTEXT_REVIEW"


normalized = []


for idx, row in enumerate(
    rows,
    start=1
):

    severity = normalize_severity(
        row.get(
            severity_col,
            ""
        )
    )

    section = normalize_section(
        row.get(
            section_col,
            ""
        )
        if section_col
        else ""
    )

    paragraph = clean(
        row.get(
            paragraph_col,
            ""
        )
        if paragraph_col
        else ""
    )

    rule = clean(
        row.get(
            rule_col,
            ""
        )
        if rule_col
        else ""
    )

    text = norm(
        row.get(
            text_col,
            ""
        )
    )

    trigger = clean(
        row.get(
            trigger_col,
            ""
        )
        if trigger_col
        else ""
    )

    reason = clean(
        row.get(
            reason_col,
            ""
        )
        if reason_col
        else ""
    )

    suggestion = clean(
        row.get(
            suggestion_col,
            ""
        )
        if suggestion_col
        else ""
    )


    item = {
        "item_id":
            f"QA2C-{idx:03d}",

        "severity":
            severity,

        "section":
            section,

        "paragraph":
            paragraph,

        "rule":
            rule,

        "trigger":
            trigger,

        "text":
            text,

        "reason":
            reason,

        "qa2_suggestion":
            suggestion,
    }

    item[
        "risk_family"
    ] = risk_family(
        item
    )

    item[
        "review_bucket"
    ] = review_bucket(
        item
    )

    normalized.append(
        item
    )


# ============================================================
# 6. Verify frozen QA2 counts
# ============================================================

print("===== 3. VERIFY QA2 COUNTS =====")
print()


severity_counts = Counter(
    x[
        "severity"
    ]
    for x in normalized
)


for severity in [
    "CRITICAL",
    "HIGH",
    "MEDIUM",
    "LOW",
    "UNKNOWN"
]:

    print(
        f"{severity:<10s}:",
        severity_counts.get(
            severity,
            0
        )
    )


print(
    "TOTAL      :",
    len(
        normalized
    )
)

print()


expected = {
    "CRITICAL":
        0,

    "HIGH":
        34,

    "MEDIUM":
        20
}


count_ok = (
    len(
        normalized
    )
    == 54
    and
    severity_counts.get(
        "CRITICAL",
        0
    )
    == 0
    and
    severity_counts.get(
        "HIGH",
        0
    )
    == 34
    and
    severity_counts.get(
        "MEDIUM",
        0
    )
    == 20
)


if not count_ok:

    fail(
        "Frozen QA2 counts do not match expected "
        "CRITICAL 0 / HIGH 34 / MEDIUM 20 / TOTAL 54.\n"
        "Stopping rather than silently summarizing a changed audit."
    )


print(
    "✓ Frozen QA2 count lock confirmed."
)

print()


# ============================================================
# 7. Sort review queue in manuscript order
# ============================================================

SECTION_ORDER = {
    "Title":
        10,

    "Abstract":
        20,

    "Introduction":
        30,

    "Results":
        40,

    "Discussion":
        50,

    "Limitations":
        55,

    "Methods":
        60,

    "Data Availability":
        70,

    "Code Availability":
        75,

    "References":
        80,

    "Figure Legends":
        90,

    "Tables":
        95,

    "Supplementary":
        100,

    "Unspecified":
        999
}


SEVERITY_ORDER = {
    "CRITICAL":
        0,

    "HIGH":
        1,

    "MEDIUM":
        2,

    "LOW":
        3,

    "UNKNOWN":
        4
}


def paragraph_number(value):

    m = re.search(
        r"\d+",
        clean(
            value
        )
    )

    return (
        int(
            m.group(0)
        )
        if m
        else 999999
    )


normalized.sort(
    key=lambda x: (
        SECTION_ORDER.get(
            x[
                "section"
            ],
            500
        ),
        SEVERITY_ORDER.get(
            x[
                "severity"
            ],
            9
        ),
        paragraph_number(
            x[
                "paragraph"
            ]
        ),
        x[
            "item_id"
        ]
    )
)


# Reassign IDs after manuscript-order sorting.
for idx, item in enumerate(
    normalized,
    start=1
):

    item[
        "review_order"
    ] = idx


# ============================================================
# 8. Build section summary
# ============================================================

print("===== 4. SECTION-LEVEL SUMMARY =====")
print()


by_section = defaultdict(
    list
)


for item in normalized:

    by_section[
        item[
            "section"
        ]
    ].append(
        item
    )


section_summary_rows = []


ordered_sections = sorted(
    by_section.keys(),
    key=lambda x:
        SECTION_ORDER.get(
            x,
            500
        )
)


for section in ordered_sections:

    items = by_section[
        section
    ]

    sev = Counter(
        x[
            "severity"
        ]
        for x in items
    )

    families = Counter(
        x[
            "risk_family"
        ]
        for x in items
    )

    buckets = Counter(
        x[
            "review_bucket"
        ]
        for x in items
    )

    triggers = Counter(
        x[
            "trigger"
        ]
        for x in items
        if clean(
            x[
                "trigger"
            ]
        )
    )

    top_triggers = "; ".join(
        f"{k} ({v})"
        for k, v in triggers.most_common(
            5
        )
    )

    section_summary_rows.append({
        "section":
            section,

        "total":
            len(
                items
            ),

        "critical":
            sev.get(
                "CRITICAL",
                0
            ),

        "high":
            sev.get(
                "HIGH",
                0
            ),

        "medium":
            sev.get(
                "MEDIUM",
                0
            ),

        "cross_cohort_generalization":
            families.get(
                "CROSS_COHORT_GENERALIZATION",
                0
            ),

        "validation_strength":
            families.get(
                "VALIDATION_STRENGTH",
                0
            ),

        "mechanism_causality":
            families.get(
                "MECHANISM_CAUSALITY",
                0
            ),

        "packaging_completeness":
            families.get(
                "PACKAGING_COMPLETENESS",
                0
            ),

        "other":
            families.get(
                "OTHER_SCIENTIFIC_WORDING",
                0
            ),

        "must_review_and_revise":
            buckets.get(
                "MUST_REVIEW_AND_REVISE",
                0
            ),

        "qualify_to_discovery_scope":
            (
                buckets.get(
                    "QUALIFY_TO_DISCOVERY_SCOPE",
                    0
                )
                +
                buckets.get(
                    "QUALIFY_OR_KEEP_DISCOVERY_SCOPE",
                    0
                )
            ),

        "likely_keep_boundary_language":
            buckets.get(
                "LIKELY_KEEP_BOUNDARY_LANGUAGE",
                0
            ),

        "context_review":
            buckets.get(
                "CONTEXT_REVIEW",
                0
            ),

        "top_triggers":
            top_triggers
    })


for row in section_summary_rows:

    print(
        f"{row['section']:<22s} "
        f"TOTAL={row['total']:<3d} "
        f"HIGH={row['high']:<3d} "
        f"MEDIUM={row['medium']:<3d} "
        f"MUST={row['must_review_and_revise']:<3d} "
        f"QUALIFY={row['qualify_to_discovery_scope']:<3d} "
        f"KEEP_BOUNDARY={row['likely_keep_boundary_language']:<3d}"
    )


print()


# ============================================================
# 9. Save full 54-item review queue
# ============================================================

write_csv(
    REVIEW_QUEUE_OUT,
    normalized,
    [
        "review_order",
        "item_id",
        "severity",
        "section",
        "paragraph",
        "risk_family",
        "review_bucket",
        "trigger",
        "rule",
        "reason",
        "qa2_suggestion",
        "text"
    ]
)


write_csv(
    SECTION_SUMMARY_OUT,
    section_summary_rows,
    [
        "section",
        "total",
        "critical",
        "high",
        "medium",
        "cross_cohort_generalization",
        "validation_strength",
        "mechanism_causality",
        "packaging_completeness",
        "other",
        "must_review_and_revise",
        "qualify_to_discovery_scope",
        "likely_keep_boundary_language",
        "context_review",
        "top_triggers"
    ]
)


# ============================================================
# 10. Global summary
# ============================================================

family_counts = Counter(
    x[
        "risk_family"
    ]
    for x in normalized
)

bucket_counts = Counter(
    x[
        "review_bucket"
    ]
    for x in normalized
)


# ============================================================
# 11. Write human-readable report
# ============================================================

report = []


report.extend([
    "STEP23C-QA2C SECTION-LEVEL SCIENTIFIC WORDING RISK REPORT",
    "=" * 78,
    "",
    "STATUS:",
    "READ_ONLY_SUMMARY_COMPLETE_NO_DOCX_MODIFICATION",
    "",
    "INPUT MANUSCRIPT:",
    relative(
        DOCX
    ),
    "",
    "INPUT MANUSCRIPT SHA256:",
    DOCX_SHA_BEFORE,
    "",
    "QA2 AUDIT:",
    relative(
        QA2_AUDIT
    ),
    "",
    "FROZEN QA2 COUNTS:",
    "  CRITICAL : 0",
    "  HIGH     : 34",
    "  MEDIUM   : 20",
    "  TOTAL    : 54",
    "",
    "IMPORTANT:",
    "The review buckets below are triage guidance only.",
    "No sentence has been edited and no editorial decision has been",
    "made automatically.",
    "",
    "=" * 78,
    "GLOBAL RISK-FAMILY SUMMARY",
    "=" * 78,
    ""
])


for family in [
    "CROSS_COHORT_GENERALIZATION",
    "VALIDATION_STRENGTH",
    "MECHANISM_CAUSALITY",
    "PACKAGING_COMPLETENESS",
    "OTHER_SCIENTIFIC_WORDING"
]:

    report.append(
        f"{family:<32s}: "
        f"{family_counts.get(family, 0)}"
    )


report.extend([
    "",
    "PROVISIONAL HUMAN-REVIEW BUCKETS",
    ""
])


for bucket in [
    "MUST_REVIEW_AND_REVISE",
    "QUALIFY_TO_DISCOVERY_SCOPE",
    "QUALIFY_OR_KEEP_DISCOVERY_SCOPE",
    "LIKELY_KEEP_BOUNDARY_LANGUAGE",
    "CONTEXT_REVIEW"
]:

    report.append(
        f"{bucket:<38s}: "
        f"{bucket_counts.get(bucket, 0)}"
    )


report.extend([
    "",
    "=" * 78,
    "SECTION SUMMARY",
    "=" * 78,
    ""
])


for row in section_summary_rows:

    report.extend([
        f"[{row['section']}]",
        (
            f"  Total {row['total']} | "
            f"HIGH {row['high']} | "
            f"MEDIUM {row['medium']}"
        ),
        (
            "  Families: "
            f"generalization={row['cross_cohort_generalization']}, "
            f"validation={row['validation_strength']}, "
            f"mechanism={row['mechanism_causality']}, "
            f"packaging={row['packaging_completeness']}, "
            f"other={row['other']}"
        ),
        (
            "  Triage: "
            f"must-review/revise={row['must_review_and_revise']}, "
            f"scope-qualify={row['qualify_to_discovery_scope']}, "
            f"likely-keep-boundary={row['likely_keep_boundary_language']}, "
            f"context-review={row['context_review']}"
        ),
        (
            "  Top triggers: "
            + (
                row[
                    "top_triggers"
                ]
                or "<not available in audit schema>"
            )
        ),
        ""
    ])


report.extend([
    "=" * 78,
    "ALL 54 FLAGGED ITEMS — MANUSCRIPT ORDER",
    "=" * 78,
    ""
])


current_section = None


for item in normalized:

    if (
        item[
            "section"
        ]
        != current_section
    ):

        current_section = item[
            "section"
        ]

        report.extend([
            "",
            "-" * 78,
            current_section.upper(),
            "-" * 78,
            ""
        ])


    report.append(
        (
            f"{item['review_order']:02d}. "
            f"[{item['severity']}] "
            f"[{item['risk_family']}] "
            f"[{item['review_bucket']}]"
        )
    )

    if item[
        "paragraph"
    ]:

        report.append(
            "    Paragraph: "
            + item[
                "paragraph"
            ]
        )

    if item[
        "trigger"
    ]:

        report.append(
            "    Trigger: "
            + item[
                "trigger"
            ]
        )

    if item[
        "rule"
    ]:

        report.append(
            "    Rule: "
            + item[
                "rule"
            ]
        )

    if item[
        "reason"
    ]:

        report.append(
            "    QA2 reason: "
            + item[
                "reason"
            ]
        )

    if item[
        "qa2_suggestion"
    ]:

        report.append(
            "    QA2 suggestion: "
            + item[
                "qa2_suggestion"
            ]
        )

    report.append(
        "    Original: "
        + item[
            "text"
        ]
    )

    report.append(
        ""
    )


report.extend([
    "=" * 78,
    "DECISION BOUNDARY FOR NEXT STEP",
    "=" * 78,
    "",
    "QA2C has NOT changed the manuscript.",
    "",
    "The next editorial step should review items in this order:",
    "",
    "1. MUST_REVIEW_AND_REVISE",
    "   Direct cross-cohort conflict or missing packaging information.",
    "",
    "2. QUALIFY_TO_DISCOVERY_SCOPE / QUALIFY_OR_KEEP_DISCOVERY_SCOPE",
    "   Determine whether the statement is valid when explicitly limited",
    "   to GSE197677 + GSE221561 / the two discovery cohorts.",
    "",
    "3. LIKELY_KEEP_BOUNDARY_LANGUAGE",
    "   These often already contain appropriate language such as",
    "   'inconclusive', 'association-only', 'not validated', or",
    "   'hypothesis-generating'. They should not be weakened merely",
    "   because a keyword triggered QA2.",
    "",
    "4. CONTEXT_REVIEW",
    "   Requires human scientific judgment before any rewrite.",
    "",
    "No global find/replace should be used for conserved, reproducible,",
    "validation, mechanism, or related terms.",
    "",
    "FINAL STATUS:",
    "STEP23C_QA2C_READ_ONLY_SECTION_REPORT_COMPLETE",
    ""
])


REPORT_OUT.write_text(
    "\n".join(
        report
    ),
    encoding="utf-8"
)


# ============================================================
# 12. Confirm DOCX unchanged
# ============================================================

DOCX_SHA_AFTER = sha256_file(
    DOCX
)


if (
    DOCX_SHA_AFTER
    != DOCX_SHA_BEFORE
):

    fail(
        "CRITICAL: QA1B DOCX SHA256 changed during read-only QA2C."
    )


# ============================================================
# 13. Status file
# ============================================================

status_rows = [
    {
        "item":
            "input_manuscript",

        "status":
            relative(
                DOCX
            )
    },

    {
        "item":
            "input_SHA256_before",

        "status":
            DOCX_SHA_BEFORE
    },

    {
        "item":
            "input_SHA256_after",

        "status":
            DOCX_SHA_AFTER
    },

    {
        "item":
            "DOCX_modified",

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
            "QA2_total_flags",

        "status":
            "54"
    },

    {
        "item":
            "QA2_CRITICAL",

        "status":
            "0"
    },

    {
        "item":
            "QA2_HIGH",

        "status":
            "34"
    },

    {
        "item":
            "QA2_MEDIUM",

        "status":
            "20"
    },

    {
        "item":
            "editorial_decisions_made",

        "status":
            "NO"
    },

    {
        "item":
            "automatic_rewrites_performed",

        "status":
            "NO"
    },

    {
        "item":
            "next_step",

        "status":
            "HUMAN_REVIEW_OF_QA2C_REPORT_BEFORE_QA2B"
    },

    {
        "item":
            "status",

        "status":
            "STEP23C_QA2C_READ_ONLY_SECTION_REPORT_COMPLETE"
    },
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

print("=" * 82)
print("STEP23C-QA2C COMPLETE")
print("=" * 82)
print()

print(
    "INPUT DOCX MODIFIED           : NO"
)

print(
    "STEP22A RECOMPUTED            : NO"
)

print(
    "AUTOMATIC REWRITES            : NO"
)

print(
    "EDITORIAL DECISIONS MADE      : NO"
)

print()

print(
    "QA2 FLAGS                     : 54"
)

print(
    "CRITICAL                      : 0"
)

print(
    "HIGH                          : 34"
)

print(
    "MEDIUM                        : 20"
)

print()

print(
    "RISK FAMILIES:"
)

for family in [
    "CROSS_COHORT_GENERALIZATION",
    "VALIDATION_STRENGTH",
    "MECHANISM_CAUSALITY",
    "PACKAGING_COMPLETENESS",
    "OTHER_SCIENTIFIC_WORDING"
]:

    print(
        f"  {family:<31s}: "
        f"{family_counts.get(family, 0)}"
    )


print()

print(
    "PROVISIONAL REVIEW BUCKETS:"
)

for bucket in [
    "MUST_REVIEW_AND_REVISE",
    "QUALIFY_TO_DISCOVERY_SCOPE",
    "QUALIFY_OR_KEEP_DISCOVERY_SCOPE",
    "LIKELY_KEEP_BOUNDARY_LANGUAGE",
    "CONTEXT_REVIEW"
]:

    print(
        f"  {bucket:<37s}: "
        f"{bucket_counts.get(bucket, 0)}"
    )


print()

print(
    "SECTION SUMMARY:"
)

for row in section_summary_rows:

    print(
        f"  {row['section']:<22s} "
        f"TOTAL={row['total']:<3d} "
        f"HIGH={row['high']:<3d} "
        f"MED={row['medium']:<3d}"
    )


print()

print(
    "OUTPUTS:"
)

print(
    " ",
    relative(
        SECTION_SUMMARY_OUT
    )
)

print(
    " ",
    relative(
        REVIEW_QUEUE_OUT
    )
)

print(
    " ",
    relative(
        REPORT_OUT
    )
)

print()

print(
    "STATUS:"
)

print(
    "STEP23C_QA2C_READ_ONLY_SECTION_REPORT_COMPLETE"
)

print()

print(
    "NEXT:"
)

print(
    "REVIEW QA2C REPORT BEFORE ANY QA2B MANUSCRIPT EDIT."
)

print()
