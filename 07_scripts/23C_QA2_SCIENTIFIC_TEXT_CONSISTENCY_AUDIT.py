#!/usr/bin/env python3

# ============================================================
# STEP23C-QA2
# SCIENTIFIC TEXT CONSISTENCY AUDIT
#
# READ-ONLY.
#
# Does NOT:
#   - modify DOCX
#   - recompute Step22A
#   - alter conclusions
#
# Purpose:
#   Identify wording that became ambiguous or overstated after
#   integration of OMIX005710 directional non-replication.
# ============================================================

from pathlib import Path
from docx import Document
import csv
import hashlib
import re


PROJECT = Path.cwd()

DOCX = (
    PROJECT
    / "04_results"
    / "manuscript_update"
    / "Step23"
    / "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"
)

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

AUDIT_OUT = (
    META
    / "STEP23C_QA2_SCIENTIFIC_WORDING_AUDIT.csv"
)

REVISION_MAP_OUT = (
    RESULT
    / "STEP23C_QA2_REVISION_MAP.txt"
)

STATUS_OUT = (
    META
    / "STEP23C_QA2_STATUS.csv"
)


def norm(x):
    return re.sub(
        r"\s+",
        " ",
        str(x or "")
    ).strip()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def write_csv(path, rows, fields):
    with open(
        path,
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields
        )
        w.writeheader()
        w.writerows(rows)


if not DOCX.exists():
    raise SystemExit(
        f"Missing QA1B manuscript:\n{DOCX}"
    )


SOURCE_SHA = sha256(
    DOCX
)

doc = Document(
    DOCX
)

paragraphs = []

for i, p in enumerate(
    doc.paragraphs,
    start=1
):
    t = norm(
        p.text
    )

    if not t:
        continue

    paragraphs.append({
        "paragraph":
            i,

        "style":
            (
                p.style.name
                if p.style
                else ""
            ),

        "text":
            t
    })


# ============================================================
# Section inference based on actual manuscript ordering:
#
# front matter -> Abstract -> Introduction -> Results ->
# Discussion -> Methods -> post-method material
# ============================================================

SECTION_PATTERNS = [
    (
        "Abstract",
        r"^\s*Abstract\s*$"
    ),
    (
        "Introduction",
        r"^\s*1\s*[\.\)]?\s*Introduction\s*$"
    ),
    (
        "Results",
        r"^\s*2\s*[\.\)]?\s*Results\s*$"
    ),
    (
        "Discussion",
        r"^\s*3\s*[\.\)]?\s*Discussion\s*$"
    ),
    (
        "Methods",
        r"^\s*4\s*[\.\)]?\s*Methods\s*$"
    ),
    (
        "Data Availability",
        r"^\s*Data\s+Availability\s*$"
    ),
    (
        "References",
        r"^\s*References\s*$"
    ),
    (
        "Supplementary",
        r"^\s*Supplementary\s+Material\s+Index\s*$"
    ),
]


section_starts = []

for row in paragraphs:
    for section, pattern in SECTION_PATTERNS:
        if re.match(
            pattern,
            row["text"],
            flags=re.I
        ):
            section_starts.append(
                (
                    row["paragraph"],
                    section
                )
            )


section_starts.sort()


def section_for_paragraph(n):
    current = "FrontMatter"

    for start, section in section_starts:
        if n > start:
            current = section
        else:
            break

    return current


# ============================================================
# Audit rules
#
# Severity:
#   CRITICAL = scientific hierarchy conflict
#   HIGH     = likely misleading without qualification
#   MEDIUM   = incompleteness / manuscript packaging issue
# ============================================================

RULES = [
    {
        "id":
            "UNQUALIFIED_CONSERVED",

        "pattern":
            r"\bconserved\b",

        "severity":
            "HIGH",

        "issue":
            (
                "Unqualified 'conserved' may imply replication "
                "across all independent cohorts."
            ),

        "guidance":
            (
                "Retain only when explicitly limited to the two "
                "discovery cohorts / original analytical context."
            )
    },

    {
        "id":
            "UNQUALIFIED_REPRODUCIBLE",

        "pattern":
            r"\breproducib(?:le|ility)\b",

        "severity":
            "HIGH",

        "issue":
            (
                "Unqualified reproducibility language may conflict "
                "with OMIX directional non-replication."
            ),

        "guidance":
            (
                "Qualify as reproducible across the two discovery "
                "cohorts or within the discovery analytical context."
            )
    },

    {
        "id":
            "VALIDATED_OR_VALIDATION",

        "pattern":
            r"\bvalidat(?:e|ed|ion)\b",

        "severity":
            "MEDIUM",

        "issue":
            (
                "Validation wording requires primary-vs-secondary "
                "qualification for OMIX005710."
            ),

        "guidance":
            (
                "Ensure OMIX primary ESCC-only validation is described "
                "as unavailable, not failed or completed."
            )
    },

    {
        "id":
            "SIGNIFICANCE_LANGUAGE",

        "pattern":
            r"\b(significant|significantly|p\s*[<=>])",

        "severity":
            "CRITICAL",

        "issue":
            (
                "Formal significance claims are incompatible with the "
                "frozen small-n OMIX inferential boundary unless tied "
                "to another valid analysis."
            ),

        "guidance":
            (
                "Audit context manually; OMIX inference should emphasize "
                "direction and robustness rather than P<0.05."
            )
    },

    {
        "id":
            "UNIFORM_GENERALIZATION",

        "pattern":
            r"\b(uniform|universally|generaliz|across independent cohorts)\b",

        "severity":
            "HIGH",

        "issue":
            (
                "Generalization wording must acknowledge OMIX discordance."
            ),

        "guidance":
            (
                "State that the CORE response was not uniformly reproduced "
                "across independent cohorts."
            )
    },

    {
        "id":
            "MECHANISM_CAUSAL",

        "pattern":
            r"\b(caused|causes|mechanism|mechanistic|driver|drives)\b",

        "severity":
            "MEDIUM",

        "issue":
            (
                "Mechanistic language must not assign a cause to OMIX "
                "directional discordance."
            ),

        "guidance":
            (
                "Keep OMIX mechanism explicitly unresolved unless discussing "
                "separately labeled exploratory hypotheses."
            )
    },
]


audit_rows = []

for row in paragraphs:

    text = row[
        "text"
    ]

    for rule in RULES:

        if re.search(
            rule["pattern"],
            text,
            flags=re.I
        ):

            audit_rows.append({
                "paragraph":
                    row[
                        "paragraph"
                    ],

                "section":
                    section_for_paragraph(
                        row[
                            "paragraph"
                        ]
                    ),

                "style":
                    row[
                        "style"
                    ],

                "rule_id":
                    rule[
                        "id"
                    ],

                "severity":
                    rule[
                        "severity"
                    ],

                "issue":
                    rule[
                        "issue"
                    ],

                "guidance":
                    rule[
                        "guidance"
                    ],

                "text":
                    text
            })


# ============================================================
# Required-content packaging checks
# ============================================================

full_text = "\n".join(
    x[
        "text"
    ]
    for x in paragraphs
)


def add_packaging_check(
    rule_id,
    severity,
    passed,
    issue,
    guidance
):

    if passed:
        return

    audit_rows.append({
        "paragraph":
            "",

        "section":
            "Packaging",

        "style":
            "",

        "rule_id":
            rule_id,

        "severity":
            severity,

        "issue":
            issue,

        "guidance":
            guidance,

        "text":
            ""
    })


# Data Availability should include OMIX.
data_rows = [
    x
    for x in paragraphs
    if section_for_paragraph(
        x[
            "paragraph"
        ]
    )
    ==
    "Data Availability"
]


data_text = " ".join(
    x[
        "text"
    ]
    for x in data_rows
)


add_packaging_check(
    "DATA_AVAILABILITY_OMIX_MISSING",
    "HIGH",
    (
        "OMIX005710"
        in data_text
    ),
    (
        "Data Availability does not currently include OMIX005710."
    ),
    (
        "Add the public OMIX/GSA-Human dataset accession/source in the "
        "final citation/data-availability pass."
    )
)


# Supplementary index should mention external validation.
supp_rows = [
    x
    for x in paragraphs
    if section_for_paragraph(
        x[
            "paragraph"
        ]
    )
    ==
    "Supplementary"
]


supp_text = " ".join(
    x[
        "text"
    ]
    for x in supp_rows
)


add_packaging_check(
    "SUPPLEMENTARY_STEP22A_MISSING",
    "HIGH",
    (
        "OMIX005710"
        in supp_text
        or
        "external validation"
        in supp_text.lower()
    ),
    (
        "Supplementary Material Index does not yet list the "
        "Step22A external-validation audit material."
    ),
    (
        "Add a dedicated supplementary table/material entry covering "
        "OMIX technical eligibility, 22 histology scenarios, CORE "
        "directional summary, and 105-run LOPO audit."
    )
)


# References remain placeholder.
add_packaging_check(
    "REFERENCES_PLACEHOLDER",
    "HIGH",
    (
        "[CITATION NEEDED]"
        not in full_text
        and
        "References require human completion"
        not in full_text
    ),
    (
        "References remain explicitly incomplete."
    ),
    (
        "Do not finalize the manuscript until citations for the "
        "source cohorts and OMIX005710 are completed."
    )
)


write_csv(
    AUDIT_OUT,
    audit_rows,
    [
        "paragraph",
        "section",
        "style",
        "rule_id",
        "severity",
        "issue",
        "guidance",
        "text"
    ]
)


# ============================================================
# Summary + revision priorities
# ============================================================

severity_order = {
    "CRITICAL": 0,
    "HIGH": 1,
    "MEDIUM": 2
}


counts = {}

for row in audit_rows:
    key = row[
        "severity"
    ]

    counts[
        key
    ] = (
        counts.get(
            key,
            0
        )
        + 1
    )


revision_map = f"""STEP23C-QA2 SCIENTIFIC TEXT CONSISTENCY AUDIT
================================================================

SOURCE:
{DOCX}

SOURCE SHA256:
{SOURCE_SHA}

DOCX MODIFIED:
NO

STEP22A RECOMPUTED:
NO

AUDIT FLAGS:
CRITICAL = {counts.get("CRITICAL", 0)}
HIGH     = {counts.get("HIGH", 0)}
MEDIUM   = {counts.get("MEDIUM", 0)}

================================================================
FROZEN SCIENTIFIC HIERARCHY
================================================================

DISCOVERY:
GSE197677 and GSE221561 support the frozen CORE framework.

OMIX005710 PRIMARY ESCC-ONLY:
UNAVAILABLE because the EASC patient identity remains unresolved.

OMIX005710 SECONDARY:
22/22 histology scenarios show direction opposite to the frozen
discovery-cohort direction.

TECHNICAL ROBUSTNESS:
No simple implementation artifact detected; 105/105 paired LOPO runs
retain the negative CORE2 direction.

BIOLOGICAL CAUSE:
NOT ESTABLISHED.

================================================================
REVISION PRIORITIES FOR NEXT EDITING STEP
================================================================

1. TITLE / RUNNING TITLE
   Qualify or remove unbounded 'conserved' language so the title does
   not imply three-cohort directional replication.

2. ABSTRACT
   Preserve the two-discovery-cohort result, but recalibrate the final
   interpretation to acknowledge OMIX directional discordance.

3. INTRODUCTION
   Clarify that the reproducibility objective concerned the discovery
   cohorts and that OMIX assessed external generalizability.

4. RESULTS
   Keep 'conserved' only when explicitly referring to GSE197677 and
   GSE221561.

5. DISCUSSION / CONCLUSION
   Replace project-wide claims of stable/conserved response with
   discovery-context wording plus independent-cohort heterogeneity.

6. DATA AVAILABILITY
   Add OMIX005710.

7. SUPPLEMENTARY MATERIAL INDEX
   Add Step22A external-validation audit material.

8. REFERENCES
   Complete real citations before finalization.

================================================================
NEXT
================================================================

STEP23C-QA2B:
CREATE A NEW SCIENTIFICALLY RECALIBRATED WORKING COPY.

Do not overwrite QA1B.
"""


REVISION_MAP_OUT.write_text(
    revision_map,
    encoding="utf-8"
)


status_rows = [
    {
        "item":
            "source_docx",
        "status":
            str(
                DOCX.relative_to(
                    PROJECT
                )
            )
    },

    {
        "item":
            "source_SHA256",
        "status":
            SOURCE_SHA
    },

    {
        "item":
            "visual_QA1B",
        "status":
            "PASS"
    },

    {
        "item":
            "docx_modified",
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
            "critical_flags",
        "status":
            str(
                counts.get(
                    "CRITICAL",
                    0
                )
            )
    },

    {
        "item":
            "high_flags",
        "status":
            str(
                counts.get(
                    "HIGH",
                    0
                )
            )
    },

    {
        "item":
            "medium_flags",
        "status":
            str(
                counts.get(
                    "MEDIUM",
                    0
                )
            )
    },

    {
        "item":
            "next_step",
        "status":
            "STEP23C_QA2B_SCIENTIFIC_RECALIBRATION_EDIT"
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


print("=" * 76)
print("STEP23C-QA2 COMPLETE")
print("=" * 76)
print()

print(
    "VISUAL QA1B                : PASS"
)

print(
    "DOCX MODIFIED              : NO"
)

print(
    "STEP22A RECOMPUTED          : NO"
)

print()

print(
    "AUDIT FLAGS:"
)

print(
    "  CRITICAL                  :",
    counts.get(
        "CRITICAL",
        0
    )
)

print(
    "  HIGH                      :",
    counts.get(
        "HIGH",
        0
    )
)

print(
    "  MEDIUM                    :",
    counts.get(
        "MEDIUM",
        0
    )
)

print()

print(
    "AUDIT FILE:"
)

print(
    " ",
    AUDIT_OUT.relative_to(
        PROJECT
    )
)

print()

print(
    "REVISION MAP:"
)

print(
    " ",
    REVISION_MAP_OUT.relative_to(
        PROJECT
    )
)

print()

print(
    "NEXT:"
)

print(
    "STEP23C_QA2B_SCIENTIFIC_RECALIBRATION_EDIT"
)

print()
