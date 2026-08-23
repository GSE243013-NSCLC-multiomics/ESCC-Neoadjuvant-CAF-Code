#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP23B
#
# MANUSCRIPT EVIDENCE-MAP FREEZE
#
# PURPOSE
#   Read the frozen Step21E manuscript and frozen Step22A
#   conclusions, then build a section-by-section change map.
#
# STRICT RULES
#   - DO NOT modify the DOCX
#   - DO NOT recompute Step22A
#   - DO NOT change Step22A conclusions
#   - DO NOT create new biological claims
#   - DO NOT promote secondary OMIX analysis to primary
#
# OUTPUTS
#   1. STEP23B_MANUSCRIPT_SECTION_AUDIT.csv
#   2. STEP23B_MANUSCRIPT_KEYWORD_ANCHORS.csv
#   3. STEP23B_FROZEN_EVIDENCE_MAP.csv
#   4. STEP23B_FROZEN_REPORTING_RULES.txt
#   5. STEP23B_PROPOSED_INSERTION_TEXT.txt
#   6. STEP23B_STATUS.csv
#
# This step is READ-ONLY with respect to the manuscript.
# ============================================================

from pathlib import Path
import csv
import hashlib
import re
import zipfile
import xml.etree.ElementTree as ET
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

STEP22_META = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

STEP22_RESULT = (
    PROJECT
    / "04_results"
    / "external_validation"
    / "Step22A"
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
# Inputs
# ============================================================

STEP23A_STATUS = (
    STEP23_META
    / "STEP23A_PREFLIGHT_STATUS.csv"
)

STEP23A_CANDIDATES = (
    STEP23_META
    / "STEP23A_MANUSCRIPT_CANDIDATE_INVENTORY.csv"
)

STEP22_LOCK = (
    STEP22_META
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_LOCK.csv"
)

STEP22_DOSSIER = (
    STEP22_META
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_DOSSIER.txt"
)

STEP22_REPORTING = (
    STEP22_RESULT
    / "STEP22A_07D_REPORTING_LANGUAGE.txt"
)

STEP22_README = (
    STEP22_META
    / "STEP22A_FINAL_PROJECT_README.md"
)

STEP22_CLOSEOUT = (
    STEP22_META
    / "STEP22A_FINAL_CLOSEOUT.txt"
)

# ============================================================
# Outputs
# ============================================================

SECTION_AUDIT_OUT = (
    STEP23_META
    / "STEP23B_MANUSCRIPT_SECTION_AUDIT.csv"
)

ANCHOR_OUT = (
    STEP23_META
    / "STEP23B_MANUSCRIPT_KEYWORD_ANCHORS.csv"
)

EVIDENCE_MAP_OUT = (
    STEP23_META
    / "STEP23B_FROZEN_EVIDENCE_MAP.csv"
)

RULES_OUT = (
    STEP23_RESULT
    / "STEP23B_FROZEN_REPORTING_RULES.txt"
)

TEXT_OUT = (
    STEP23_RESULT
    / "STEP23B_PROPOSED_INSERTION_TEXT.txt"
)

STATUS_OUT = (
    STEP23_META
    / "STEP23B_STATUS.csv"
)


print()
print("=" * 80)
print("STEP23B: MANUSCRIPT EVIDENCE-MAP FREEZE")
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


def clean(x):

    return (
        str(x or "")
        .replace("\r", "")
        .strip()
    )


def read_csv(path):

    if not path.exists():

        fail(
            "Missing required file:\n"
            + str(path)
        )

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
        writer.writerows(rows)


def sha256_file(path):

    h = hashlib.sha256()

    with open(
        path,
        "rb"
    ) as f:

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
            path.relative_to(
                PROJECT
            )
        )

    except Exception:
        return str(path)


def truth(x):

    return clean(x).upper() in {
        "YES",
        "TRUE",
        "T",
        "1",
        "PASS"
    }


# ============================================================
# DOCX read-only parser
# stdlib only: no document modification
# ============================================================

W_NS = (
    "http://schemas.openxmlformats.org/"
    "wordprocessingml/2006/main"
)

NS = {
    "w": W_NS
}


def qn(tag):

    return (
        "{"
        + W_NS
        + "}"
        + tag
    )


def read_docx_paragraphs(path):

    if not path.exists():

        fail(
            "Manuscript DOCX missing:\n"
            + str(path)
        )

    paragraphs = []

    with zipfile.ZipFile(
        path,
        "r"
    ) as z:

        if "word/document.xml" not in z.namelist():

            fail(
                "Invalid DOCX: word/document.xml missing."
            )

        xml_bytes = z.read(
            "word/document.xml"
        )

    root = ET.fromstring(
        xml_bytes
    )

    body = root.find(
        "w:body",
        NS
    )

    if body is None:

        fail(
            "DOCX body not found."
        )

    paragraph_index = 0

    for child in body:

        if child.tag != qn("p"):
            continue

        paragraph_index += 1

        texts = []

        for node in child.iter():

            if node.tag == qn("t"):

                texts.append(
                    node.text or ""
                )

            elif node.tag == qn("tab"):

                texts.append(
                    "\t"
                )

            elif node.tag in {
                qn("br"),
                qn("cr")
            }:

                texts.append(
                    "\n"
                )

        text = "".join(
            texts
        ).strip()

        ppr = child.find(
            "w:pPr",
            NS
        )

        style = ""

        if ppr is not None:

            pstyle = ppr.find(
                "w:pStyle",
                NS
            )

            if pstyle is not None:

                style = (
                    pstyle.attrib.get(
                        qn("val"),
                        ""
                    )
                )

        paragraphs.append({
            "paragraph_index":
                paragraph_index,

            "style":
                style,

            "text":
                text
        })

    return paragraphs


# ============================================================
# 1. Verify Step23A
# ============================================================

print("===== 1. VERIFY STEP23A =====")
print()

status23a = read_csv(
    STEP23A_STATUS
)

status23a_map = {
    clean(
        x.get(
            "item"
        )
    ):
    clean(
        x.get(
            "status"
        )
    )
    for x in status23a
}


if (
    status23a_map.get(
        "Step22A_archived_and_closed"
    )
    != "PASS"
):

    fail(
        "Step23A does not confirm Step22A closure."
    )


if (
    status23a_map.get(
        "original_manuscript_modified"
    )
    != "NO"
):

    fail(
        "Original manuscript was unexpectedly modified."
    )


print(
    "Step22A archived/closed : PASS"
)

print(
    "Original manuscript edited: NO"
)

print()


# ============================================================
# 2. Resolve top manuscript candidate
# ============================================================

print("===== 2. RESOLVE MANUSCRIPT SOURCE =====")
print()

candidate_rows = read_csv(
    STEP23A_CANDIDATES
)

docx_candidates = [
    x
    for x in candidate_rows
    if clean(
        x.get(
            "suffix"
        )
    ).lower()
    == ".docx"
]


if not docx_candidates:

    fail(
        "No DOCX manuscript candidate found."
    )


docx_candidates.sort(
    key=lambda x:
        int(
            clean(
                x.get(
                    "rank"
                )
            )
            or "999999"
        )
)


top = docx_candidates[0]

MANUSCRIPT = (
    PROJECT
    / clean(
        top.get(
            "path"
        )
    )
)


expected_sha = clean(
    top.get(
        "sha256"
    )
)

observed_sha = sha256_file(
    MANUSCRIPT
)


print(
    "Manuscript:"
)

print(
    " ",
    relative(
        MANUSCRIPT
    )
)

print()

print(
    "Recorded SHA256:",
    expected_sha
)

print(
    "Current  SHA256:",
    observed_sha
)

print()


if expected_sha != observed_sha:

    fail(
        "MANUSCRIPT_CHANGED_SINCE_STEP23A"
    )


print(
    "✓ Manuscript source unchanged since Step23A."
)

print()


# ============================================================
# 3. Verify Step22A final frozen state
# ============================================================

print("===== 3. VERIFY STEP22A FROZEN CONCLUSION =====")
print()

for path in [
    STEP22_LOCK,
    STEP22_DOSSIER,
    STEP22_REPORTING,
    STEP22_README,
    STEP22_CLOSEOUT
]:

    if not path.exists():

        fail(
            "Missing Step22A frozen artifact:\n"
            + str(path)
        )


lock_rows = read_csv(
    STEP22_LOCK
)

lock = {
    clean(
        x.get(
            "item"
        )
    ):
    clean(
        x.get(
            "value"
        )
    )
    for x in lock_rows
}


expected_lock = {
    "primary_external_validation":
        "PRIMARY_ESCC_ONLY_VALIDATION_UNAVAILABLE_HISTOLOGY_UNRESOLVED",

    "secondary_histology_sensitivity":
        "SECONDARY_HISTOLOGY_SENSITIVITY_DIRECTION_ROBUSTLY_OPPOSITE",

    "CORE2_mean_negative_scenarios":
        "22_OF_22",

    "CORE4_mean_negative_scenarios":
        "22_OF_22",

    "paired_LOPO_runs":
        "105",

    "paired_LOPO_mean_negative":
        "105_OF_105",

    "posthoc_CAF_subtype_explanation_performed":
        "NO"
}


for key, expected in expected_lock.items():

    observed = lock.get(
        key,
        "<missing>"
    )

    print(
        f"{key:<44s}: "
        f"{observed}"
    )

    if observed != expected:

        fail(
            "Step22A frozen lock mismatch:\n"
            + key
        )


print()
print(
    "✓ Step22A evidence state confirmed."
)
print()


# ============================================================
# 4. Parse manuscript paragraphs
# ============================================================

print("===== 4. READ MANUSCRIPT STRUCTURE =====")
print()

paragraphs = read_docx_paragraphs(
    MANUSCRIPT
)


nonempty = [
    x
    for x in paragraphs
    if clean(
        x.get(
            "text"
        )
    )
]


print(
    "Total paragraphs:",
    len(
        paragraphs
    )
)

print(
    "Non-empty paragraphs:",
    len(
        nonempty
    )
)

print()


# ============================================================
# 5. Detect section headings
# ============================================================

print("===== 5. DETECT MANUSCRIPT SECTIONS =====")
print()


SECTION_PATTERNS = [
    (
        "Abstract",
        re.compile(
            r"^\s*abstract\s*$",
            re.I
        )
    ),

    (
        "Introduction",
        re.compile(
            r"^\s*introduction\s*$",
            re.I
        )
    ),

    (
        "Methods",
        re.compile(
            r"^\s*(materials\s+and\s+methods|methods)\s*$",
            re.I
        )
    ),

    (
        "Results",
        re.compile(
            r"^\s*results\s*$",
            re.I
        )
    ),

    (
        "Discussion",
        re.compile(
            r"^\s*discussion\s*$",
            re.I
        )
    ),

    (
        "Limitations",
        re.compile(
            r"^\s*(limitations|study\s+limitations)\s*$",
            re.I
        )
    ),

    (
        "Conclusion",
        re.compile(
            r"^\s*(conclusion|conclusions)\s*$",
            re.I
        )
    ),

    (
        "References",
        re.compile(
            r"^\s*references\s*$",
            re.I
        )
    )
]


section_hits = []


for p in paragraphs:

    text = clean(
        p.get(
            "text"
        )
    )

    if not text:
        continue

    for section_name, pattern in SECTION_PATTERNS:

        if pattern.match(
            text
        ):

            section_hits.append({
                "section":
                    section_name,

                "paragraph_index":
                    p[
                        "paragraph_index"
                    ],

                "style":
                    p[
                        "style"
                    ],

                "heading_text":
                    text
            })


for row in section_hits:

    print(
        f"{row['section']:<14s} "
        f"paragraph={row['paragraph_index']:<4d} "
        f"style={row['style'] or '<none>'}"
    )


print()


# ============================================================
# 6. Assign paragraphs to detected sections
# ============================================================

section_hits_sorted = sorted(
    section_hits,
    key=lambda x:
        x[
            "paragraph_index"
        ]
)


def section_for_paragraph(index):

    current = "FrontMatter"

    for hit in section_hits_sorted:

        if index > hit[
            "paragraph_index"
        ]:

            current = hit[
                "section"
            ]

        else:
            break

    return current


section_counts = {}

for p in paragraphs:

    section = section_for_paragraph(
        p[
            "paragraph_index"
        ]
    )

    if clean(
        p.get(
            "text"
        )
    ):

        section_counts[
            section
        ] = (
            section_counts.get(
                section,
                0
            )
            + 1
        )


section_audit_rows = []

expected_sections = [
    "Abstract",
    "Introduction",
    "Methods",
    "Results",
    "Discussion",
    "Limitations",
    "Conclusion"
]


for section in expected_sections:

    heading_rows = [
        x
        for x in section_hits
        if x[
            "section"
        ]
        == section
    ]

    found = (
        len(
            heading_rows
        )
        > 0
    )

    section_audit_rows.append({
        "section":
            section,

        "heading_detected":
            "YES"
            if found
            else "NO",

        "heading_paragraph_index":
            (
                heading_rows[0][
                    "paragraph_index"
                ]
                if found
                else ""
            ),

        "body_nonempty_paragraphs":
            section_counts.get(
                section,
                0
            )
    })


write_csv(
    SECTION_AUDIT_OUT,
    section_audit_rows,
    [
        "section",
        "heading_detected",
        "heading_paragraph_index",
        "body_nonempty_paragraphs"
    ]
)


# ============================================================
# 7. Keyword / anchor audit
# ============================================================

print("===== 6. KEYWORD / ANCHOR AUDIT =====")
print()


KEYWORDS = [
    "GSE197677",
    "GSE221561",
    "OMIX005710",
    "external validation",
    "validation cohort",
    "CORE2",
    "CORE4",
    "fibroblast",
    "CAF",
    "neoadjuvant",
    "treatment",
    "before",
    "after",
    "pseudobulk",
    "sensitivity",
    "EASC",
    "adenosquamous"
]


anchor_rows = []


for keyword in KEYWORDS:

    regex = re.compile(
        re.escape(
            keyword
        ),
        re.I
    )

    hits = []

    for p in paragraphs:

        text = clean(
            p.get(
                "text"
            )
        )

        if not text:
            continue

        if regex.search(
            text
        ):

            hits.append(
                p
            )

    print(
        f"{keyword:<22s}: "
        f"{len(hits)}"
    )

    for p in hits:

        context_start = max(
            1,
            p[
                "paragraph_index"
            ]
            - 1
        )

        context_end = min(
            len(
                paragraphs
            ),
            p[
                "paragraph_index"
            ]
            + 1
        )

        context_text = []

        for q in paragraphs:

            if (
                context_start
                <=
                q[
                    "paragraph_index"
                ]
                <=
                context_end
            ):

                t = clean(
                    q.get(
                        "text"
                    )
                )

                if t:

                    context_text.append(
                        t
                    )

        anchor_rows.append({
            "keyword":
                keyword,

            "paragraph_index":
                p[
                    "paragraph_index"
                ],

            "section":
                section_for_paragraph(
                    p[
                        "paragraph_index"
                    ]
                ),

            "style":
                p[
                    "style"
                ],

            "paragraph_text":
                clean(
                    p.get(
                        "text"
                    )
                ),

            "local_context":
                " || ".join(
                    context_text
                )
        })


write_csv(
    ANCHOR_OUT,
    anchor_rows,
    [
        "keyword",
        "paragraph_index",
        "section",
        "style",
        "paragraph_text",
        "local_context"
    ]
)


print()


# ============================================================
# 8. Freeze manuscript evidence map
# ============================================================

print("===== 7. FREEZE SECTION-BY-SECTION EVIDENCE MAP =====")
print()


evidence_map = [
    {
        "section":
            "Abstract",

        "priority":
            "HIGH",

        "action":
            "UPDATE",

        "evidence_level":
            "SECONDARY_EXTERNAL_VALIDATION",

        "required_content":
            (
                "Add one concise sentence stating that definitive "
                "ESCC-only validation in OMIX005710 was unavailable "
                "because the single EASC case could not be resolved "
                "at patient level, while the pre-specified 22-scenario "
                "sensitivity analysis showed a consistently opposite "
                "treatment-associated CORE direction."
            ),

        "must_not_claim":
            (
                "Do not say OMIX005710 provided definitive ESCC-only "
                "validation; do not say statistically significant "
                "reversal; do not claim opposite biological mechanism."
            )
    },

    {
        "section":
            "Methods",

        "priority":
            "HIGH",

        "action":
            "ADD_SUBSECTION",

        "evidence_level":
            "METHODS",

        "required_content":
            (
                "Describe OMIX005710 as an external neoadjuvant "
                "single-cell cohort; state technical fibroblast QC "
                "(>=20 cells), P15_Before exclusion, raw-count "
                "pseudobulk construction, log2(CPM+1), scenario-specific "
                "gene-wise z-scoring, frozen CORE2/CORE4 formulas, "
                "unresolved EASC histology gate, 22 hypothetical EASC "
                "sensitivity scenarios, and paired LOPO influence audit."
            ),

        "must_not_claim":
            (
                "Do not describe the 22 scenarios as independent cohorts; "
                "do not describe LOPO as a significance test; do not "
                "describe P15 as biologically excluded."
            )
    },

    {
        "section":
            "Results",

        "priority":
            "CRITICAL",

        "action":
            "ADD_SUBSECTION",

        "evidence_level":
            "SECONDARY_EXTERNAL_VALIDATION",

        "required_content":
            (
                "State primary ESCC-only validation was unavailable. "
                "Then report the pre-specified secondary sensitivity: "
                "CORE2_STRONG mean and median After-Before changes were "
                "negative in 22/22 scenarios; CORE4 mean change was "
                "negative in 22/22; CDKN1B and GSN both contributed "
                "negatively in 22/22; BGN/TIMP1 was negative in 18/22 "
                "and positive in 4/22; independent technical "
                "recalculation reproduced scores within floating-point "
                "precision; 105/105 paired LOPO runs retained negative "
                "CORE2 direction."
            ),

        "must_not_claim":
            (
                "Do not call this a successful replication; do not call "
                "it definitive primary external validation; do not "
                "describe the directional result as statistically "
                "significant."
            )
    },

    {
        "section":
            "Discussion",

        "priority":
            "CRITICAL",

        "action":
            "UPDATE",

        "evidence_level":
            "INTERPRETATION",

        "required_content":
            (
                "Interpret OMIX005710 as a robust secondary directional "
                "non-replication / discordance. Emphasize that the "
                "treatment-associated CORE response was not uniformly "
                "replicated across independent cohorts and may be "
                "context-dependent. Note that technical audits argue "
                "against simple sign, pairing, reference-set, or "
                "single-patient explanations."
            ),

        "must_not_claim":
            (
                "Do not infer a reverse biological mechanism. Do not "
                "attribute the discordance to a CAF subtype, treatment "
                "regimen, pathology response, or platform without a "
                "separate exploratory analysis."
            )
    },

    {
        "section":
            "Limitations",

        "priority":
            "CRITICAL",

        "action":
            "ADD_OR_UPDATE",

        "evidence_level":
            "LIMITATION",

        "required_content":
            (
                "Explicitly state that the EASC patient's identity was "
                "not resolvable from public patient-level metadata, "
                "preventing a definitive ESCC-only primary analysis. "
                "Also note the small paired sample size (5 technically "
                "evaluable pairs before hypothetical EASC removal; "
                "4 in scenarios excluding a paired patient), limiting "
                "exact-test resolution."
            ),

        "must_not_claim":
            (
                "Do not imply missing histology was imputed; do not "
                "present lack of P<0.05 as evidence of no effect."
            )
    },

    {
        "section":
            "Conclusion",

        "priority":
            "MEDIUM",

        "action":
            "RECALIBRATE",

        "evidence_level":
            "SYNTHESIS",

        "required_content":
            (
                "Ensure the overall conclusion no longer implies uniform "
                "external reproducibility of the treatment-associated "
                "CORE direction. Preserve the discovery-cohort findings "
                "while acknowledging independent-cohort heterogeneity."
            ),

        "must_not_claim":
            (
                "Do not state that CORE response is universally "
                "reproducible across cohorts."
            )
    },

    {
        "section":
            "Supplementary",

        "priority":
            "HIGH",

        "action":
            "ADD",

        "evidence_level":
            "AUDIT_TRAIL",

        "required_content":
            (
                "Provide technical eligibility/QC table, 22-scenario "
                "histology sensitivity summary, CORE2/CORE4 directional "
                "summary, component decomposition, 105-run LOPO summary, "
                "and governance/provenance statement."
            ),

        "must_not_claim":
            (
                "Do not create 22 pseudo-replication P values or "
                "meta-analyze the scenarios."
            )
    }
]


write_csv(
    EVIDENCE_MAP_OUT,
    evidence_map,
    [
        "section",
        "priority",
        "action",
        "evidence_level",
        "required_content",
        "must_not_claim"
    ]
)


for row in evidence_map:

    print(
        f"{row['section']:<14s} "
        f"{row['priority']:<8s} "
        f"{row['action']}"
    )


print()


# ============================================================
# 9. Freeze precise reporting rules
# ============================================================

print("===== 8. FREEZE REPORTING RULES =====")
print()


rules_text = """STEP23B FROZEN REPORTING RULES
======================================================================

STATUS:
FROZEN_BEFORE_MANUSCRIPT_EDITING

SOURCE MANUSCRIPT:
{manuscript}

SOURCE SHA256:
{sha}

----------------------------------------------------------------------
A. EVIDENCE HIERARCHY
----------------------------------------------------------------------

1. GSE197677 / GSE221561
   Discovery evidence remains frozen.

2. OMIX005710 primary ESCC-only external validation
   UNAVAILABLE.

   Reason:
   the single EASC patient could not be mapped to an individual
   patient using the available public metadata.

3. OMIX005710 secondary histology sensitivity
   PRE-SPECIFIED / FROZEN BEFORE CORE SCORING.

   22 hypothetical EASC assignments were evaluated.

   CORE2_STRONG:
   negative treatment-related mean delta in 22/22 scenarios.

   CORE2_STRONG:
   negative treatment-related median delta in 22/22 scenarios.

   CORE4:
   negative mean delta in 22/22 scenarios.

   CDKN1B contribution:
   negative in 22/22 scenarios.

   GSN contribution:
   negative in 22/22 scenarios.

   BGN/TIMP1 arm:
   negative in 18/22 scenarios;
   positive in 4/22 scenarios.

4. Technical replication
   Independent reconstruction reproduced Step22A within
   floating-point precision.

5. Paired influence audit
   105/105 full LOPO runs retained negative CORE2 direction.

6. Biological explanation
   NOT ESTABLISHED.

----------------------------------------------------------------------
B. WORDING THAT IS ALLOWED
----------------------------------------------------------------------

Preferred terms:

- "pre-specified secondary histology-sensitivity analysis"
- "directional non-replication"
- "directional discordance"
- "opposite to the frozen discovery-cohort replication direction"
- "robust to hypothetical EASC assignment"
- "robust to paired leave-one-patient-out analysis"
- "primary ESCC-only validation was unavailable"
- "may be context-dependent"
- "not uniformly replicated across independent cohorts"

----------------------------------------------------------------------
C. WORDING THAT IS NOT ALLOWED
----------------------------------------------------------------------

Do NOT write:

- "OMIX005710 definitively validated the opposite effect"
- "OMIX005710 proved the opposite biological mechanism"
- "the ESCC-only external validation failed"
- "CORE2 significantly decreased after treatment"
  unless an explicitly valid inferential analysis supports that exact
  statement; the frozen exact sign-flip resolution at n=5/n=4 cannot
  attain a two-sided P<0.05.
- "the third cohort confirmed the findings"
- "CORE4 rescued the external validation"
- "a specific CAF subtype caused the discordance"
- "EASC was identified as Pxx"
- "P15 was biologically excluded"

----------------------------------------------------------------------
D. PRIMARY VS SECONDARY BOUNDARY
----------------------------------------------------------------------

Primary:
  unavailable.

Secondary:
  complete and robustly opposite.

The secondary analysis must never be relabeled as primary.

The 22 scenarios are a sensitivity envelope, not 22 independent
datasets, cohorts, or replications.

----------------------------------------------------------------------
E. STATISTICAL BOUNDARY
----------------------------------------------------------------------

Paired technical candidates:
  P6, P9, P16, P18, P19.

Depending on hypothetical EASC identity:
  paired n = 5 or n = 4.

Exact two-sided sign-flip minimum possible P:
  n=5 -> 0.0625
  n=4 -> 0.125

Therefore inferential emphasis should be on:
  direction,
  effect size,
  patient-level consistency,
  scenario robustness,
  and LOPO robustness,

not on a claim of P<0.05.

----------------------------------------------------------------------
F. MECHANISM BOUNDARY
----------------------------------------------------------------------

No confirmatory biological mechanism was tested for the OMIX
discordance.

Any future analysis of:
  CAF subtypes,
  treatment regimen,
  response class,
  platform effects,
  sampling time,
  cell composition,
  tumor microenvironment context

must be explicitly labeled:

EXPLORATORY / HYPOTHESIS-GENERATING.

----------------------------------------------------------------------
G. MANUSCRIPT EDITING RULE
----------------------------------------------------------------------

The Step21E manuscript source must never be overwritten.

Step23C must create a NEW working copy before any content edit.

Step22A archival files must remain untouched.

======================================================================
""".format(
    manuscript=relative(
        MANUSCRIPT
    ),
    sha=observed_sha
)


RULES_OUT.write_text(
    rules_text,
    encoding="utf-8"
)


# ============================================================
# 10. Freeze proposed insertion text
#
# Draft only. Not yet inserted into DOCX.
# ============================================================

print("===== 9. FREEZE PROPOSED INSERTION TEXT =====")
print()


proposed_text = """STEP23B PROPOSED MANUSCRIPT INSERTION TEXT
======================================================================

IMPORTANT:
These paragraphs are DRAFT INSERTIONS ONLY.
They have NOT been inserted into the manuscript.
Step23C will create a new working copy before editing.

======================================================================
1. ABSTRACT — RESULTS SENTENCE
======================================================================

In an independent neoadjuvant single-cell cohort (OMIX005710),
definitive ESCC-only validation was precluded because the single
patient with esophageal adenosquamous carcinoma could not be
identified from the available patient-level metadata. In a
pre-specified secondary histology-sensitivity analysis spanning all
22 possible EASC assignments, the treatment-associated CORE
direction was consistently opposite to that observed in the
discovery cohorts.

======================================================================
2. METHODS — EXTERNAL VALIDATION
======================================================================

External validation was additionally assessed using the publicly
available OMIX005710 neoadjuvant single-cell RNA-sequencing cohort.
Fibroblast pseudobulk profiles were constructed from raw counts after
broad cell-type annotation using a marker framework independent of
the CORE genes. A pre-specified minimum of 20 fibroblast cells per
patient-timepoint was required for technical eligibility; one sample
(P15 before treatment, 18 fibroblasts) therefore remained in raw
provenance but was excluded from the technical reference and paired
analysis. Raw pseudobulk counts were transformed to library-size CPM
and log2(CPM+1), followed by gene-wise z-scoring within the eligible
OMIX tumor reference. CORE2_STRONG and CORE4 were calculated using
the previously frozen formulas without refitting gene weights.

The source cohort contained one patient with esophageal
adenosquamous carcinoma (EASC), but the corresponding individual
patient could not be resolved from the publicly available metadata.
Accordingly, definitive ESCC-only primary validation was considered
unavailable. Before examining CORE scores, we therefore pre-specified
a secondary histology-sensitivity framework comprising all 22
possible assignments of the EASC identity. Each scenario excluded
the corresponding patient's technically eligible tumor observations
from the z-score reference and, when applicable, from the paired
analysis. Paired treatment effects were defined as after-treatment
minus before-treatment score. Robustness to individual paired
patients was assessed by full leave-one-patient-out recomputation,
with the omitted patient removed from both the z-score reference and
the paired effect calculation.

======================================================================
3. RESULTS — OMIX005710
======================================================================

Definitive ESCC-only external validation in OMIX005710 was not
possible because the single EASC case could not be mapped to an
individual patient using the available public metadata. We therefore
evaluated the pre-specified secondary histology-sensitivity analysis
across all 22 possible EASC assignments, without selecting or
promoting any individual scenario.

The treatment-associated change in CORE2_STRONG was negative in all
22 scenarios, both for the mean and median paired change, opposite to
the frozen replication direction from the discovery cohorts. CORE4
showed the same negative mean direction in all 22 scenarios. Both
CDKN1B and GSN contributed to the negative CORE2_STRONG direction in
all 22 scenarios, whereas the BGN/TIMP1 component was negative in 18
and positive in 4 scenarios. Independent reconstruction of the
normalization, z-scoring, and score calculations reproduced the
original analysis within floating-point precision.

The result was also insensitive to any single paired patient. Across
105 full leave-one-paired-patient-out recomputations, the mean and
median CORE2_STRONG treatment-associated changes remained negative
in every run, and all remaining paired patients retained negative
changes in every run. Thus, within this pre-specified secondary
framework, OMIX005710 showed robust directional discordance rather
than replication of the treatment-associated CORE response observed
in the discovery cohorts.

======================================================================
4. DISCUSSION
======================================================================

The OMIX005710 analysis provides an important negative external
finding. Although unresolved patient-level histology prevented a
definitive ESCC-only primary validation, the pre-specified
histology-sensitivity analysis consistently yielded a
treatment-associated CORE direction opposite to that observed in
the discovery cohorts. This directional discordance was robust to
all hypothetical EASC assignments, independent reconstruction of the
scoring pipeline, and paired leave-one-patient-out analyses, arguing
against a simple sign, pairing, reference-set, or single-patient
artifact.

These findings indicate that the treatment-associated CORE response
is not uniformly replicated across independent cohorts and may be
context-dependent. However, the present analysis does not establish
the biological cause of the discordance. In particular, differences
in CAF state composition, treatment exposure, sampling context,
clinical response, or technical platform were not tested as
confirmatory explanations and would require separately designated
exploratory analyses.

======================================================================
5. LIMITATIONS
======================================================================

A key limitation of the OMIX005710 analysis was the inability to
identify the single EASC patient at the individual-patient level from
the available public metadata, which prevented construction of a
definitive ESCC-only primary validation set. We therefore relied on a
pre-specified secondary sensitivity analysis across all 22 possible
EASC assignments. In addition, only five paired patients satisfied
the frozen fibroblast technical criteria before hypothetical EASC
exclusion, and scenarios in which a paired patient was assigned as
EASC contained four evaluable pairs. This small paired sample size
limits formal inferential resolution; accordingly, interpretation
focused on direction, patient-level consistency, scenario
robustness, and leave-one-patient-out influence rather than
statistical significance.

======================================================================
6. CONCLUSION — CALIBRATION SENTENCE
======================================================================

Together, the discovery cohorts support a reproducible fibroblast
CORE framework within the original analytical context, whereas the
independent OMIX005710 sensitivity analysis showed robust directional
discordance, indicating that the treatment-associated CORE response
is not uniformly reproduced across cohorts.

======================================================================
7. SUPPLEMENTARY MATERIAL — REQUIRED CONTENT
======================================================================

Supplementary external-validation material should include:

1. OMIX005710 sample and technical eligibility table.
2. Frozen >=20-fibroblast criterion and P15-before provenance status.
3. All 22 hypothetical EASC scenario memberships.
4. Scenario-level CORE2_STRONG and CORE4 directional summaries.
5. CDKN1B / GSN / BGN-TIMP1 component summary.
6. Full 105-run paired LOPO influence summary.
7. Technical independent-recalculation audit.
8. Explicit statement that the 22 scenarios are a sensitivity
   envelope and were not meta-analyzed.

======================================================================
"""


TEXT_OUT.write_text(
    proposed_text,
    encoding="utf-8"
)


# ============================================================
# 11. Final status
# ============================================================

print("===== 10. FINAL STEP23B STATUS =====")
print()


status_rows = [
    {
        "item":
            "manuscript_source",

        "status":
            relative(
                MANUSCRIPT
            )
    },

    {
        "item":
            "manuscript_source_SHA256",

        "status":
            observed_sha
    },

    {
        "item":
            "manuscript_modified",

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
            "evidence_hierarchy_frozen",

        "status":
            "YES"
    },

    {
        "item":
            "reporting_rules_frozen",

        "status":
            "YES"
    },

    {
        "item":
            "proposed_text_generated",

        "status":
            "YES_NOT_INSERTED"
    },

    {
        "item":
            "primary_OMIX_status",

        "status":
            "UNAVAILABLE"
    },

    {
        "item":
            "secondary_OMIX_status",

        "status":
            "ROBUSTLY_OPPOSITE"
    },

    {
        "item":
            "biological_mechanism_status",

        "status":
            "NOT_ESTABLISHED"
    },

    {
        "item":
            "next_step",

        "status":
            "STEP23C_CREATE_NEW_WORKING_COPY_AND_INTEGRATE"
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


print("=" * 80)
print("STEP23B COMPLETE")
print("=" * 80)
print()

print(
    "SOURCE MANUSCRIPT:"
)

print(
    " ",
    relative(
        MANUSCRIPT
    )
)

print()

print(
    "SOURCE SHA256:"
)

print(
    " ",
    observed_sha
)

print()

print(
    "MANUSCRIPT MODIFIED           : NO"
)

print(
    "STEP22A RECOMPUTED            : NO"
)

print()

print(
    "SECTION EVIDENCE MAP          : FROZEN"
)

print(
    "REPORTING RULES               : FROZEN"
)

print(
    "PROPOSED INSERTION TEXT       : GENERATED / NOT INSERTED"
)

print()

print(
    "OMIX PRIMARY                  : UNAVAILABLE"
)

print(
    "OMIX SECONDARY                : ROBUSTLY OPPOSITE"
)

print(
    "BIOLOGICAL MECHANISM          : NOT ESTABLISHED"
)

print()

print(
    "OUTPUTS:"
)

print(
    " ",
    relative(
        SECTION_AUDIT_OUT
    )
)

print(
    " ",
    relative(
        ANCHOR_OUT
    )
)

print(
    " ",
    relative(
        EVIDENCE_MAP_OUT
    )
)

print(
    " ",
    relative(
        RULES_OUT
    )
)

print(
    " ",
    relative(
        TEXT_OUT
    )
)

print()

print(
    "NEXT:"
)

print(
    "STEP23C_CREATE_NEW_WORKING_COPY_AND_INTEGRATE"
)

print()

print(
    "DO NOT OVERWRITE THE STEP21E SOURCE MANUSCRIPT."
)

print()
