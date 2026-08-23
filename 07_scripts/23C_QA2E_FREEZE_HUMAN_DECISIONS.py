#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP23C-QA2E
#
# FREEZE HUMAN SCIENTIFIC / EDITORIAL DECISIONS
#
# READ-ONLY WITH RESPECT TO MANUSCRIPT.
#
# DOES NOT:
#   - modify QA1B DOCX
#   - recompute Step22A
#   - perform global find/replace
#   - invent bibliographic details
#
# PURPOSE:
#   Freeze the 15 manually reviewed QA2D decisions into a
#   machine-readable decision lock before QA2B manuscript editing.
#
# IMPORTANT:
#   Remaining QA2C items are NOT automatically decided here.
# ============================================================

from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import re
import sys


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

META.mkdir(
    parents=True,
    exist_ok=True
)

RESULT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Inputs
# ============================================================

DOCX = (
    RESULT
    / "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"
)

QA2C_QUEUE = (
    META
    / "STEP23C_QA2C_54_ITEM_REVIEW_QUEUE.csv"
)

QA2D_PACKET_CSV = (
    META
    / "STEP23C_QA2D_MUST_REVIEW_PACKET.csv"
)


# ============================================================
# Outputs
# ============================================================

LOCK_CSV = (
    META
    / "STEP23C_QA2E_DECISION_LOCK.csv"
)

DOSSIER_TXT = (
    RESULT
    / "STEP23C_QA2E_DECISION_DOSSIER.txt"
)

REMAINING_CSV = (
    META
    / "STEP23C_QA2E_REMAINING_REVIEW_QUEUE.csv"
)

STATUS_CSV = (
    META
    / "STEP23C_QA2E_STATUS.csv"
)

CHECKSUMS = (
    META
    / "STEP23C_QA2E_LOCK_CHECKSUMS.sha256"
)


print()
print("=" * 84)
print("STEP23C-QA2E: FREEZE HUMAN DECISION LOCK")
print("=" * 84)
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

            block = f.read(
                1024 * 1024
            )

            if not block:
                break

            h.update(
                block
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


def relative(path):

    try:

        return str(
            path.relative_to(
                PROJECT
            )
        )

    except Exception:

        return str(path)


# ============================================================
# 1. Verify frozen inputs
# ============================================================

print("===== 1. VERIFY INPUTS =====")
print()


for path in [
    DOCX,
    QA2C_QUEUE,
    QA2D_PACKET_CSV
]:

    if not path.exists():

        fail(
            "Missing required input:\n"
            + str(path)
        )


DOCX_SHA_BEFORE = sha256_file(
    DOCX
)


queue = read_csv(
    QA2C_QUEUE
)

packet = read_csv(
    QA2D_PACKET_CSV
)


if len(
    queue
) != 54:

    fail(
        f"Expected frozen QA2C total = 54; "
        f"found {len(queue)}."
    )


if len(
    packet
) != 15:

    fail(
        f"Expected frozen QA2D packet = 15; "
        f"found {len(packet)}."
    )


queue_by_id = {
    clean(
        x.get(
            "item_id"
        )
    ):
    x
    for x in queue
}


print(
    "QA1B DOCX:"
)

print(
    " ",
    relative(
        DOCX
    )
)

print(
    "SHA256:"
)

print(
    " ",
    DOCX_SHA_BEFORE
)

print()

print(
    "QA2C frozen items:",
    len(
        queue
    )
)

print(
    "QA2D reviewed items:",
    len(
        packet
    )
)

print()


# ============================================================
# 2. Expected reviewed item IDs
# ============================================================

EXPECTED_IDS = [
    "QA2C-004",
    "QA2C-006",
    "QA2C-009",
    "QA2C-010",
    "QA2C-011",
    "QA2C-029",
    "QA2C-040",
    "QA2C-042",
    "QA2C-048",
    "QA2C-049",
    "QA2C-050",
    "QA2C-051",
    "QA2C-052",
    "QA2C-053",
    "QA2C-054",
]


packet_ids = {
    clean(
        x.get(
            "item_id"
        )
    )
    for x in packet
}


if packet_ids != set(
    EXPECTED_IDS
):

    print(
        "Expected IDs:"
    )

    for x in EXPECTED_IDS:
        print(
            " ",
            x
        )

    print()

    print(
        "Observed QA2D IDs:"
    )

    for x in sorted(
        packet_ids
    ):
        print(
            " ",
            x
        )

    fail(
        "QA2D reviewed item set does not match frozen 15-item decision set."
    )


print(
    "✓ 15/15 reviewed item IDs confirmed."
)

print()


# ============================================================
# 3. Human-approved decisions
#
# These are NOT inferred by the script.
# They encode the manually reviewed decisions made after QA2D.
# ============================================================

DECISIONS = {


# ----------------------------------------------------------------
# ABSTRACT
# ----------------------------------------------------------------

"QA2C-004": {
    "decision_group":
        "ABSTRACT_SCOPE_01",

    "section":
        "Abstract",

    "decision":
        "REVISE",

    "execution_mode":
        "CONTEXTUAL_SENTENCE_REWRITE",

    "approved_text":
        (
            "We asked whether a treatment-associated fibroblast "
            "response could be reproducibly identified across two "
            "discovery ESCC single-cell RNA-seq cohorts at the "
            "sample/patient level."
        ),

    "scientific_rationale":
        (
            "Retain the reproducibility question but explicitly limit "
            "its scope to the two discovery cohorts; do not imply "
            "OMIX005710 replication."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-006": {
    "decision_group":
        "ABSTRACT_SCOPE_02",

    "section":
        "Abstract",

    "decision":
        "REVISE",

    "execution_mode":
        "CONTEXTUAL_SENTENCE_REWRITE",

    "approved_text":
        (
            "Within the two discovery cohorts, a fibroblast CORE2 "
            "response showed reproducible same-direction "
            "treatment-associated effects."
        ),

    "scientific_rationale":
        (
            "The original 'across cohorts' wording was too broad after "
            "OMIX005710 showed robust directional discordance."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-009": {
    "decision_group":
        "ABSTRACT_CONCLUSION_03",

    "section":
        "Abstract",

    "decision":
        "MERGE_REWRITE_WITH_QA2C_010",

    "execution_mode":
        "REPLACE_COMBINED_LOGICAL_SENTENCE",

    "approved_text":
        (
            "These results nominate a fibroblast CORE2 response that "
            "was reproducible across the two discovery cohorts, while "
            "the independent OMIX005710 secondary sensitivity analysis "
            "showed directional discordance; upstream regulator and "
            "extracellular signaling candidates remain "
            "hypothesis-generating."
        ),

    "scientific_rationale":
        (
            "Discovery-cohort reproducibility and independent-cohort "
            "directional non-replication must appear in the same "
            "summary conclusion."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-010": {
    "decision_group":
        "ABSTRACT_CONCLUSION_03",

    "section":
        "Abstract",

    "decision":
        "MERGE_REWRITE_WITH_QA2C_009",

    "execution_mode":
        "REPLACE_COMBINED_LOGICAL_SENTENCE",

    "approved_text":
        (
            "These results nominate a fibroblast CORE2 response that "
            "was reproducible across the two discovery cohorts, while "
            "the independent OMIX005710 secondary sensitivity analysis "
            "showed directional discordance; upstream regulator and "
            "extracellular signaling candidates remain "
            "hypothesis-generating."
        ),

    "scientific_rationale":
        (
            "This item is the continuation of QA2C-009 and must be "
            "edited as one logical sentence, not by global token "
            "replacement."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-011": {
    "decision_group":
        "ABSTRACT_OMIX_BOUNDARY_04",

    "section":
        "Abstract",

    "decision":
        "REVISE",

    "execution_mode":
        "REPLACE_FULL_PARAGRAPH",

    "approved_text":
        (
            "In an independent neoadjuvant single-cell cohort "
            "(OMIX005710), definitive primary ESCC-only validation "
            "was precluded because the single patient with esophageal "
            "adenosquamous carcinoma could not be identified from the "
            "available patient-level metadata. In a pre-specified "
            "secondary histology-sensitivity analysis spanning all 22 "
            "possible EASC assignments, the treatment-associated CORE "
            "direction was consistently opposite to that observed in "
            "the discovery cohorts."
        ),

    "scientific_rationale":
        (
            "Make the primary-vs-secondary distinction explicit while "
            "preserving the frozen Step22A conclusion."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# RESULTS
# ----------------------------------------------------------------

"QA2C-029": {
    "decision_group":
        "RESULTS_OMIX_BOUNDARY_01",

    "section":
        "Results",

    "decision":
        "REVISE",

    "execution_mode":
        "REPLACE_FULL_PARAGRAPH",

    "approved_text":
        (
            "Definitive primary ESCC-only external validation in "
            "OMIX005710 was unavailable because the single EASC case "
            "could not be mapped to an individual patient using the "
            "available public metadata. We therefore evaluated the "
            "pre-specified secondary histology-sensitivity analysis "
            "across all 22 possible EASC assignments, without "
            "selecting or promoting any individual scenario."
        ),

    "scientific_rationale":
        (
            "Use the frozen governance term 'unavailable' and state "
            "that only the secondary sensitivity analysis was "
            "evaluable."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# DISCUSSION
# ----------------------------------------------------------------

"QA2C-040": {
    "decision_group":
        "DISCUSSION_OMIX_01",

    "section":
        "Discussion",

    "decision":
        "REVISE",

    "execution_mode":
        "REPLACE_FULL_PARAGRAPH",

    "approved_text":
        (
            "The OMIX005710 analysis provides an important secondary "
            "external finding. Although unresolved patient-level "
            "histology precluded definitive primary ESCC-only "
            "validation, the pre-specified secondary "
            "histology-sensitivity analysis consistently yielded a "
            "treatment-associated CORE direction opposite to that "
            "observed in the two discovery cohorts. This directional "
            "discordance was robust to all hypothetical EASC "
            "assignments, independent reconstruction of the scoring "
            "pipeline, and paired leave-one-patient-out analyses, "
            "arguing against a simple sign, pairing, reference-set, "
            "or single-patient artifact."
        ),

    "scientific_rationale":
        (
            "Avoid calling the unavailable primary analysis a "
            "'negative validation'; identify this as a secondary "
            "external finding and explicitly scope the comparator to "
            "the two discovery cohorts."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-042": {
    "decision_group":
        "DISCUSSION_LIMITATION_02",

    "section":
        "Discussion",

    "decision":
        "REVISE_MINOR",

    "execution_mode":
        "REPLACE_FULL_PARAGRAPH",

    "approved_text":
        (
            "A key limitation of the OMIX005710 analysis was the "
            "inability to identify the single EASC patient at the "
            "individual-patient level from the available public "
            "metadata, which precluded construction of a definitive "
            "primary ESCC-only validation set. We therefore relied on "
            "a pre-specified secondary sensitivity analysis across all "
            "22 possible EASC assignments. In addition, only five "
            "paired patients satisfied the frozen fibroblast technical "
            "criteria before hypothetical EASC exclusion, and "
            "scenarios in which a paired patient was assigned as EASC "
            "contained four evaluable pairs. This small paired sample "
            "size limits formal inferential resolution; accordingly, "
            "interpretation focused on direction, patient-level "
            "consistency, scenario robustness, and "
            "leave-one-patient-out influence rather than statistical "
            "significance."
        ),

    "scientific_rationale":
        (
            "The existing paragraph was scientifically appropriate; "
            "only the primary-validation phrasing is standardized."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# METHODS
# ----------------------------------------------------------------

"QA2C-048": {
    "decision_group":
        "METHODS_HEADING_01",

    "section":
        "Methods",

    "decision":
        "REVISE",

    "execution_mode":
        "REPLACE_HEADING",

    "approved_text":
        (
            "14. OMIX005710 external-validation framework and "
            "secondary histology-sensitivity analysis"
        ),

    "scientific_rationale":
        (
            "The heading must not imply that definitive primary "
            "external validation was successfully completed."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-049": {
    "decision_group":
        "METHODS_FRAMEWORK_02",

    "section":
        "Methods",

    "decision":
        "REVISE",

    "execution_mode":
        "REPLACE_FULL_PARAGRAPH",

    "approved_text":
        (
            "A third independent neoadjuvant single-cell "
            "RNA-sequencing cohort, OMIX005710, was evaluated under a "
            "pre-specified external-validation framework. Fibroblast "
            "pseudobulk profiles were constructed from raw counts "
            "after broad cell-type annotation using a marker framework "
            "independent of the CORE genes. A pre-specified minimum of "
            "20 fibroblast cells per patient-timepoint was required "
            "for technical eligibility; one sample (P15 before "
            "treatment, 18 fibroblasts) therefore remained in raw "
            "provenance but was excluded from the technical reference "
            "and paired analysis. Raw pseudobulk counts were "
            "transformed to library-size CPM and log2(CPM+1), followed "
            "by gene-wise z-scoring within the eligible OMIX tumor "
            "reference. CORE2_STRONG and CORE4 were calculated using "
            "the previously frozen formulas without refitting gene "
            "weights."
        ),

    "scientific_rationale":
        (
            "Use 'external-validation framework' rather than wording "
            "that implies completion of definitive primary validation."
        ),

    "qa2b_ready":
        "YES"
},


"QA2C-050": {
    "decision_group":
        "METHODS_SENSITIVITY_03",

    "section":
        "Methods",

    "decision":
        "REVISE_MINOR",

    "execution_mode":
        "REPLACE_FULL_PARAGRAPH",

    "approved_text":
        (
            "The source cohort contained one patient with esophageal "
            "adenosquamous carcinoma (EASC), but the corresponding "
            "individual patient could not be resolved from the "
            "publicly available metadata. Accordingly, definitive "
            "primary ESCC-only validation was classified as "
            "unavailable. Before reading any CORE scores, we therefore "
            "pre-specified a secondary histology-sensitivity framework "
            "comprising all 22 possible assignments of the EASC "
            "identity. Each scenario excluded the corresponding "
            "patient's technically eligible tumor observations from "
            "the z-score reference and, when applicable, from the "
            "paired analysis. Paired treatment effects were defined as "
            "after-treatment minus before-treatment score. Robustness "
            "to individual paired patients was assessed by full "
            "leave-one-patient-out recomputation, with the omitted "
            "patient removed from both the z-score reference and the "
            "paired effect calculation."
        ),

    "scientific_rationale":
        (
            "Make the frozen pre-score governance sequence and "
            "primary-unavailable classification explicit."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# SUPPLEMENTARY
# ----------------------------------------------------------------

"QA2C-051": {
    "decision_group":
        "SUPPLEMENTARY_LABEL_01",

    "section":
        "Supplementary",

    "decision":
        "REVISE",

    "execution_mode":
        "REPLACE_INDEX_ENTRY",

    "approved_text":
        (
            "Supplementary Figure S2: Discordant signaling candidates "
            "across the two discovery cohorts"
        ),

    "scientific_rationale":
        (
            "Avoid ambiguity between discovery-cohort signaling "
            "discordance and OMIX005710 external directional "
            "non-replication."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# DATA AVAILABILITY
# ----------------------------------------------------------------

"QA2C-052": {
    "decision_group":
        "DATA_AVAILABILITY_01",

    "section":
        "Data Availability",

    "decision":
        "ADD_OR_REPLACE",

    "execution_mode":
        "REPLACE_DATA_AVAILABILITY_BODY",

    "approved_text":
        (
            "GSE197677 and GSE221561 are publicly available "
            "single-cell RNA-seq datasets. OMIX005710, used for the "
            "secondary external histology-sensitivity analysis, is "
            "publicly available under accession OMIX005710. No new "
            "sequencing data were generated in this study. Accession "
            "and repository details will be provided in the final "
            "submission."
        ),

    "scientific_rationale":
        (
            "The manuscript must disclose the third public dataset "
            "used in the Step22A secondary external analysis."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# SUPPLEMENTARY MATERIAL INDEX ADDITIONS
# ----------------------------------------------------------------

"QA2C-053": {
    "decision_group":
        "SUPPLEMENTARY_STEP22A_INDEX_02",

    "section":
        "Supplementary",

    "decision":
        "ADD",

    "execution_mode":
        "APPEND_SUPPLEMENTARY_INDEX_ENTRIES",

    "approved_text":
        (
            "Supplementary Table S13: OMIX005710 sample metadata, "
            "technical eligibility, and patient/timepoint composition.\n"
            "Supplementary Table S14: Pre-specified 22-scenario "
            "histology-sensitivity results for CORE2_STRONG, CORE4, "
            "and component directions.\n"
            "Supplementary Table S15: Paired leave-one-patient-out "
            "robustness analysis and independent score-recalculation "
            "audit.\n"
            "Supplementary Note S1: OMIX005710 histology-resolution "
            "limitation, primary-validation gate, and analysis "
            "governance."
        ),

    "scientific_rationale":
        (
            "Expose the already frozen Step22A external-validation "
            "audit trail without creating new biological analyses."
        ),

    "qa2b_ready":
        "YES"
},


# ----------------------------------------------------------------
# REFERENCES
# ----------------------------------------------------------------

"QA2C-054": {
    "decision_group":
        "REFERENCES_01",

    "section":
        "References",

    "decision":
        "REQUIRE_VERIFIED_CITATIONS",

    "execution_mode":
        "HOLD_FOR_BIBLIOGRAPHIC_VERIFICATION",

    "approved_text":
        "",

    "scientific_rationale":
        (
            "References must include the real source publications for "
            "GSE197677, GSE221561, and OMIX005710. Author names, title, "
            "journal metadata, and final reference formatting must be "
            "verified from primary publication records before insertion; "
            "the decision lock explicitly forbids guessing bibliographic "
            "details."
        ),

    "qa2b_ready":
        "NO_REFERENCE_VERIFICATION_REQUIRED"
},

}


# ============================================================
# 4. Verify decision completeness
# ============================================================

print("===== 2. VERIFY HUMAN DECISION SET =====")
print()


if set(
    DECISIONS.keys()
) != set(
    EXPECTED_IDS
):

    fail(
        "Decision dictionary does not contain exactly the frozen "
        "15 QA2D items."
    )


for item_id in EXPECTED_IDS:

    if item_id not in queue_by_id:

        fail(
            "Frozen QA2C item missing from queue: "
            + item_id
        )


print(
    "✓ Human decisions defined:",
    len(
        DECISIONS
    ),
    "/ 15"
)

print()


# ============================================================
# 5. Build lock rows
# ============================================================

lock_rows = []


for item_id in EXPECTED_IDS:

    source = queue_by_id[
        item_id
    ]

    decision = DECISIONS[
        item_id
    ]


    lock_rows.append({

        "item_id":
            item_id,

        "decision_group":
            decision[
                "decision_group"
            ],

        "section":
            decision[
                "section"
            ],

        "severity":
            clean(
                source.get(
                    "severity"
                )
            ),

        "risk_family":
            clean(
                source.get(
                    "risk_family"
                )
            ),

        "original_text":
            clean(
                source.get(
                    "text"
                )
            ),

        "human_decision":
            decision[
                "decision"
            ],

        "execution_mode":
            decision[
                "execution_mode"
            ],

        "approved_text":
            decision[
                "approved_text"
            ],

        "scientific_rationale":
            decision[
                "scientific_rationale"
            ],

        "qa2b_ready":
            decision[
                "qa2b_ready"
            ],

        "lock_status":
            "LOCKED"
    })


write_csv(
    LOCK_CSV,
    lock_rows,
    [
        "item_id",
        "decision_group",
        "section",
        "severity",
        "risk_family",
        "original_text",
        "human_decision",
        "execution_mode",
        "approved_text",
        "scientific_rationale",
        "qa2b_ready",
        "lock_status"
    ]
)


# ============================================================
# 6. Preserve all remaining QA2C items as explicitly UNDECIDED
# ============================================================

remaining = [
    x
    for x in queue
    if clean(
        x.get(
            "item_id"
        )
    )
    not in set(
        EXPECTED_IDS
    )
]


if len(
    remaining
) != 39:

    fail(
        f"Expected 39 remaining QA2C items; "
        f"found {len(remaining)}."
    )


for row in remaining:

    row[
        "qa2e_status"
    ] = (
        "PENDING_HUMAN_REVIEW_NOT_AUTHORIZED_FOR_QA2B"
    )


remaining.sort(
    key=lambda x:
        int(
            clean(
                x.get(
                    "review_order"
                )
            )
            or 9999
        )
)


remaining_fields = list(
    remaining[0].keys()
)


write_csv(
    REMAINING_CSV,
    remaining,
    remaining_fields
)


remaining_buckets = Counter(
    clean(
        x.get(
            "review_bucket"
        )
    )
    for x in remaining
)


remaining_sections = Counter(
    clean(
        x.get(
            "section"
        )
    )
    for x in remaining
)


# ============================================================
# 7. Human-readable dossier
# ============================================================

lines = [
    "STEP23C-QA2E HUMAN DECISION LOCK DOSSIER",
    "=" * 84,
    "",
    "STATUS:",
    "HUMAN_DECISION_LOCK_FROZEN",
    "",
    "QA1B MANUSCRIPT:",
    relative(
        DOCX
    ),
    "",
    "QA1B SHA256:",
    DOCX_SHA_BEFORE,
    "",
    "LOCKED QA2D ITEMS:",
    "15 / 15",
    "",
    "REMAINING QA2C ITEMS:",
    "39 — explicitly NOT authorized for QA2B editing yet.",
    "",
    "GLOBAL EDITING RULE:",
    "NO global find/replace for conserved, reproducible, validation,",
    "mechanism, replicated, uniform, or related terminology.",
    "",
    "SCIENTIFIC BOUNDARIES LOCKED:",
    "",
    "1. Discovery evidence may use reproducibility/conservation language",
    "   only when explicitly limited to the two discovery cohorts.",
    "",
    "2. OMIX005710 definitive primary ESCC-only validation is UNAVAILABLE,",
    "   not 'failed'.",
    "",
    "3. OMIX005710 evaluated a pre-specified SECONDARY",
    "   histology-sensitivity analysis.",
    "",
    "4. The OMIX treatment-associated CORE direction showed robust",
    "   directional discordance / non-replication relative to the",
    "   discovery-cohort direction.",
    "",
    "5. The 22 histology scenarios are a sensitivity envelope, not",
    "   independent cohorts or independent replications.",
    "",
    "6. The 105 LOPO runs are robustness/influence audits, not independent",
    "   significance tests.",
    "",
    "7. No reverse biological mechanism has been established.",
    "",
    "8. Reference metadata must be verified from primary publication",
    "   records before bibliography insertion.",
    "",
    "=" * 84,
    "LOCKED DECISIONS",
    "=" * 84,
    ""
]


for n, row in enumerate(
    lock_rows,
    start=1
):

    lines.extend([
        f"{n:02d}. {row['item_id']} — {row['section']}",
        f"Decision group : {row['decision_group']}",
        f"Decision       : {row['human_decision']}",
        f"Execution mode : {row['execution_mode']}",
        f"QA2B ready     : {row['qa2b_ready']}",
        "",
        "ORIGINAL:",
        row[
            "original_text"
        ]
        or "<packaging item without paragraph text>",
        "",
        "APPROVED TEXT:"
    ])


    if row[
        "approved_text"
    ]:

        lines.append(
            row[
                "approved_text"
            ]
        )

    else:

        lines.append(
            "<NO TEXT LOCKED — BIBLIOGRAPHIC VERIFICATION REQUIRED>"
        )


    lines.extend([
        "",
        "RATIONALE:",
        row[
            "scientific_rationale"
        ],
        "",
        "-" * 84,
        ""
    ])


lines.extend([
    "=" * 84,
    "REMAINING REVIEW QUEUE",
    "=" * 84,
    "",
    f"Total remaining: {len(remaining)}",
    "",
    "By review bucket:"
])


for key, value in remaining_buckets.most_common():

    lines.append(
        f"  {key}: {value}"
    )


lines.extend([
    "",
    "By manuscript section:"
])


for key, value in remaining_sections.most_common():

    lines.append(
        f"  {key or '<unspecified>'}: {value}"
    )


lines.extend([
    "",
    "IMPORTANT:",
    "These 39 items remain PENDING and must not be changed by QA2B until",
    "they are separately reviewed and locked.",
    "",
    "TITLE / RUNNING-TITLE NOTE:",
    "Any remaining Title/Running-title conserved/reproducible language",
    "remains pending discovery-scope review unless explicitly represented",
    "in the 15-item lock above.",
    "",
    "REFERENCES NOTE:",
    "QA2C-054 locks the REQUIREMENT to complete real citations, but does",
    "not authorize invented or unverified bibliography strings.",
    "",
    "FINAL STATUS:",
    "STEP23C_QA2E_DECISION_LOCK_COMPLETE",
    ""
])


DOSSIER_TXT.write_text(
    "\n".join(
        lines
    ),
    encoding="utf-8"
)


# ============================================================
# 8. Confirm manuscript untouched
# ============================================================

DOCX_SHA_AFTER = sha256_file(
    DOCX
)


if (
    DOCX_SHA_AFTER
    != DOCX_SHA_BEFORE
):

    fail(
        "CRITICAL: QA1B manuscript changed during QA2E."
    )


# ============================================================
# 9. Lock checksums
# ============================================================

lock_sha = sha256_file(
    LOCK_CSV
)

dossier_sha = sha256_file(
    DOSSIER_TXT
)

remaining_sha = sha256_file(
    REMAINING_CSV
)


CHECKSUMS.write_text(
    "\n".join([
        f"{lock_sha}  {LOCK_CSV.name}",
        f"{dossier_sha}  {DOSSIER_TXT.name}",
        f"{remaining_sha}  {REMAINING_CSV.name}",
        ""
    ]),
    encoding="utf-8"
)


# ============================================================
# 10. Status file
# ============================================================

qa2b_ready_count = sum(
    1
    for row in lock_rows
    if row[
        "qa2b_ready"
    ]
    == "YES"
)


reference_hold_count = sum(
    1
    for row in lock_rows
    if row[
        "qa2b_ready"
    ]
    != "YES"
)


status_rows = [
    {
        "item":
            "QA1B_DOCX",

        "status":
            relative(
                DOCX
            )
    },

    {
        "item":
            "QA1B_SHA256_before",

        "status":
            DOCX_SHA_BEFORE
    },

    {
        "item":
            "QA1B_SHA256_after",

        "status":
            DOCX_SHA_AFTER
    },

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
            "QA2D_decisions_locked",

        "status":
            "15"
    },

    {
        "item":
            "QA2B_ready_locked_items",

        "status":
            str(
                qa2b_ready_count
            )
    },

    {
        "item":
            "reference_verification_hold",

        "status":
            str(
                reference_hold_count
            )
    },

    {
        "item":
            "remaining_QA2C_items",

        "status":
            str(
                len(
                    remaining
                )
            )
    },

    {
        "item":
            "global_find_replace_authorized",

        "status":
            "NO"
    },

    {
        "item":
            "automatic_manuscript_rewrite_performed",

        "status":
            "NO"
    },

    {
        "item":
            "decision_lock_SHA256",

        "status":
            lock_sha
    },

    {
        "item":
            "status",

        "status":
            "STEP23C_QA2E_DECISION_LOCK_COMPLETE"
    },

    {
        "item":
            "next_step",

        "status":
            "REVIEW_REMAINING_DISCOVERY_SCOPE_ITEMS_BEFORE_QA2B"
    }
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
print("STEP23C-QA2E COMPLETE")
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

print(
    "GLOBAL FIND/REPLACE           : NOT AUTHORIZED"
)

print()

print(
    "HUMAN DECISIONS LOCKED        :",
    len(
        lock_rows
    ),
    "/ 15"
)

print(
    "QA2B-READY LOCKED ITEMS       :",
    qa2b_ready_count
)

print(
    "REFERENCE VERIFICATION HOLD   :",
    reference_hold_count
)

print()

print(
    "REMAINING QA2C ITEMS          :",
    len(
        remaining
    )
)

print()

print(
    "REMAINING REVIEW BUCKETS:"
)

for key, value in remaining_buckets.most_common():

    print(
        f"  {key:<42s}: {value}"
    )


print()

print(
    "OUTPUTS:"
)

print(
    " ",
    relative(
        LOCK_CSV
    )
)

print(
    " ",
    relative(
        DOSSIER_TXT
    )
)

print(
    " ",
    relative(
        REMAINING_CSV
    )
)

print(
    " ",
    relative(
        CHECKSUMS
    )
)

print()

print(
    "STATUS:"
)

print(
    "STEP23C_QA2E_DECISION_LOCK_COMPLETE"
)

print()

print(
    "NEXT:"
)

print(
    "REVIEW REMAINING DISCOVERY-SCOPE / CONTEXT ITEMS BEFORE QA2B."
)

print()

