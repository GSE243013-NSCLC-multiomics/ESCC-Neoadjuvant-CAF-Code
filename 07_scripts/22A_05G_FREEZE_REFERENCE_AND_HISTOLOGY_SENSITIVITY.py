#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5G
#
# FREEZE:
#   1. exact OMIX z-score reference framework
#   2. histology-unresolved sensitivity design
#   3. remaining all-tumor technical target inventory
#
# IMPORTANT:
#   NO expression values are read.
#   NO pseudobulk count matrix is read.
#   NO normalization is executed.
#   NO CORE2 / CORE4 score is calculated.
#   NO Before/After inference is performed.
#   NO files are downloaded.
# ============================================================

from pathlib import Path
import csv
import hashlib
import sys


PROJECT = Path.cwd()

META_DIR = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

GATE_FILE = (
    META_DIR
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

IMPLEMENTATION_CHECK = (
    META_DIR
    / "STEP22A_05F_IMPLEMENTATION_RESOLUTION_CHECK.csv"
)

SCORING_LOCK = (
    META_DIR
    / "STEP22A_05F_EXACT_SCORING_IMPLEMENTATION_LOCK.txt"
)

CANONICAL_META = (
    META_DIR
    / "STEP22A_SAMPLE_METADATA_CANONICAL.csv"
)

DOWNLOAD_PLAN = (
    META_DIR
    / "STEP22A_TARGETED_TUMOR_DOWNLOAD_PLAN.csv"
)

CURRENT_PT_META = (
    META_DIR
    / "STEP22A_05C_FIBROBLAST_PATIENT_TIMEPOINT_METADATA.csv"
)

FROZEN_INV = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

META_DIR.mkdir(
    parents=True,
    exist_ok=True
)


print()
print("=" * 80)
print("STEP22A-5G: FREEZE Z-SCORE REFERENCE + HISTOLOGY SENSITIVITY")
print("=" * 80)
print()


# ============================================================
# Helpers
# ============================================================

def read_csv(path, delimiter=","):

    if not path.exists():

        raise SystemExit(
            "ERROR: missing required file:\n"
            + str(path)
        )

    with open(
        path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        return list(
            csv.DictReader(
                f,
                delimiter=delimiter
            )
        )


def write_csv(path, rows, fields, delimiter=","):

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter=delimiter,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def patient_num(x):

    try:
        return int(
            str(x)
            .upper()
            .replace(
                "P",
                ""
            )
        )
    except Exception:
        return 999999


def sha256_file(path):

    h = hashlib.sha256()

    with open(
        path,
        "rb"
    ) as f:

        while True:

            chunk = f.read(
                1024 * 1024
            )

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


# ============================================================
# 1. Verify gates
# ============================================================

print("===== 1. ANALYSIS GATE CHECK =====")
print()

gate_rows = read_csv(
    GATE_FILE
)

gates = {
    x["gate"]:
        x["allowed"]
    for x in gate_rows
}


required_gates = [
    "technical_preprocessing",
    "paired_CORE2_test",
    "manuscript_update"
]

for g in required_gates:

    if g not in gates:

        raise SystemExit(
            "ERROR: missing gate: "
            + g
        )


print(
    "Technical preprocessing :",
    gates[
        "technical_preprocessing"
    ]
)

print(
    "Paired CORE2 inference  :",
    gates[
        "paired_CORE2_test"
    ]
)

print(
    "Manuscript update       :",
    gates[
        "manuscript_update"
    ]
)

print()


if (
    gates[
        "technical_preprocessing"
    ]
    != "YES"
):

    raise SystemExit(
        "TECHNICAL_PREPROCESSING_BLOCKED"
    )


if (
    gates[
        "paired_CORE2_test"
    ]
    == "YES"
):

    raise SystemExit(
        "ERROR: CORE2 inference gate unexpectedly open."
    )


if (
    gates[
        "manuscript_update"
    ]
    == "YES"
):

    raise SystemExit(
        "ERROR: manuscript gate unexpectedly open."
    )


print(
    "✓ Technical preparation may continue."
)

print(
    "✓ CORE2 primary inference remains blocked."
)

print()


# ============================================================
# 2. Verify Step22A-5F implementation lock
# ============================================================

print("===== 2. VERIFY STEP22A-5F IMPLEMENTATION LOCK =====")
print()


if not SCORING_LOCK.exists():

    raise SystemExit(
        "ERROR: missing Step22A-5F scoring lock."
    )


checks = read_csv(
    IMPLEMENTATION_CHECK
)

review = [
    x
    for x in checks
    if x.get(
        "status",
        ""
    ) != "PASS"
]


print(
    "5F checks:",
    len(
        checks
    )
)

print(
    "PASS:",
    len(
        checks
    )
    -
    len(
        review
    )
)

print(
    "REVIEW:",
    len(
        review
    )
)

print()


if review:

    for x in review:

        print(
            x.get(
                "check",
                ""
            ),
            "->",
            x.get(
                "detail",
                ""
            )
        )

    raise SystemExit(
        "STEP22A_5F_NOT_FULLY_RESOLVED"
    )


lock_sha = sha256_file(
    SCORING_LOCK
)

print(
    "5F scoring-lock SHA256:",
    lock_sha
)

print()

print(
    "✓ Exact scoring implementation already locked."
)

print()


# ============================================================
# 3. Verify canonical cohort structure
#
# We use patient/tissue/timepoint fields ONLY.
# We DO NOT trust the old histology column as final,
# because EASC remains unresolved.
# ============================================================

print("===== 3. COHORT STRUCTURE CHECK =====")
print()


canonical = read_csv(
    CANONICAL_META
)


patients = sorted(
    {
        x.get(
            "patient_id",
            ""
        )
        for x in canonical
        if x.get(
            "patient_id",
            ""
        )
    },
    key=patient_num
)


print(
    "Canonical metadata rows:",
    len(
        canonical
    )
)

print(
    "Source patients:",
    len(
        patients
    )
)

print(
    "Patients:",
    ", ".join(
        patients
    )
)

print()


if len(patients) != 22:

    raise SystemExit(
        "ERROR: expected 22 source patients."
    )


tumor_rows = [
    x
    for x in canonical
    if (
        x.get(
            "canonical_tissue",
            ""
        )
        == "Tumor"
    )
]


tumor_patients = sorted(
    {
        x[
            "patient_id"
        ]
        for x in tumor_rows
    },
    key=patient_num
)


print(
    "Tumor metadata rows:",
    len(
        tumor_rows
    )
)

print(
    "Patients with tumor observations:",
    len(
        tumor_patients
    )
)

print()


# ============================================================
# 4. Verify current paired technical subset
# ============================================================

print("===== 4. CURRENT TECHNICAL SUBSET =====")
print()


pt_rows = read_csv(
    CURRENT_PT_META
)


current_paired_patients = sorted(
    {
        x.get(
            "patient_id",
            ""
        )
        for x in pt_rows
        if x.get(
            "patient_id",
            ""
        )
    },
    key=patient_num
)


print(
    "Current patient-timepoint pseudobulks:",
    len(
        pt_rows
    )
)

print(
    "Current paired patients:",
    len(
        current_paired_patients
    )
)

print(
    "Paired IDs:",
    ", ".join(
        current_paired_patients
    )
)

print()


if (
    len(
        pt_rows
    )
    != 12
):

    raise SystemExit(
        "ERROR: expected 12 current paired patient-timepoints."
    )


if (
    len(
        current_paired_patients
    )
    != 6
):

    raise SystemExit(
        "ERROR: expected 6 paired patients."
    )


for p in current_paired_patients:

    times = {
        x.get(
            "timepoint",
            ""
        )
        for x in pt_rows
        if x.get(
            "patient_id",
            ""
        )
        == p
    }

    if times != {
        "Before",
        "After"
    }:

        raise SystemExit(
            "ERROR: incomplete pair for "
            + p
        )


print(
    "✓ Current 12 pseudobulks are a valid paired technical subset."
)

print()


# ============================================================
# 5. Lock exact scoring-reference semantics
# ============================================================

print("===== 5. FREEZE Z-SCORE REFERENCE FRAMEWORK =====")
print()


reference_rules = [

    {
        "rule":
            "expression_transform",

        "value":
            "log2((raw_count/library_size)*1e6 + 1)"
    },

    {
        "rule":
            "library_size_definition",

        "value":
            "column_sum_raw_pseudobulk_counts"
    },

    {
        "rule":
            "zero_library_guard",

        "value":
            "library_size_zero_replaced_by_1_as_frozen_Step19L"
    },

    {
        "rule":
            "zscore_axis",

        "value":
            "gene_wise_across_tumor_patient_timepoints"
    },

    {
        "rule":
            "zscore_mean",

        "value":
            "mean_within_current_eligible_OMIX_tumor_reference"
    },

    {
        "rule":
            "zscore_sd",

        "value":
            "R_sample_sd_within_current_eligible_OMIX_tumor_reference"
    },

    {
        "rule":
            "zscore_sd_floor",

        "value":
            "1e-10_as_frozen_Step19L"
    },

    {
        "rule":
            "zscore_reference_dataset",

        "value":
            "OMIX005710_INTERNAL"
    },

    {
        "rule":
            "zscore_reference_universe",

        "value":
            "ALL_ELIGIBLE_ESCC_TUMOR_PATIENT_TIMEPOINTS"
    },

    {
        "rule":
            "paired_only_reference",

        "value":
            "NO"
    },

    {
        "rule":
            "reuse_GSE197677_mean_sd",

        "value":
            "NO"
    },

    {
        "rule":
            "reuse_GSE221561_mean_sd",

        "value":
            "NO"
    },

    {
        "rule":
            "CORE4_formula",

        "value":
            "(BGN_z-CDKN1B_z-GSN_z+TIMP1_z)/4"
    },

    {
        "rule":
            "CORE2_STRONG_formula",

        "value":
            "(-CDKN1B_z-GSN_z)/2"
    },

    {
        "rule":
            "CORE2_WEAK_formula",

        "value":
            "(BGN_z+TIMP1_z)/2"
    },

    {
        "rule":
            "weight_refitting",

        "value":
            "FORBIDDEN"
    },

    {
        "rule":
            "cutoff_optimization",

        "value":
            "FORBIDDEN"
    },

    {
        "rule":
            "normalization_selection_by_external_result",

        "value":
            "FORBIDDEN"
    },

    {
        "rule":
            "histology_requirement_primary",

        "value":
            "ESCC_ONLY"
    },

    {
        "rule":
            "current_histology_status",

        "value":
            "INCOMPLETE"
    },

    {
        "rule":
            "primary_score_calculation_now",

        "value":
            "BLOCKED"
    },

    {
        "rule":
            "primary_inference_now",

        "value":
            "BLOCKED"
    }
]


REFERENCE_OUT = (
    META_DIR
    / "STEP22A_05G_ZSCORE_REFERENCE_LOCK.csv"
)


write_csv(
    REFERENCE_OUT,
    reference_rules,
    [
        "rule",
        "value"
    ]
)


for x in reference_rules:

    if x["rule"] in {
        "zscore_reference_dataset",
        "zscore_reference_universe",
        "paired_only_reference",
        "reuse_GSE197677_mean_sd",
        "reuse_GSE221561_mean_sd",
        "primary_score_calculation_now"
    }:

        print(
            f"{x['rule']:<34s}: "
            f"{x['value']}"
        )


print()

print(
    "✓ OMIX internal z-score reference locked."
)

print(
    "✓ Prior-cohort mean/SD reuse explicitly forbidden."
)

print(
    "✓ Paired-only 12-sample reference explicitly forbidden."
)

print()


# ============================================================
# 6. Why current 12 samples are NOT the final reference
# ============================================================

print("===== 6. FINAL REFERENCE COMPLETENESS CHECK =====")
print()


download_plan = read_csv(
    DOWNLOAD_PLAN
)


planned_tumor = [
    x
    for x in download_plan
    if x.get(
        "patient_id",
        ""
    )
]


planned_sample_names = {
    x.get(
        "sample_name",
        ""
    )
    for x in planned_tumor
}


current_sample_names = set()


for row in pt_rows:

    source_samples = (
        row.get(
            "source_samples",
            ""
        )
        or ""
    )

    for x in source_samples.split(
        ";"
    ):

        x = x.strip()

        if x:
            current_sample_names.add(
                x
            )


remaining_sample_names = sorted(
    planned_sample_names
    -
    current_sample_names
)


print(
    "Tumor files in frozen download plan:",
    len(
        planned_sample_names
    )
)

print(
    "Tumor samples technically processed now:",
    len(
        current_sample_names
    )
)

print(
    "Remaining tumor samples:",
    len(
        remaining_sample_names
    )
)

print()


if (
    len(
        planned_sample_names
    )
    != 27
):

    raise SystemExit(
        "ERROR: expected frozen tumor download plan to contain 27 samples."
    )


if (
    len(
        current_sample_names
    )
    != 12
):

    raise SystemExit(
        "ERROR: expected 12 currently processed tumor samples."
    )


if (
    len(
        remaining_sample_names
    )
    != 15
):

    raise SystemExit(
        "ERROR: expected 15 remaining tumor samples."
    )


print(
    "FINAL OMIX Z-SCORE REFERENCE COMPLETE : NO"
)

print(
    "Reason:"
)

print(
    "  paired-only technical subset is 12/27 planned tumor samples"
)

print(
    "  and EASC identity remains unresolved"
)

print()


# ============================================================
# 7. Generate remaining technical target list
#
# NO DOWNLOAD in this step.
# ============================================================

remaining_rows = [
    x
    for x in planned_tumor
    if x.get(
        "sample_name",
        ""
    )
    in remaining_sample_names
]


REMAINING_OUT = (
    META_DIR
    / "STEP22A_05G_REMAINING_15_TUMOR_TECHNICAL_TARGETS.tsv"
)


remaining_fields = [
    "priority",
    "patient_id",
    "timepoint",
    "sample_name",
    "hrs_accession",
    "omix_file_id",
    "listed_size_MB",
    "download_url_https",
    "download_now",
    "reason"
]


write_csv(
    REMAINING_OUT,
    remaining_rows,
    remaining_fields,
    delimiter="\t"
)


remaining_mb = 0.0


for x in remaining_rows:

    try:

        remaining_mb += float(
            x.get(
                "listed_size_MB",
                "0"
            )
            or 0
        )

    except Exception:
        pass


print(
    "Remaining technical tumor files:",
    len(
        remaining_rows
    )
)

print(
    "Listed remaining size:",
    round(
        remaining_mb,
        2
    ),
    "MB"
)

print(
    "Listed remaining size:",
    round(
        remaining_mb / 1024,
        3
    ),
    "GiB"
)

print()

print(
    "Downloaded now: 0"
)

print()


# ============================================================
# 8. Freeze histology sensitivity scenarios
#
# The unknown EASC could be any source-study patient.
#
# We pre-specify ALL possible patient identities now,
# before reading any CORE expression or scores.
#
# If a hypothetical EASC patient has tumor observations:
#   exclude all that patient's tumor observations and
#   recompute OMIX-internal mean/SD in that hypothetical
#   eligible tumor reference.
#
# If that patient belongs to the paired primary set:
#   that patient is also removed from the paired effect.
#
# These are SENSITIVITY scenarios only.
#
# They DO NOT bypass the frozen ESCC-only primary gate.
# ============================================================

print("===== 7. FREEZE UNKNOWN-EASC SENSITIVITY DESIGN =====")
print()


tumor_patient_set = set(
    tumor_patients
)

paired_set = set(
    current_paired_patients
)


sensitivity_rows = []


# Scenario 0 = unresolved, no exclusion.
sensitivity_rows.append({

    "scenario_id":
        "S00_UNRESOLVED_NO_EXCLUSION",

    "hypothetical_EASC_patient":
        "UNKNOWN",

    "patient_has_tumor_observation":
        "NA",

    "patient_is_current_paired_primary":
        "NA",

    "reference_action":
        "NO_EXCLUSION_UNRESOLVED_REFERENCE",

    "zscore_action":
        "DO_NOT_CALCULATE_WHILE_GATE_BLOCKED",

    "paired_effect_action":
        "DO_NOT_CALCULATE_WHILE_GATE_BLOCKED",

    "purpose":
        "STATUS_REFERENCE_ONLY",

    "may_be_called_primary_ESCC_validation":
        "NO"
})


for p in patients:

    has_tumor = (
        p in tumor_patient_set
    )

    is_paired = (
        p in paired_set
    )

    sensitivity_rows.append({

        "scenario_id":
            "S_HYPOTHETICAL_EASC_"
            + p,

        "hypothetical_EASC_patient":
            p,

        "patient_has_tumor_observation":
            (
                "YES"
                if has_tumor
                else "NO"
            ),

        "patient_is_current_paired_primary":
            (
                "YES"
                if is_paired
                else "NO"
            ),

        "reference_action":
            (
                "EXCLUDE_ALL_TUMOR_OBSERVATIONS_FOR_"
                + p
                + "_THEN_REBUILD_ELIGIBLE_OMIX_REFERENCE"
                if has_tumor
                else
                "NO_TUMOR_REFERENCE_CHANGE"
            ),

        "zscore_action":
            (
                "RECOMPUTE_LOG2CPM_MEAN_SD_AND_Z_WITHIN_REMAINING_OMIX_TUMORS"
                if has_tumor
                else
                "UNCHANGED_REFERENCE"
            ),

        "paired_effect_action":
            (
                "REMOVE_"
                + p
                + "_FROM_PAIRED_EFFECT"
                if is_paired
                else
                "PAIRED_PATIENT_SET_UNCHANGED"
            ),

        "purpose":
            "PRE_SPECIFIED_HISTOLOGY_SENSITIVITY_ONLY",

        "may_be_called_primary_ESCC_validation":
            "NO_UNLESS_TRUE_HISTOLOGY_IDENTITY_RESOLVED"
    })


SENSITIVITY_OUT = (
    META_DIR
    / "STEP22A_05G_HISTOLOGY_SENSITIVITY_SCENARIOS.csv"
)


write_csv(
    SENSITIVITY_OUT,
    sensitivity_rows,
    [
        "scenario_id",
        "hypothetical_EASC_patient",
        "patient_has_tumor_observation",
        "patient_is_current_paired_primary",
        "reference_action",
        "zscore_action",
        "paired_effect_action",
        "purpose",
        "may_be_called_primary_ESCC_validation"
    ]
)


print(
    "Sensitivity scenarios frozen:",
    len(
        sensitivity_rows
    )
)

print(
    "  unresolved baseline:",
    1
)

print(
    "  hypothetical EASC identities:",
    len(
        patients
    )
)

print()


print(
    "Current paired patients potentially affected if EASC:"
)

for p in current_paired_patients:

    print(
        " ",
        p
    )


print()

print(
    "IMPORTANT:"
)

print(
    "Sensitivity scenarios DO NOT reopen primary inference."
)

print(
    "They are frozen now only to prevent post-result redesign."
)

print()


# ============================================================
# 9. Separate histology sensitivity from statistical LOPO
# ============================================================

distinction_rows = [

    {
        "analysis":
            "HISTOLOGY_EASC_SENSITIVITY",

        "purpose":
            "unknown_histology_identity",

        "changes_eligible_reference_universe":
            "YES_IF_PATIENT_HAS_TUMOR",

        "recompute_zscore_reference":
            "YES",

        "paired_patient_may_be_removed":
            "YES_IF_HYPOTHETICAL_EASC_IS_PAIRED",

        "allowed_now":
            "DESIGN_ONLY_NO_SCORE"
    },

    {
        "analysis":
            "FROZEN_PRIMARY_LOPO",

        "purpose":
            "patient_level_effect_stability",

        "changes_eligible_reference_universe":
            "NO_AFTER_PRIMARY_SCORES_ARE_FROZEN",

        "recompute_zscore_reference":
            "NO_IN_THIS_5G_DESIGN",

        "paired_patient_may_be_removed":
            "YES_ONE_AT_A_TIME",

        "allowed_now":
            "NO_PRIMARY_SCORE_NOT_AVAILABLE"
    }
]


DISTINCTION_OUT = (
    META_DIR
    / "STEP22A_05G_HISTOLOGY_VS_LOPO_RULES.csv"
)


write_csv(
    DISTINCTION_OUT,
    distinction_rows,
    [
        "analysis",
        "purpose",
        "changes_eligible_reference_universe",
        "recompute_zscore_reference",
        "paired_patient_may_be_removed",
        "allowed_now"
    ]
)


# ============================================================
# 10. Human-readable decision lock
# ============================================================

LOCK_OUT = (
    META_DIR
    / "STEP22A_05G_REFERENCE_AND_HISTOLOGY_DECISION_LOCK.txt"
)


with open(
    LOCK_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-5G REFERENCE + HISTOLOGY DECISION LOCK\n"
    )

    f.write(
        "=" * 70
        + "\n\n"
    )

    f.write(
        "SCORING REFERENCE DECISION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Use OMIX005710-INTERNAL normalization and z-scoring.\n"
    )

    f.write(
        "Do NOT reuse GSE197677 mean/SD.\n"
    )

    f.write(
        "Do NOT reuse GSE221561 mean/SD.\n"
    )

    f.write(
        "Do NOT use paired-only 12 samples as the final z-score reference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Final frozen reference universe:\n"
    )

    f.write(
        "ALL ELIGIBLE ESCC TUMOR PATIENT-TIMEPOINT PSEUDOBULKS\n"
    )

    f.write(
        "within OMIX005710 after frozen fibroblast QC.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Frozen transformation:\n"
    )

    f.write(
        "raw pseudobulk counts -> CPM by library size -> log2(CPM+1)\n"
    )

    f.write(
        "-> per-gene mean/SD across eligible OMIX tumor reference\n"
    )

    f.write(
        "-> z-score -> fixed CORE formulas.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CURRENT REFERENCE STATUS\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Current paired technical samples: 12\n"
    )

    f.write(
        "Frozen tumor download plan: 27\n"
    )

    f.write(
        "Remaining technical tumor samples: 15\n"
    )

    f.write(
        "Final z-score reference complete: NO\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "HISTOLOGY STATUS\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Source study reports one EASC patient.\n"
    )

    f.write(
        "EASC identity: UNRESOLVED\n"
    )

    f.write(
        "Primary ESCC-only score calculation: BLOCKED\n"
    )

    f.write(
        "Primary CORE2 inference: BLOCKED\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "SENSITIVITY DESIGN\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "All 22 possible source-patient EASC identities are frozen\n"
    )

    f.write(
        "as hypothetical sensitivity scenarios BEFORE any external\n"
    )

    f.write(
        "CORE expression or score is read.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "If the hypothetical EASC has tumor observations:\n"
    )

    f.write(
        "exclude that patient's tumor observations and recompute\n"
    )

    f.write(
        "OMIX-internal log2(CPM+1) mean/SD/z-score reference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "If the hypothetical EASC is one of the paired patients:\n"
    )

    f.write(
        "also remove that patient from the paired effect analysis.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "These sensitivity scenarios NEVER substitute for the\n"
    )

    f.write(
        "frozen ESCC-only primary histology requirement.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "NO SCORE CALCULATED IN STEP22A-5G.\n"
    )


# ============================================================
# 11. Frozen project integrity
# ============================================================

print("===== 8. FROZEN PROJECT CHECK =====")
print()


frozen = read_csv(
    FROZEN_INV
)

changes = []


for row in frozen:

    path = Path(
        row.get(
            "path",
            ""
        )
    )

    if not path.is_absolute():

        path = PROJECT / path


    if not path.exists():

        changes.append(
            (
                str(path),
                "MISSING"
            )
        )

        continue


    expected = row.get(
        "size_bytes",
        ""
    )


    if expected:

        try:

            expected = int(
                float(
                    expected
                )
            )

            actual = (
                path.stat().st_size
            )

            if actual != expected:

                changes.append(
                    (
                        str(path),
                        "SIZE_CHANGED"
                    )
                )

        except Exception:
            pass


print(
    "Frozen original files changed:",
    len(
        changes
    )
)


if changes:

    for x in changes[:20]:
        print(x)

    raise SystemExit(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )


print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# FINAL
# ============================================================

print("=" * 80)
print("STEP22A-5G COMPLETE")
print("=" * 80)
print()

print(
    "SCORING IMPLEMENTATION       : LOCKED"
)

print(
    "EXPRESSION TRANSFORM         : log2(CPM+1)"
)

print(
    "ZSCORE REFERENCE             : OMIX005710 INTERNAL"
)

print(
    "ZSCORE UNIVERSE              : ALL ELIGIBLE ESCC TUMOR PTs"
)

print(
    "REUSE PRIOR COHORT MEAN/SD   : NO"
)

print(
    "PAIRED-ONLY ZSCORE REFERENCE : NO"
)

print()

print(
    "CURRENT TECHNICAL TUMORS     :",
    len(
        current_sample_names
    )
)

print(
    "FROZEN TUMOR PLAN            :",
    len(
        planned_sample_names
    )
)

print(
    "REMAINING TECHNICAL TUMORS   :",
    len(
        remaining_sample_names
    )
)

print(
    "FINAL REFERENCE COMPLETE     : NO"
)

print()

print(
    "EASC IDENTITY                : UNRESOLVED"
)

print(
    "HISTOLOGY SENSITIVITY DESIGN : FROZEN"
)

print(
    "HYPOTHETICAL EASC SCENARIOS  :",
    len(
        patients
    )
)

print()

print(
    "CORE GENE VALUES READ        : NO"
)

print(
    "NORMALIZATION EXECUTED       : NO"
)

print(
    "CORE2 SCORING                : NOT PERFORMED"
)

print(
    "CORE4 SCORING                : NOT PERFORMED"
)

print(
    "PRIMARY CORE2 INFERENCE      : BLOCKED"
)

print()

print(
    "NEW FILES DOWNLOADED         : 0"
)

print(
    "FROZEN STEP19-21 MODIFIED    : NO"
)

print(
    "MANUSCRIPT MODIFIED          : NO"
)

print()

print(
    "NEXT REQUIRED STEP:"
)

print(
    "STEP22A-6A = COMPLETE THE REMAINING"
)

print(
    "15 TUMOR SAMPLES AS TECHNICAL DATA ONLY"
)

print(
    "SO THE FULL OMIX TUMOR REFERENCE CAN BE BUILT."
)

print()

print(
    "HISTOLOGY GATE REMAINS CLOSED."
)

print(
    "DO NOT SCORE CORE2."
)

print()
