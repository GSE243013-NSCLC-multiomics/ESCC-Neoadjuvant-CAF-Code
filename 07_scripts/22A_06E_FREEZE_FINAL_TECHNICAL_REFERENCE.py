#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-6E
#
# FREEZE FINAL TECHNICAL REFERENCE
# BEFORE ANY NORMALIZATION / CORE SCORE
#
# Frozen consequence of >=20 fibroblast rule:
#
#   P15 Before = 18 cells -> TECHNICAL FAIL
#
# Therefore:
#   - P15_Before remains stored in raw pseudobulk provenance
#   - P15_Before is NOT in primary technical reference
#   - P15 is NOT in technically evaluable paired set
#   - threshold is NOT relaxed
#   - NO "borderline rescue" sensitivity arm is invented
#
# Histology sensitivity scenarios are intersected with the
# technical eligibility mask BEFORE any CORE values are read.
#
# This script DOES NOT:
#   - open the pseudobulk RDS
#   - read any gene expression values
#   - calculate CPM
#   - calculate log2(CPM+1)
#   - calculate z-scores
#   - calculate CORE2 / CORE4
#   - perform Before/After inference
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

OBJECT_DIR = (
    PROJECT
    / "03_objects"
    / "OMIX005710"
)

PT_META_FILE = (
    META_DIR
    / "STEP22A_06D_ALL27_FIBROBLAST_PATIENT_TIMEPOINT_METADATA.csv"
)

PAIR_STATUS_FILE = (
    META_DIR
    / "STEP22A_06D_PAIRED_FIBROBLAST_TECHNICAL_STATUS.csv"
)

ZSCORE_LOCK_FILE = (
    META_DIR
    / "STEP22A_05G_ZSCORE_REFERENCE_LOCK.csv"
)

OLD_HIST_SENS_FILE = (
    META_DIR
    / "STEP22A_05G_HISTOLOGY_SENSITIVITY_SCENARIOS.csv"
)

CANONICAL_META_FILE = (
    META_DIR
    / "STEP22A_SAMPLE_METADATA_CANONICAL.csv"
)

GATE_FILE = (
    META_DIR
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

RAW_BUNDLE = (
    OBJECT_DIR
    / "OMIX005710_ALL27_FIBROBLAST_RAW_PSEUDOBULK_HISTOLOGY_UNRESOLVED.rds"
)

TECH_REFERENCE_OUT = (
    META_DIR
    / "STEP22A_06E_FINAL_TECHNICAL_REFERENCE_LOCK.csv"
)

PAIRED_LOCK_OUT = (
    META_DIR
    / "STEP22A_06E_FINAL_PAIRED_TECHNICAL_LOCK.csv"
)

HIST_SENS_OUT = (
    META_DIR
    / "STEP22A_06E_HISTOLOGY_SENSITIVITY_TECHNICAL_INTERSECTION.csv"
)

MANIFEST_OUT = (
    META_DIR
    / "STEP22A_06E_REFERENCE_MEMBERSHIP_MANIFEST.csv"
)

DECISION_OUT = (
    META_DIR
    / "STEP22A_06E_FINAL_TECHNICAL_DECISION_LOCK.txt"
)

SUMMARY_OUT = (
    META_DIR
    / "STEP22A_06E_FINAL_TECHNICAL_REFERENCE_SUMMARY.txt"
)


print()
print("=" * 80)
print("STEP22A-6E: FREEZE FINAL TECHNICAL REFERENCE")
print("=" * 80)
print()


# ============================================================
# Helpers
# ============================================================

def fail(message):

    print()
    print("ERROR:")
    print(message)
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


def patient_number(x):

    try:

        return int(
            clean(x)
            .upper()
            .replace(
                "P",
                ""
            )
        )

    except Exception:

        return 999999


def truth(x):

    return (
        clean(x).upper()
        in {
            "TRUE",
            "T",
            "YES",
            "1"
        }
    )


def membership_hash(items):

    text = "\n".join(
        sorted(
            clean(x)
            for x in items
        )
    ) + "\n"

    return hashlib.sha256(
        text.encode(
            "utf-8"
        )
    ).hexdigest()


def frozen_integrity():

    rows = read_csv(
        FROZEN_INV
    )

    changes = []

    for row in rows:

        raw_path = clean(
            row.get(
                "path"
            )
        )

        if not raw_path:
            continue

        path = Path(
            raw_path
        )

        if not path.is_absolute():

            path = (
                PROJECT
                / path
            )

        if not path.exists():

            changes.append(
                (
                    str(path),
                    "MISSING"
                )
            )

            continue

        expected = clean(
            row.get(
                "size_bytes"
            )
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

                if (
                    expected
                    != actual
                ):

                    changes.append(
                        (
                            str(path),
                            "SIZE_CHANGED"
                        )
                    )

            except Exception:
                pass

    return changes


# ============================================================
# 1. Analysis gate
# ============================================================

print("===== 1. ANALYSIS GATE =====")
print()

gate_rows = read_csv(
    GATE_FILE
)

gates = {
    clean(
        x.get(
            "gate"
        )
    ):
    clean(
        x.get(
            "allowed"
        )
    )
    for x in gate_rows
}


for required in [
    "technical_preprocessing",
    "paired_CORE2_test",
    "manuscript_update"
]:

    if required not in gates:

        fail(
            "Missing analysis gate: "
            + required
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

    fail(
        "Technical preparation gate is closed."
    )


if (
    gates[
        "paired_CORE2_test"
    ]
    == "YES"
):

    fail(
        "CORE2 inference gate unexpectedly open."
    )


if (
    gates[
        "manuscript_update"
    ]
    == "YES"
):

    fail(
        "Manuscript gate unexpectedly open."
    )


print(
    "✓ Technical-reference freezing allowed."
)

print(
    "✓ CORE2 inference remains blocked."
)

print()


# ============================================================
# 2. Frozen Step19-21 integrity
# ============================================================

print("===== 2. FROZEN PROJECT CHECK =====")
print()

changes = frozen_integrity()

print(
    "Frozen original files changed:",
    len(
        changes
    )
)

if changes:

    for item in changes[:20]:
        print(item)

    fail(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )


print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# 3. Verify raw bundle exists — DO NOT OPEN IT
# ============================================================

print("===== 3. RAW PSEUDOBULK BUNDLE FIREWALL =====")
print()


if not RAW_BUNDLE.exists():

    fail(
        "Missing Step22A-6D raw pseudobulk bundle:\n"
        + str(
            RAW_BUNDLE
        )
    )


print(
    "Raw pseudobulk bundle exists : YES"
)

print(
    "Raw pseudobulk bundle opened : NO"
)

print(
    "Gene-level values read       : NO"
)

print()


# ============================================================
# 4. Verify frozen z-score framework
# ============================================================

print("===== 4. VERIFY STEP22A-5G Z-SCORE REFERENCE LOCK =====")
print()


zrows = read_csv(
    ZSCORE_LOCK_FILE
)

zlock = {
    clean(
        x.get(
            "rule"
        )
    ):
    clean(
        x.get(
            "value"
        )
    )
    for x in zrows
}


required_rules = {
    "zscore_reference_dataset":
        "OMIX005710_INTERNAL",

    "zscore_reference_universe":
        "ALL_ELIGIBLE_ESCC_TUMOR_PATIENT_TIMEPOINTS",

    "paired_only_reference":
        "NO",

    "reuse_GSE197677_mean_sd":
        "NO",

    "reuse_GSE221561_mean_sd":
        "NO",

    "primary_score_calculation_now":
        "BLOCKED"
}


for rule, expected in (
    required_rules.items()
):

    observed = zlock.get(
        rule,
        ""
    )

    print(
        f"{rule:<36s}: "
        f"{observed}"
    )

    if observed != expected:

        fail(
            "Frozen z-score rule mismatch:\n"
            + rule
            + "\nExpected: "
            + expected
            + "\nObserved: "
            + observed
        )


print()

print(
    "✓ Frozen normalization/z-score framework unchanged."
)

print()


# ============================================================
# 5. Load final Step22A-6D technical metadata
# ============================================================

print("===== 5. LOAD FINAL TECHNICAL METADATA =====")
print()


pt_rows = read_csv(
    PT_META_FILE
)


if len(
    pt_rows
) != 27:

    fail(
        "Expected exactly 27 patient-timepoint pseudobulks."
    )


pt_ids = [
    clean(
        x.get(
            "patient_timepoint"
        )
    )
    for x in pt_rows
]


if (
    len(
        set(
            pt_ids
        )
    )
    != 27
):

    fail(
        "Duplicate patient-timepoint IDs."
    )


patients_with_tumor = sorted(
    {
        clean(
            x.get(
                "patient_id"
            )
        )
        for x in pt_rows
    },
    key=patient_number
)


print(
    "Patient-timepoint pseudobulks:",
    len(
        pt_rows
    )
)

print(
    "Patients with tumor:",
    len(
        patients_with_tumor
    )
)

print()


# ============================================================
# 6. Freeze >=20 technical reference membership
# ============================================================

print("===== 6. FREEZE >=20 TECHNICAL REFERENCE =====")
print()


technical_pass = []

technical_fail = []


for row in pt_rows:

    status = clean(
        row.get(
            "frozen_technical_status"
        )
    )

    n_fib = int(
        float(
            clean(
                row.get(
                    "n_fibroblast_cells"
                )
            )
        )
    )

    if status == "PASS_GE20":

        if n_fib < 20:

            fail(
                "PASS_GE20 row contains <20 cells."
            )

        technical_pass.append(
            row
        )

    elif status == "FAIL_LT20":

        if n_fib >= 20:

            fail(
                "FAIL_LT20 row contains >=20 cells."
            )

        technical_fail.append(
            row
        )

    else:

        fail(
            "Unexpected technical status: "
            + status
        )


print(
    "Technical PASS >=20:",
    len(
        technical_pass
    ),
    "/ 27"
)

print(
    "Technical FAIL <20:",
    len(
        technical_fail
    )
)

print()


if (
    len(
        technical_pass
    )
    != 26
    or
    len(
        technical_fail
    )
    != 1
):

    fail(
        "Expected final technical structure = 26 PASS / 1 FAIL."
    )


fail_row = technical_fail[0]


fail_patient = clean(
    fail_row.get(
        "patient_id"
    )
)

fail_timepoint = clean(
    fail_row.get(
        "timepoint"
    )
)

fail_pt = clean(
    fail_row.get(
        "patient_timepoint"
    )
)

fail_n = int(
    float(
        clean(
            fail_row.get(
                "n_fibroblast_cells"
            )
        )
    )
)


print(
    "Technical FAIL patient :",
    fail_patient
)

print(
    "Technical FAIL timepoint:",
    fail_timepoint
)

print(
    "Technical FAIL cells   :",
    fail_n
)

print()


if not (
    fail_patient == "P15"
    and
    fail_timepoint == "Before"
    and
    fail_n == 18
):

    fail(
        "Final technical failure is not the expected P15 Before = 18."
    )


print(
    "✓ Frozen threshold consequence confirmed."
)

print()


# ============================================================
# 7. Write final technical-reference lock
#
# P15 Before is retained in provenance/raw data but is NOT
# part of technical reference.
# ============================================================

technical_reference_rows = []


for row in sorted(
    pt_rows,
    key=lambda x: (
        patient_number(
            x.get(
                "patient_id"
            )
        ),
        0
        if clean(
            x.get(
                "timepoint"
            )
        )
        == "Before"
        else 1
    )
):

    n_fib = int(
        float(
            clean(
                row.get(
                    "n_fibroblast_cells"
                )
            )
        )
    )

    pass_qc = (
        clean(
            row.get(
                "frozen_technical_status"
            )
        )
        == "PASS_GE20"
    )


    if pass_qc:

        reference_status = (
            "TECHNICAL_REFERENCE_CANDIDATE"
        )

        exclusion_reason = ""

        retained_for_provenance = "YES"

        normalization_allowed_now = "NO"

        score_allowed_now = "NO"

    else:

        reference_status = (
            "EXCLUDED_FROM_TECHNICAL_REFERENCE"
        )

        exclusion_reason = (
            "FROZEN_FIBROBLAST_QC_LT20"
        )

        retained_for_provenance = "YES"

        normalization_allowed_now = "NO"

        score_allowed_now = "NO"


    technical_reference_rows.append({
        "patient_id":
            clean(
                row.get(
                    "patient_id"
                )
            ),

        "timepoint":
            clean(
                row.get(
                    "timepoint"
                )
            ),

        "patient_timepoint":
            clean(
                row.get(
                    "patient_timepoint"
                )
            ),

        "n_fibroblast_cells":
            n_fib,

        "frozen_technical_status":
            clean(
                row.get(
                    "frozen_technical_status"
                )
            ),

        "technical_reference_status":
            reference_status,

        "exclusion_reason":
            exclusion_reason,

        "retained_in_raw_bundle":
            retained_for_provenance,

        "borderline_rescue_branch":
            "NO",

        "threshold_relaxed":
            "NO",

        "histology_status":
            "UNRESOLVED",

        "final_ESCC_reference_eligible":
            "NO_HISTOLOGY_GATE_BLOCKED",

        "normalization_allowed_now":
            normalization_allowed_now,

        "CORE_scoring_allowed_now":
            score_allowed_now
    })


write_csv(
    TECH_REFERENCE_OUT,
    technical_reference_rows,
    [
        "patient_id",
        "timepoint",
        "patient_timepoint",
        "n_fibroblast_cells",
        "frozen_technical_status",
        "technical_reference_status",
        "exclusion_reason",
        "retained_in_raw_bundle",
        "borderline_rescue_branch",
        "threshold_relaxed",
        "histology_status",
        "final_ESCC_reference_eligible",
        "normalization_allowed_now",
        "CORE_scoring_allowed_now"
    ]
)


technical_reference_ids = [
    clean(
        x.get(
            "patient_timepoint"
        )
    )
    for x in technical_pass
]


technical_reference_sha = membership_hash(
    technical_reference_ids
)


print(
    "Technical-reference candidate PTs:",
    len(
        technical_reference_ids
    )
)

print(
    "Membership SHA256:",
    technical_reference_sha
)

print()


# ============================================================
# 8. Freeze paired technical candidate set
# ============================================================

print("===== 7. FREEZE PAIRED TECHNICAL SET =====")
print()


pair_rows = read_csv(
    PAIR_STATUS_FILE
)


complete_pairs = [
    x
    for x in pair_rows
    if truth(
        x.get(
            "complete_tumor_pair"
        )
    )
]


technical_pairs = [
    x
    for x in complete_pairs
    if truth(
        x.get(
            "technically_evaluable_pair"
        )
    )
]


complete_pair_ids = sorted(
    [
        clean(
            x.get(
                "patient_id"
            )
        )
        for x in complete_pairs
    ],
    key=patient_number
)


technical_pair_ids = sorted(
    [
        clean(
            x.get(
                "patient_id"
            )
        )
        for x in technical_pairs
    ],
    key=patient_number
)


expected_complete_pairs = [
    "P6",
    "P9",
    "P15",
    "P16",
    "P18",
    "P19"
]


expected_technical_pairs = [
    "P6",
    "P9",
    "P16",
    "P18",
    "P19"
]


print(
    "Complete paired patients:",
    len(
        complete_pair_ids
    )
)

print(
    "IDs:",
    ", ".join(
        complete_pair_ids
    )
)

print()

print(
    "Technically evaluable paired patients:",
    len(
        technical_pair_ids
    )
)

print(
    "IDs:",
    ", ".join(
        technical_pair_ids
    )
)

print()


if (
    complete_pair_ids
    != expected_complete_pairs
):

    fail(
        "Unexpected complete-pair membership."
    )


if (
    technical_pair_ids
    != expected_technical_pairs
):

    fail(
        "Unexpected technically evaluable paired membership."
    )


paired_lock_rows = []


for patient in complete_pair_ids:

    eligible = (
        patient
        in technical_pair_ids
    )

    if eligible:

        status = (
            "PAIRED_TECHNICAL_CANDIDATE"
        )

        reason = ""

    else:

        status = (
            "EXCLUDED_FROM_PAIRED_TECHNICAL_SET"
        )

        reason = (
            "BEFORE_FIBROBLAST_QC_LT20"
            if patient == "P15"
            else
            "TECHNICAL_QC_FAILURE"
        )


    paired_lock_rows.append({
        "patient_id":
            patient,

        "complete_tumor_pair":
            "YES",

        "technically_evaluable_pair":
            (
                "YES"
                if eligible
                else "NO"
            ),

        "paired_candidate_status":
            status,

        "exclusion_reason":
            reason,

        "threshold_relaxed":
            "NO",

        "borderline_rescue_branch":
            "NO",

        "histology_status":
            "UNRESOLVED",

        "primary_paired_inference_allowed_now":
            "NO"
    })


write_csv(
    PAIRED_LOCK_OUT,
    paired_lock_rows,
    [
        "patient_id",
        "complete_tumor_pair",
        "technically_evaluable_pair",
        "paired_candidate_status",
        "exclusion_reason",
        "threshold_relaxed",
        "borderline_rescue_branch",
        "histology_status",
        "primary_paired_inference_allowed_now"
    ]
)


paired_membership_sha = membership_hash(
    technical_pair_ids
)


print(
    "Paired technical membership SHA256:",
    paired_membership_sha
)

print()

print(
    "✓ P15 excluded from paired technical candidate set."
)

print(
    "✓ No threshold relaxation."
)

print(
    "✓ No post-hoc borderline rescue branch."
)

print()


# ============================================================
# 9. Verify all 22 source-study patients
# ============================================================

print("===== 8. SOURCE PATIENT SET =====")
print()


canonical_rows = read_csv(
    CANONICAL_META_FILE
)


all_source_patients = sorted(
    {
        clean(
            x.get(
                "patient_id"
            )
        )
        for x in canonical_rows
        if clean(
            x.get(
                "patient_id"
            )
        )
    },
    key=patient_number
)


print(
    "Source-study patients:",
    len(
        all_source_patients
    )
)

print(
    "IDs:",
    ", ".join(
        all_source_patients
    )
)

print()


if len(
    all_source_patients
) != 22:

    fail(
        "Expected exactly 22 source-study patients."
    )


# ============================================================
# 10. Verify prior histology sensitivity design exists
# ============================================================

print("===== 9. VERIFY PRE-FROZEN HISTOLOGY SENSITIVITY DESIGN =====")
print()


old_sens = read_csv(
    OLD_HIST_SENS_FILE
)


hypothetical_rows = [
    x
    for x in old_sens
    if clean(
        x.get(
            "hypothetical_EASC_patient"
        )
    )
    not in {
        "",
        "UNKNOWN"
    }
]


old_hypothetical_patients = sorted(
    {
        clean(
            x.get(
                "hypothetical_EASC_patient"
            )
        )
        for x in hypothetical_rows
    },
    key=patient_number
)


if (
    old_hypothetical_patients
    != all_source_patients
):

    fail(
        "Step22A-5G histology sensitivity design "
        "does not contain the same 22 patients."
    )


print(
    "Pre-frozen hypothetical EASC identities:",
    len(
        old_hypothetical_patients
    )
)

print(
    "✓ All 22 identities were frozen before CORE values."
)

print()


# ============================================================
# 11. Intersect histology scenarios with technical eligibility
#
# Baseline technical reference:
#   26 PASS patient-timepoints.
#
# For each hypothetical EASC identity:
#   remove ALL technically-passing tumor observations of that
#   patient from the technical reference candidate.
#
# For paired analysis:
#   begin with the fixed 5 technically evaluable pairs.
#   if hypothetical EASC is one of those 5, remove that pair.
#
# P15 is already absent from paired candidate due to QC.
#
# No scores are calculated.
# ============================================================

print("===== 10. FREEZE HISTOLOGY × TECHNICAL INTERSECTION =====")
print()


baseline_reference = sorted(
    technical_reference_ids
)

baseline_pairs = sorted(
    technical_pair_ids,
    key=patient_number
)


scenario_rows = []


# Unresolved baseline
scenario_rows.append({
    "scenario_id":
        "S00_UNRESOLVED",

    "hypothetical_EASC_patient":
        "UNKNOWN",

    "technical_reference_n":
        len(
            baseline_reference
        ),

    "technical_reference_membership_SHA256":
        membership_hash(
            baseline_reference
        ),

    "technical_reference_removed":
        "",

    "paired_candidate_n":
        len(
            baseline_pairs
        ),

    "paired_candidate_membership_SHA256":
        membership_hash(
            baseline_pairs
        ),

    "paired_candidate_removed":
        "",

    "P15_before_status":
        "ALREADY_EXCLUDED_TECHNICAL_QC_LT20",

    "normalization_allowed":
        "NO",

    "CORE_scoring_allowed":
        "NO",

    "primary_inference_allowed":
        "NO",

    "reason":
        "EASC_IDENTITY_UNRESOLVED"
})


for patient in all_source_patients:

    patient_ref_members = sorted(
        [
            clean(
                x.get(
                    "patient_timepoint"
                )
            )
            for x in technical_pass
            if clean(
                x.get(
                    "patient_id"
                )
            )
            == patient
        ]
    )


    scenario_reference = sorted(
        [
            x
            for x in baseline_reference
            if x not in patient_ref_members
        ]
    )


    if patient in baseline_pairs:

        scenario_pairs = [
            x
            for x in baseline_pairs
            if x != patient
        ]

        removed_pair = patient

    else:

        scenario_pairs = list(
            baseline_pairs
        )

        removed_pair = ""


    scenario_rows.append({
        "scenario_id":
            "S_HYPOTHETICAL_EASC_"
            + patient,

        "hypothetical_EASC_patient":
            patient,

        "technical_reference_n":
            len(
                scenario_reference
            ),

        "technical_reference_membership_SHA256":
            membership_hash(
                scenario_reference
            ),

        "technical_reference_removed":
            ";".join(
                patient_ref_members
            ),

        "paired_candidate_n":
            len(
                scenario_pairs
            ),

        "paired_candidate_membership_SHA256":
            membership_hash(
                scenario_pairs
            ),

        "paired_candidate_removed":
            removed_pair,

        "P15_before_status":
            "ALREADY_EXCLUDED_TECHNICAL_QC_LT20",

        "normalization_allowed":
            "NO",

        "CORE_scoring_allowed":
            "NO",

        "primary_inference_allowed":
            "NO",

        "reason":
            "PRE_SPECIFIED_HISTOLOGY_SENSITIVITY_ONLY"
    })


write_csv(
    HIST_SENS_OUT,
    scenario_rows,
    [
        "scenario_id",
        "hypothetical_EASC_patient",
        "technical_reference_n",
        "technical_reference_membership_SHA256",
        "technical_reference_removed",
        "paired_candidate_n",
        "paired_candidate_membership_SHA256",
        "paired_candidate_removed",
        "P15_before_status",
        "normalization_allowed",
        "CORE_scoring_allowed",
        "primary_inference_allowed",
        "reason"
    ]
)


print(
    "Histology × technical scenarios:",
    len(
        scenario_rows
    )
)

print(
    "  unresolved baseline : 1"
)

print(
    "  hypothetical EASC   : 22"
)

print()


# ============================================================
# 12. Special audit: what if EASC = P15?
#
# P15 Before is already excluded technically.
# P15 After, if technically PASS, would be removed from
# histology-sensitive reference if P15 were the EASC patient.
#
# P15 does NOT re-enter paired set under any scenario.
# ============================================================

print("===== 11. P15 SPECIAL-CASE AUDIT =====")
print()


p15_pass_pt = sorted(
    [
        clean(
            x.get(
                "patient_timepoint"
            )
        )
        for x in technical_pass
        if clean(
            x.get(
                "patient_id"
            )
        )
        == "P15"
    ]
)


print(
    "P15 technical-reference PASS observations:",
    (
        ", ".join(
            p15_pass_pt
        )
        if p15_pass_pt
        else "<none>"
    )
)

print(
    "P15 Before technical status:",
    "EXCLUDED_LT20"
)

print(
    "P15 paired technical status:",
    "EXCLUDED"
)

print(
    "P15 borderline-rescue sensitivity:",
    "NONE"
)

print()


# ============================================================
# 13. Membership manifest
# ============================================================

manifest_rows = [
    {
        "item":
            "technical_threshold",

        "value":
            "FIBROBLAST_CELLS_GE20"
    },

    {
        "item":
            "technical_reference_candidate_n",

        "value":
            str(
                len(
                    baseline_reference
                )
            )
    },

    {
        "item":
            "technical_reference_membership_SHA256",

        "value":
            technical_reference_sha
    },

    {
        "item":
            "technical_fail_n",

        "value":
            "1"
    },

    {
        "item":
            "technical_fail_patient_timepoint",

        "value":
            fail_pt
    },

    {
        "item":
            "technical_fail_fibroblast_cells",

        "value":
            str(
                fail_n
            )
    },

    {
        "item":
            "complete_paired_patients_n",

        "value":
            str(
                len(
                    complete_pair_ids
                )
            )
    },

    {
        "item":
            "technical_paired_candidate_n",

        "value":
            str(
                len(
                    baseline_pairs
                )
            )
    },

    {
        "item":
            "technical_paired_candidate_ids",

        "value":
            ";".join(
                baseline_pairs
            )
    },

    {
        "item":
            "technical_paired_membership_SHA256",

        "value":
            paired_membership_sha
    },

    {
        "item":
            "P15_paired_status",

        "value":
            "EXCLUDED_BEFORE_LT20"
    },

    {
        "item":
            "threshold_relaxed",

        "value":
            "NO"
    },

    {
        "item":
            "borderline_rescue_branch_created",

        "value":
            "NO"
    },

    {
        "item":
            "raw_P15_before_retained",

        "value":
            "YES_PROVENANCE_ONLY"
    },

    {
        "item":
            "histology_scenarios",

        "value":
            "23_TOTAL_1_UNRESOLVED_PLUS_22_HYPOTHETICAL"
    },

    {
        "item":
            "histology_mapping",

        "value":
            "INCOMPLETE"
    },

    {
        "item":
            "normalization_executed",

        "value":
            "NO"
    },

    {
        "item":
            "CORE_gene_values_read",

        "value":
            "NO"
    },

    {
        "item":
            "CORE2_scored",

        "value":
            "NO"
    },

    {
        "item":
            "CORE4_scored",

        "value":
            "NO"
    },

    {
        "item":
            "primary_inference_allowed",

        "value":
            "NO"
    },

    {
        "item":
            "lock_status",

        "value":
            "FINAL_TECHNICAL_REFERENCE_FROZEN"
    }
]


write_csv(
    MANIFEST_OUT,
    manifest_rows,
    [
        "item",
        "value"
    ]
)


# ============================================================
# 14. Human-readable decision lock
# ============================================================

with open(
    DECISION_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-6E FINAL TECHNICAL REFERENCE DECISION LOCK\n"
    )

    f.write(
        "=" * 72
        + "\n\n"
    )

    f.write(
        "FROZEN TECHNICAL QC RULE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Fibroblast_CAF cells >=20 : PASS\n"
    )

    f.write(
        "Fibroblast_CAF cells <20  : FAIL\n"
    )

    f.write(
        "Threshold relaxation       : FORBIDDEN\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "P15 BEFORE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Fibroblast_CAF cells : 18\n"
    )

    f.write(
        "Technical status     : FAIL_LT20\n"
    )

    f.write(
        "Raw pseudobulk kept  : YES\n"
    )

    f.write(
        "Technical reference  : EXCLUDED\n"
    )

    f.write(
        "Paired primary set   : EXCLUDED\n"
    )

    f.write(
        "Borderline rescue    : NONE\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "RATIONALE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "The >=20 threshold was frozen before the final all-27\n"
    )

    f.write(
        "annotation result. It is not changed because P15_Before\n"
    )

    f.write(
        "happens to contain 18 fibroblast cells.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Creating a special 18-cell rescue branch after observing\n"
    )

    f.write(
        "the failure would weaken the pre-specified QC firewall.\n"
    )

    f.write(
        "Therefore P15_Before remains available for provenance\n"
    )

    f.write(
        "but does not enter normalization, z-scoring or inference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "FINAL TECHNICAL REFERENCE CANDIDATE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        f"Patient-timepoints: {len(baseline_reference)} / 27\n"
    )

    f.write(
        f"Membership SHA256: {technical_reference_sha}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "FINAL PAIRED TECHNICAL CANDIDATE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Complete source-study pairs: 6\n"
    )

    f.write(
        "Technically evaluable pairs: 5\n"
    )

    f.write(
        "Patients: "
        + ", ".join(
            baseline_pairs
        )
        + "\n"
    )

    f.write(
        "P15 excluded because Before has 18 fibroblasts.\n"
    )

    f.write(
        f"Membership SHA256: {paired_membership_sha}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "HISTOLOGY GATE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "EASC identity: UNRESOLVED\n"
    )

    f.write(
        "The 26 technical-pass observations are only candidates.\n"
    )

    f.write(
        "They are NOT yet a final ESCC reference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "For each of the 22 pre-frozen hypothetical EASC identities,\n"
    )

    f.write(
        "technical eligibility is applied first, then all technical-\n"
    )

    f.write(
        "pass tumor observations of that hypothetical EASC patient\n"
    )

    f.write(
        "are removed from that sensitivity reference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "The paired candidate begins with the fixed 5 technical pairs.\n"
    )

    f.write(
        "If hypothetical EASC is one of those 5, that pair is removed.\n"
    )

    f.write(
        "P15 never re-enters the paired set under any scenario.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CURRENT ANALYTIC FIREWALL\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Normalization executed : NO\n"
    )

    f.write(
        "CORE gene values read  : NO\n"
    )

    f.write(
        "CORE2 calculated       : NO\n"
    )

    f.write(
        "CORE4 calculated       : NO\n"
    )

    f.write(
        "Paired inference       : NO\n"
    )

    f.write(
        "Primary inference gate : BLOCKED\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "LOCK STATUS:\n"
    )

    f.write(
        "FINAL_TECHNICAL_REFERENCE_FROZEN\n"
    )


# ============================================================
# 15. Final frozen project check
# ============================================================

print("===== 12. FINAL FROZEN PROJECT CHECK =====")
print()

changes_after = frozen_integrity()

print(
    "Frozen original files changed:",
    len(
        changes_after
    )
)

if changes_after:

    for item in changes_after[:20]:
        print(item)

    fail(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )


print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# 16. Summary
# ============================================================

with open(
    SUMMARY_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-6E FINAL TECHNICAL REFERENCE SUMMARY\n"
    )

    f.write(
        "=" * 68
        + "\n\n"
    )

    f.write(
        "Frozen threshold: >=20 fibroblast cells\n"
    )

    f.write(
        "Technical PASS patient-timepoints: 26/27\n"
    )

    f.write(
        "Technical FAIL patient-timepoints: 1/27\n"
    )

    f.write(
        "Technical failure: P15 Before = 18 cells\n"
    )

    f.write(
        "P15 Before raw pseudobulk retained: YES\n"
    )

    f.write(
        "P15 Before in technical reference: NO\n"
    )

    f.write(
        "Threshold relaxed: NO\n"
    )

    f.write(
        "Borderline rescue branch: NO\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Complete paired patients: 6\n"
    )

    f.write(
        "Technically evaluable paired patients: 5\n"
    )

    f.write(
        "Paired technical candidates: "
        + ", ".join(
            baseline_pairs
        )
        + "\n"
    )

    f.write(
        "P15 paired status: EXCLUDED\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Histology sensitivity scenarios: 23\n"
    )

    f.write(
        "EASC identity: UNRESOLVED\n"
    )

    f.write(
        "Final ESCC reference: BLOCKED\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Normalization executed: NO\n"
    )

    f.write(
        "CORE gene values read: NO\n"
    )

    f.write(
        "CORE2 scoring: NO\n"
    )

    f.write(
        "CORE4 scoring: NO\n"
    )

    f.write(
        "Primary paired inference: BLOCKED\n"
    )


# ============================================================
# FINAL
# ============================================================

print("=" * 80)
print("STEP22A-6E COMPLETE")
print("=" * 80)
print()

print(
    "FROZEN FIBROBLAST THRESHOLD  : >=20"
)

print()

print(
    "TECHNICAL REFERENCE CANDIDATE : 26 / 27"
)

print(
    "TECHNICAL FAIL                : 1 / 27"
)

print(
    "FAILED PATIENT-TIMEPOINT      : P15__Before"
)

print(
    "FAILED FIBROBLAST CELLS       : 18"
)

print()

print(
    "P15 BEFORE RAW DATA RETAINED  : YES"
)

print(
    "P15 BEFORE IN REFERENCE       : NO"
)

print(
    "THRESHOLD RELAXED             : NO"
)

print(
    "BORDERLINE RESCUE BRANCH      : NO"
)

print()

print(
    "COMPLETE PAIRED PATIENTS      : 6"
)

print(
    "PAIRED TECHNICAL CANDIDATES   : 5"
)

print(
    "PAIRED IDS                    : "
    + ", ".join(
        baseline_pairs
    )
)

print(
    "P15 IN PAIRED CANDIDATE       : NO"
)

print()

print(
    "TECHNICAL REFERENCE SHA256:"
)

print(
    technical_reference_sha
)

print()

print(
    "PAIRED MEMBERSHIP SHA256:"
)

print(
    paired_membership_sha
)

print()

print(
    "HISTOLOGY SCENARIOS           : 23"
)

print(
    "  UNRESOLVED BASELINE          : 1"
)

print(
    "  HYPOTHETICAL EASC            : 22"
)

print()

print(
    "EASC IDENTITY                 : UNRESOLVED"
)

print(
    "FINAL ESCC REFERENCE          : BLOCKED"
)

print()

print(
    "RAW PSEUDOBULK BUNDLE OPENED  : NO"
)

print(
    "CORE GENE VALUES READ         : NO"
)

print(
    "NORMALIZATION                 : NOT PERFORMED"
)

print(
    "CPM                           : NOT CALCULATED"
)

print(
    "log2(CPM+1)                   : NOT CALCULATED"
)

print(
    "Z-SCORE                       : NOT CALCULATED"
)

print(
    "CORE2 SCORING                 : NOT PERFORMED"
)

print(
    "CORE4 SCORING                 : NOT PERFORMED"
)

print(
    "BEFORE/AFTER INFERENCE        : NOT PERFORMED"
)

print()

print(
    "FROZEN STEP19-21 MODIFIED     : NO"
)

print(
    "MANUSCRIPT MODIFIED           : NO"
)

print()

print(
    "LOCK STATUS:"
)

print(
    "FINAL_TECHNICAL_REFERENCE_FROZEN"
)

print()

print(
    "NEXT:"
)

print(
    "DO NOT SCORE CORE2 WHILE EASC IDENTITY IS UNRESOLVED."
)

print()

print(
    "THE TECHNICAL PIPELINE IS NOW FROZEN."
)

print(
    "THE REMAINING BLOCKER IS HISTOLOGY, NOT COMPUTATION."
)

print()
