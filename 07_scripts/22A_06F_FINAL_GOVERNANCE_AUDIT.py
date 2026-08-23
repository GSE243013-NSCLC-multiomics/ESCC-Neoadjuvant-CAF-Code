#!/usr/bin/env python3

# ============================================================
# STEP22A-6F
# FINAL GOVERNANCE / LOCK AUDIT
#
# PURPOSE:
#   audit all frozen Step22A decisions before any CORE scoring
#
# STRICT FIREWALL:
#   - DO NOT open any RDS
#   - DO NOT read expression matrices
#   - DO NOT read CORE values
#   - DO NOT normalize
#   - DO NOT score CORE2 / CORE4
#   - DO NOT run inference
#
# Expected frozen state:
#   - technical tumor candidate PTs = 26/27
#   - only failure = P15 Before, 18 fibroblasts
#   - technical paired candidates = P6/P9/P16/P18/P19
#   - EASC identity unresolved
#   - histology sensitivity scenarios = 23
#   - primary inference blocked
# ============================================================

from pathlib import Path
import csv
import hashlib
import sys


PROJECT = Path.cwd()

META = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

OBJECT_DIR = (
    PROJECT
    / "03_objects"
    / "OMIX005710"
)

GATE_FILE = (
    META
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

TECH_REFERENCE = (
    META
    / "STEP22A_06E_FINAL_TECHNICAL_REFERENCE_LOCK.csv"
)

PAIRED_LOCK = (
    META
    / "STEP22A_06E_FINAL_PAIRED_TECHNICAL_LOCK.csv"
)

TECH_MANIFEST = (
    META
    / "STEP22A_06E_REFERENCE_MEMBERSHIP_MANIFEST.csv"
)

HIST_TECH = (
    META
    / "STEP22A_06E_HISTOLOGY_SENSITIVITY_TECHNICAL_INTERSECTION.csv"
)

ZSCORE_LOCK = (
    META
    / "STEP22A_05G_ZSCORE_REFERENCE_LOCK.csv"
)

SCORING_LOCK = (
    META
    / "STEP22A_05F_EXACT_SCORING_IMPLEMENTATION_LOCK.txt"
)

MARKER_STATUS = (
    META
    / "STEP22A_06B1_ALIGNMENT_FINAL_STATUS.csv"
)

ANTI_CIRC = (
    META
    / "STEP22A_06C_ANTI_CIRCULARITY_AUDIT.csv"
)

PB_PROVENANCE = (
    META
    / "STEP22A_06D_PSEUDOBULK_PROVENANCE.csv"
)

RAW_BUNDLE = (
    OBJECT_DIR
    / "OMIX005710_ALL27_FIBROBLAST_RAW_PSEUDOBULK_HISTOLOGY_UNRESOLVED.rds"
)

FROZEN_INV = (
    META
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

AUDIT_OUT = (
    META
    / "STEP22A_06F_FINAL_GOVERNANCE_AUDIT.csv"
)

DECISION_OUT = (
    META
    / "STEP22A_06F_FINAL_DECISION_DOSSIER.txt"
)


print()
print("=" * 80)
print("STEP22A-6F: FINAL GOVERNANCE / LOCK AUDIT")
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
    return str(x or "").replace("\r", "").strip()


def read_csv(path):

    if not path.exists():
        fail(
            "Missing required lock file:\n"
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


def write_csv(path, rows, fields):

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


def dict_from_rows(rows, key, value):

    return {
        clean(x.get(key)):
            clean(x.get(value))
        for x in rows
    }


def truth(x):

    return clean(x).upper() in {
        "YES",
        "TRUE",
        "T",
        "1"
    }


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


checks = []


def add_check(
    name,
    passed,
    observed,
    expected,
    severity="HARD"
):

    checks.append({
        "check":
            name,

        "status":
            (
                "PASS"
                if passed
                else "FAIL"
            ),

        "observed":
            str(observed),

        "expected":
            str(expected),

        "severity":
            severity
    })


# ============================================================
# 1. Required artifacts
# ============================================================

print("===== 1. REQUIRED LOCK ARTIFACTS =====")
print()

required_files = [
    GATE_FILE,
    TECH_REFERENCE,
    PAIRED_LOCK,
    TECH_MANIFEST,
    HIST_TECH,
    ZSCORE_LOCK,
    SCORING_LOCK,
    MARKER_STATUS,
    ANTI_CIRC,
    PB_PROVENANCE,
    FROZEN_INV,
    RAW_BUNDLE
]


for path in required_files:

    exists = path.exists()

    print(
        f"{path.name:<72s}",
        "FOUND"
        if exists
        else "MISSING"
    )

    add_check(
        "artifact_"
        + path.name,
        exists,
        "FOUND"
        if exists
        else "MISSING",
        "FOUND"
    )


if not all(
    x["status"] == "PASS"
    for x in checks
):

    fail(
        "One or more required Step22A lock artifacts are missing."
    )


print()

print(
    "IMPORTANT: raw pseudobulk RDS existence verified only."
)

print(
    "RDS opened: NO"
)

print()


# ============================================================
# 2. Gate audit
# ============================================================

print("===== 2. ANALYSIS GATES =====")
print()

gates = dict_from_rows(
    read_csv(
        GATE_FILE
    ),
    "gate",
    "allowed"
)


for key in [
    "technical_preprocessing",
    "paired_CORE2_test",
    "manuscript_update"
]:

    print(
        f"{key:<32s}:",
        gates.get(
            key,
            "<missing>"
        )
    )


add_check(
    "technical_preprocessing_gate",
    gates.get(
        "technical_preprocessing"
    )
    == "YES",
    gates.get(
        "technical_preprocessing"
    ),
    "YES"
)


add_check(
    "paired_CORE2_gate",
    gates.get(
        "paired_CORE2_test"
    )
    != "YES",
    gates.get(
        "paired_CORE2_test"
    ),
    "NOT YES"
)


add_check(
    "manuscript_gate",
    gates.get(
        "manuscript_update"
    )
    != "YES",
    gates.get(
        "manuscript_update"
    ),
    "NOT YES"
)


print()


# ============================================================
# 3. Technical reference audit
# ============================================================

print("===== 3. FINAL TECHNICAL REFERENCE =====")
print()

tech = read_csv(
    TECH_REFERENCE
)


candidates = [
    x
    for x in tech
    if clean(
        x.get(
            "technical_reference_status"
        )
    )
    == "TECHNICAL_REFERENCE_CANDIDATE"
]


excluded = [
    x
    for x in tech
    if clean(
        x.get(
            "technical_reference_status"
        )
    )
    == "EXCLUDED_FROM_TECHNICAL_REFERENCE"
]


print(
    "Technical candidates:",
    len(
        candidates
    )
)

print(
    "Technical exclusions:",
    len(
        excluded
    )
)

print()


add_check(
    "technical_candidate_n",
    len(
        candidates
    )
    == 26,
    len(
        candidates
    ),
    26
)


add_check(
    "technical_exclusion_n",
    len(
        excluded
    )
    == 1,
    len(
        excluded
    ),
    1
)


if len(
    excluded
) == 1:

    ex = excluded[0]

    ex_patient = clean(
        ex.get(
            "patient_id"
        )
    )

    ex_time = clean(
        ex.get(
            "timepoint"
        )
    )

    ex_cells = clean(
        ex.get(
            "n_fibroblast_cells"
        )
    )

    ex_reason = clean(
        ex.get(
            "exclusion_reason"
        )
    )

    print(
        "Excluded:",
        ex_patient,
        ex_time,
        ex_cells,
        "cells"
    )

    print(
        "Reason:",
        ex_reason
    )

    print()


    add_check(
        "technical_failure_identity",
        (
            ex_patient == "P15"
            and
            ex_time == "Before"
            and
            ex_cells == "18"
        ),
        f"{ex_patient}/{ex_time}/{ex_cells}",
        "P15/Before/18"
    )


    add_check(
        "technical_failure_reason",
        ex_reason
        == "FROZEN_FIBROBLAST_QC_LT20",
        ex_reason,
        "FROZEN_FIBROBLAST_QC_LT20"
    )


    add_check(
        "P15_raw_retained",
        clean(
            ex.get(
                "retained_in_raw_bundle"
            )
        )
        == "YES",
        clean(
            ex.get(
                "retained_in_raw_bundle"
            )
        ),
        "YES"
    )


    add_check(
        "threshold_not_relaxed",
        clean(
            ex.get(
                "threshold_relaxed"
            )
        )
        == "NO",
        clean(
            ex.get(
                "threshold_relaxed"
            )
        ),
        "NO"
    )


    add_check(
        "no_borderline_rescue",
        clean(
            ex.get(
                "borderline_rescue_branch"
            )
        )
        == "NO",
        clean(
            ex.get(
                "borderline_rescue_branch"
            )
        ),
        "NO"
    )


# ============================================================
# 4. Paired-set audit
# ============================================================

print("===== 4. FINAL PAIRED TECHNICAL SET =====")
print()

pair_rows = read_csv(
    PAIRED_LOCK
)


pair_candidates = sorted(
    [
        clean(
            x.get(
                "patient_id"
            )
        )
        for x in pair_rows
        if clean(
            x.get(
                "paired_candidate_status"
            )
        )
        == "PAIRED_TECHNICAL_CANDIDATE"
    ],
    key=lambda x: int(
        x[1:]
    )
)


expected_pairs = [
    "P6",
    "P9",
    "P16",
    "P18",
    "P19"
]


print(
    "Paired candidates:",
    len(
        pair_candidates
    )
)

print(
    "IDs:",
    ", ".join(
        pair_candidates
    )
)

print()


add_check(
    "paired_candidate_n",
    len(
        pair_candidates
    )
    == 5,
    len(
        pair_candidates
    ),
    5
)


add_check(
    "paired_candidate_membership",
    pair_candidates
    == expected_pairs,
    ",".join(
        pair_candidates
    ),
    ",".join(
        expected_pairs
    )
)


p15_pair = [
    x
    for x in pair_rows
    if clean(
        x.get(
            "patient_id"
        )
    )
    == "P15"
]


add_check(
    "P15_pair_excluded",
    (
        len(
            p15_pair
        )
        == 1
        and
        clean(
            p15_pair[0].get(
                "paired_candidate_status"
            )
        )
        == "EXCLUDED_FROM_PAIRED_TECHNICAL_SET"
    ),
    (
        clean(
            p15_pair[0].get(
                "paired_candidate_status"
            )
        )
        if p15_pair
        else "<missing>"
    ),
    "EXCLUDED_FROM_PAIRED_TECHNICAL_SET"
)


# ============================================================
# 5. Histology sensitivity audit
# ============================================================

print("===== 5. HISTOLOGY SENSITIVITY LOCK =====")
print()

hist = read_csv(
    HIST_TECH
)


baseline = [
    x
    for x in hist
    if clean(
        x.get(
            "scenario_id"
        )
    )
    == "S00_UNRESOLVED"
]


hyp = [
    x
    for x in hist
    if clean(
        x.get(
            "scenario_id"
        )
    ).startswith(
        "S_HYPOTHETICAL_EASC_"
    )
]


print(
    "Unresolved baseline:",
    len(
        baseline
    )
)

print(
    "Hypothetical EASC scenarios:",
    len(
        hyp
    )
)

print(
    "Total:",
    len(
        hist
    )
)

print()


add_check(
    "histology_scenario_total",
    len(
        hist
    )
    == 23,
    len(
        hist
    ),
    23
)


add_check(
    "histology_unresolved_baseline_n",
    len(
        baseline
    )
    == 1,
    len(
        baseline
    ),
    1
)


add_check(
    "histology_hypothetical_n",
    len(
        hyp
    )
    == 22,
    len(
        hyp
    ),
    22
)


all_hist_blocked = all(
    clean(
        x.get(
            "primary_inference_allowed"
        )
    )
    == "NO"
    for x in hist
)


add_check(
    "histology_scenarios_not_primary",
    all_hist_blocked,
    "ALL NO"
    if all_hist_blocked
    else "ONE OR MORE NOT NO",
    "ALL NO"
)


# ============================================================
# 6. Z-score rule audit
# ============================================================

print("===== 6. FROZEN SCORE-REFERENCE RULES =====")
print()

zlock = dict_from_rows(
    read_csv(
        ZSCORE_LOCK
    ),
    "rule",
    "value"
)


z_expected = {
    "expression_transform":
        "log2((raw_count/library_size)*1e6 + 1)",

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

    "weight_refitting":
        "FORBIDDEN",

    "cutoff_optimization":
        "FORBIDDEN",

    "current_histology_status":
        "INCOMPLETE",

    "primary_score_calculation_now":
        "BLOCKED",

    "primary_inference_now":
        "BLOCKED"
}


for key, expected in (
    z_expected.items()
):

    observed = zlock.get(
        key,
        "<missing>"
    )

    print(
        f"{key:<36s}: "
        f"{observed}"
    )

    add_check(
        "zlock_" + key,
        observed == expected,
        observed,
        expected
    )


print()


# ============================================================
# 7. Structural marker audit
# ============================================================

print("===== 7. STRUCTURAL MARKER LOCK =====")
print()

marker = dict_from_rows(
    read_csv(
        MARKER_STATUS
    ),
    "item",
    "value"
)


marker_expected = {
    "frozen_broad_marker_dictionary":
        "33_RETAINED",

    "joint_measurable_broad_markers":
        "31",

    "structural_absence_interpretation":
        "NOT_MEASURED_NA_NEVER_ZERO",

    "sample_exclusion_due_to_KRT_missingness":
        "NONE",

    "CORE_genes_all27":
        "4_OF_4_PRESENT",

    "CORE_scoring_impact":
        "NONE_TECHNICALLY"
}


for key, expected in (
    marker_expected.items()
):

    observed = marker.get(
        key,
        "<missing>"
    )

    print(
        f"{key:<44s}: "
        f"{observed}"
    )

    add_check(
        "marker_" + key,
        observed == expected,
        observed,
        expected
    )


print()


# ============================================================
# 8. Anti-circularity audit
# ============================================================

print("===== 8. ANTI-CIRCULARITY LOCK =====")
print()

anti = dict_from_rows(
    read_csv(
        ANTI_CIRC
    ),
    "item",
    "status"
)


anti_expected = {
    "CORE_genes_removed_before_QC":
        "YES",

    "CORE_genes_removed_before_normalization":
        "YES",

    "CORE_genes_removed_before_HVG":
        "YES",

    "CORE_genes_removed_before_PCA":
        "YES",

    "CORE_genes_removed_before_clustering":
        "YES",

    "CORE_genes_removed_before_annotation":
        "YES",

    "CORE2_scoring_performed":
        "NO",

    "CORE4_scoring_performed":
        "NO",

    "treatment_used_for_QC_thresholds":
        "NO",

    "treatment_used_for_clustering":
        "NO",

    "treatment_used_for_annotation":
        "NO",

    "FindAllMarkers_used":
        "NO",

    "KRT5_KRT14_treated_as_zero":
        "NO_STRUCTURAL_NA",

    "histology_used_for_cell_annotation":
        "NO"
}


for key, expected in (
    anti_expected.items()
):

    observed = anti.get(
        key,
        "<missing>"
    )

    add_check(
        "anti_" + key,
        observed == expected,
        observed,
        expected
    )


print(
    "Anti-circularity checks loaded:",
    len(
        anti_expected
    )
)

print()


# ============================================================
# 9. Pseudobulk provenance firewall
# ============================================================

print("===== 9. PSEUDOBULK / CORE FIREWALL =====")
print()

pb = dict_from_rows(
    read_csv(
        PB_PROVENANCE
    ),
    "item",
    "value"
)


pb_expected = {
    "technical_PASS_GE20":
        "26",

    "technical_FAIL_LT20":
        "1",

    "complete_paired_patients":
        "6",

    "technically_evaluable_paired_patients":
        "5",

    "CORE_rows_retained_in_raw_structure":
        "YES_4_OF_4",

    "CORE_gene_specific_values_printed":
        "NO",

    "normalization_performed":
        "NO",

    "CPM_calculated":
        "NO",

    "log2CPM_calculated":
        "NO",

    "zscore_calculated":
        "NO",

    "CORE2_calculated":
        "NO",

    "CORE4_calculated":
        "NO",

    "Before_After_inference":
        "NO",

    "histology_mapping":
        "INCOMPLETE",

    "EASC_identity":
        "UNRESOLVED",

    "primary_inference_allowed":
        "NO"
}


for key, expected in (
    pb_expected.items()
):

    observed = pb.get(
        key,
        "<missing>"
    )

    print(
        f"{key:<40s}: "
        f"{observed}"
    )

    add_check(
        "pb_" + key,
        observed == expected,
        observed,
        expected
    )


print()


# ============================================================
# 10. Lock manifest audit
# ============================================================

print("===== 10. STEP22A-6E MEMBERSHIP MANIFEST =====")
print()

manifest = dict_from_rows(
    read_csv(
        TECH_MANIFEST
    ),
    "item",
    "value"
)


manifest_expected = {
    "technical_threshold":
        "FIBROBLAST_CELLS_GE20",

    "technical_reference_candidate_n":
        "26",

    "technical_fail_n":
        "1",

    "technical_fail_fibroblast_cells":
        "18",

    "complete_paired_patients_n":
        "6",

    "technical_paired_candidate_n":
        "5",

    "technical_paired_candidate_ids":
        "P6;P9;P16;P18;P19",

    "P15_paired_status":
        "EXCLUDED_BEFORE_LT20",

    "threshold_relaxed":
        "NO",

    "borderline_rescue_branch_created":
        "NO",

    "raw_P15_before_retained":
        "YES_PROVENANCE_ONLY",

    "histology_scenarios":
        "23_TOTAL_1_UNRESOLVED_PLUS_22_HYPOTHETICAL",

    "histology_mapping":
        "INCOMPLETE",

    "normalization_executed":
        "NO",

    "CORE_gene_values_read":
        "NO",

    "CORE2_scored":
        "NO",

    "CORE4_scored":
        "NO",

    "primary_inference_allowed":
        "NO",

    "lock_status":
        "FINAL_TECHNICAL_REFERENCE_FROZEN"
}


for key, expected in (
    manifest_expected.items()
):

    observed = manifest.get(
        key,
        "<missing>"
    )

    add_check(
        "manifest_" + key,
        observed == expected,
        observed,
        expected
    )


print(
    "Manifest checks:",
    len(
        manifest_expected
    )
)

print()


# ============================================================
# 11. Frozen original project integrity
# ============================================================

print("===== 11. FROZEN STEP19-21 INTEGRITY =====")
print()

frozen = read_csv(
    FROZEN_INV
)

changes = []


for row in frozen:

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
            str(path)
            + ":MISSING"
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

            if expected != actual:

                changes.append(
                    str(path)
                    + ":SIZE_CHANGED"
                )

        except Exception:
            pass


print(
    "Frozen original files changed:",
    len(
        changes
    )
)


add_check(
    "frozen_Step19_21_integrity",
    len(
        changes
    )
    == 0,
    len(
        changes
    ),
    0
)


print()


# ============================================================
# 12. Final contradiction audit
# ============================================================

print("===== 12. CROSS-LOCK CONTRADICTION AUDIT =====")
print()


hard_fails = [
    x
    for x in checks
    if (
        x["status"] == "FAIL"
        and
        x["severity"] == "HARD"
    )
]


print(
    "Total governance checks:",
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
        hard_fails
    )
)

print(
    "FAIL:",
    len(
        hard_fails
    )
)

print()


write_csv(
    AUDIT_OUT,
    checks,
    [
        "check",
        "status",
        "observed",
        "expected",
        "severity"
    ]
)


if hard_fails:

    print(
        "FAILED CHECKS:"
    )

    for x in hard_fails:

        print(
            " ",
            x[
                "check"
            ]
        )

        print(
            "    observed:",
            x[
                "observed"
            ]
        )

        print(
            "    expected:",
            x[
                "expected"
            ]
        )

    print()

    fail(
        "STEP22A_GOVERNANCE_CONTRADICTION_FOUND"
    )


print(
    "✓ No contradiction detected across frozen Step22A locks."
)

print()


# ============================================================
# 13. Scientific decision dossier
# ============================================================

technical_sha = manifest.get(
    "technical_reference_membership_SHA256",
    "<missing>"
)

paired_sha = manifest.get(
    "technical_paired_membership_SHA256",
    "<missing>"
)

scoring_lock_sha = sha256_file(
    SCORING_LOCK
)


with open(
    DECISION_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-6F FINAL SCIENTIFIC DECISION DOSSIER\n"
    )

    f.write(
        "=" * 72
        + "\n\n"
    )

    f.write(
        "1. TECHNICAL PIPELINE\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Status: FROZEN AND INTERNALLY CONSISTENT\n"
    )

    f.write(
        "Tumor patient-timepoints: 27\n"
    )

    f.write(
        "Technical reference candidates: 26\n"
    )

    f.write(
        "Technical failure: P15 Before, 18 fibroblast cells\n"
    )

    f.write(
        "P15 Before raw data retained: YES\n"
    )

    f.write(
        "P15 Before technical reference: NO\n"
    )

    f.write(
        "Threshold relaxation: NO\n"
    )

    f.write(
        "Borderline rescue branch: NO\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "2. PAIRED TECHNICAL SET\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Complete paired patients: 6\n"
    )

    f.write(
        "Technically evaluable candidates: 5\n"
    )

    f.write(
        "P6, P9, P16, P18, P19\n"
    )

    f.write(
        "P15 excluded because Before failed frozen >=20 rule.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "3. HISTOLOGY\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "EASC identity: UNRESOLVED\n"
    )

    f.write(
        "Primary ESCC-only reference: NOT YET IDENTIFIABLE\n"
    )

    f.write(
        "Pre-frozen histology scenarios: 23\n"
    )

    f.write(
        "  1 unresolved baseline\n"
    )

    f.write(
        "  22 hypothetical EASC identities\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "4. SCORING FIREWALL\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Frozen transform: raw counts -> CPM -> log2(CPM+1)\n"
    )

    f.write(
        "Frozen z-score framework: within eligible OMIX tumor reference\n"
    )

    f.write(
        "Prior-cohort mean/SD reuse: FORBIDDEN\n"
    )

    f.write(
        "Paired-only z-score reference: FORBIDDEN\n"
    )

    f.write(
        "CORE gene values read: NO\n"
    )

    f.write(
        "Normalization executed: NO\n"
    )

    f.write(
        "CORE2 calculated: NO\n"
    )

    f.write(
        "CORE4 calculated: NO\n"
    )

    f.write(
        "Paired inference: NO\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "5. CURRENT SCIENTIFIC INTERPRETATION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "The computational pipeline is no longer the limiting factor.\n"
    )

    f.write(
        "The remaining blocker is the unresolved identity of the one\n"
    )

    f.write(
        "EASC patient in a 22-patient source cohort.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Primary inference should NOT be redefined as a hypothetical\n"
    )

    f.write(
        "histology scenario after seeing this limitation.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Recommended hierarchy:\n"
    )

    f.write(
        "  PRIMARY: remain ESCC-only and BLOCKED until histology resolves.\n"
    )

    f.write(
        "  SECONDARY: pre-frozen 22-scenario histology sensitivity only.\n"
    )

    f.write(
        "  FALLBACK: report technical feasibility / methodological limit\n"
    )

    f.write(
        "            without claiming external CORE2 inference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "6. RECOMMENDED NEXT ACTION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "First attempt to resolve EASC identity from the original authors\n"
    )

    f.write(
        "or an authoritative patient-level pathology mapping.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Do not open the CORE2 scoring gate merely because the public\n"
    )

    f.write(
        "metadata are insufficient.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "If histology cannot be resolved after documented attempts,\n"
    )

    f.write(
        "retain the ESCC-only primary analysis as unavailable and decide\n"
    )

    f.write(
        "separately whether to run the already pre-specified scenarios\n"
    )

    f.write(
        "as explicitly secondary sensitivity analyses.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "7. MEMBERSHIP HASHES\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Technical reference SHA256:\n"
        + technical_sha
        + "\n\n"
    )

    f.write(
        "Paired technical set SHA256:\n"
        + paired_sha
        + "\n\n"
    )

    f.write(
        "Frozen scoring-lock file SHA256:\n"
        + scoring_lock_sha
        + "\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "FINAL GOVERNANCE STATUS:\n"
    )

    f.write(
        "PASS_NO_CONTRADICTIONS_PRIMARY_HISTOLOGY_BLOCKED\n"
    )


# ============================================================
# FINAL
# ============================================================

print("=" * 80)
print("STEP22A-6F COMPLETE")
print("=" * 80)
print()

print(
    "GOVERNANCE CHECKS             :",
    len(
        checks
    )
)

print(
    "GOVERNANCE FAILURES           : 0"
)

print()

print(
    "TECHNICAL REFERENCE           : 26 / 27"
)

print(
    "TECHNICAL FAILURE             : P15 Before = 18 cells"
)

print(
    "THRESHOLD RELAXED             : NO"
)

print(
    "BORDERLINE RESCUE             : NO"
)

print()

print(
    "PAIRED TECHNICAL CANDIDATES   : 5"
)

print(
    "PAIRED IDS                    : P6, P9, P16, P18, P19"
)

print()

print(
    "HISTOLOGY SCENARIOS           : 23"
)

print(
    "EASC IDENTITY                 : UNRESOLVED"
)

print(
    "PRIMARY ESCC REFERENCE        : BLOCKED"
)

print()

print(
    "RAW PSEUDOBULK RDS OPENED     : NO"
)

print(
    "CORE GENE VALUES READ         : NO"
)

print(
    "NORMALIZATION                 : NOT PERFORMED"
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
    "SCIENTIFIC DECISION:"
)

print(
    "  PRIMARY remains ESCC-only and BLOCKED."
)

print(
    "  Do NOT promote a hypothetical EASC scenario to primary."
)

print(
    "  Next preferred action = resolve histology externally."
)

print()

print(
    "FINAL GOVERNANCE STATUS:"
)

print(
    "PASS_NO_CONTRADICTIONS_PRIMARY_HISTOLOGY_BLOCKED"
)

print()

print(
    "DECISION DOSSIER:"
)

print(
    "  00_metadata/OMIX005710/"
)

print(
    "  STEP22A_06F_FINAL_DECISION_DOSSIER.txt"
)

print()
