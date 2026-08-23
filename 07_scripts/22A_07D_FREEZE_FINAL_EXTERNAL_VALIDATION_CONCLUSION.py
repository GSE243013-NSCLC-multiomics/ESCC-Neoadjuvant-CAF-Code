#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-7D
#
# FREEZE FINAL OMIX005710 EXTERNAL-VALIDATION CONCLUSION
#
# This is a GOVERNANCE / CONCLUSION step only.
#
# It does NOT:
#   - open raw pseudobulk RDS
#   - read expression matrices
#   - recompute CPM
#   - recompute z-scores
#   - recompute CORE2 / CORE4
#   - perform any new statistical test
#   - perform subtype discovery
#   - perform post-hoc biological explanation
#
# Final scientific hierarchy:
#
# PRIMARY:
#   ESCC-only validation unavailable because the identity
#   of the single EASC patient remains unresolved.
#
# SECONDARY:
#   Pre-frozen 22-scenario histology sensitivity shows
#   robust direction opposite to the discovery direction.
#
# This is NON-REPLICATION / DIRECTIONAL DISCORDANCE,
# not a "failed experiment" and not evidence that the
# original biology is necessarily false.
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

TABLE = (
    PROJECT
    / "06_tables"
    / "external_validation"
    / "Step22A"
)

RESULT = (
    PROJECT
    / "04_results"
    / "external_validation"
    / "Step22A"
)

LOG = (
    PROJECT
    / "08_logs"
)

RESULT.mkdir(
    parents=True,
    exist_ok=True
)

FINAL_LOCK_OUT = (
    META
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_LOCK.csv"
)

FINAL_DOSSIER_OUT = (
    META
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_DOSSIER.txt"
)

FINAL_SUMMARY_OUT = (
    RESULT
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_SUMMARY.txt"
)

MANUSCRIPT_TEXT_OUT = (
    RESULT
    / "STEP22A_07D_REPORTING_LANGUAGE.txt"
)

GATE_FILE = (
    META
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

GOVERNANCE_6F = (
    META
    / "STEP22A_06F_FINAL_GOVERNANCE_AUDIT.csv"
)

AUTH_7A = (
    META
    / "STEP22A_07A_SECONDARY_SENSITIVITY_AUTHORIZATION.csv"
)

RESULTS_7A = (
    TABLE
    / "STEP22A_07A_SECONDARY_PAIRED_SCENARIO_RESULTS.csv"
)

ROBUSTNESS_7A = (
    TABLE
    / "STEP22A_07A_SECONDARY_HISTOLOGY_ROBUSTNESS_SUMMARY.csv"
)

AUDIT_7B = (
    META
    / "STEP22A_07B_CORE2_TECHNICAL_REPLICATION_AUDIT.csv"
)

SCENARIO_7C = (
    TABLE
    / "STEP22A_07C_SCENARIO_COMPARATOR_SUMMARY.csv"
)

LOPO_7C = (
    TABLE
    / "STEP22A_07C_PAIRED_LOPO_RESULTS.csv"
)

ROBUSTNESS_7C = (
    TABLE
    / "STEP22A_07C_GLOBAL_ROBUSTNESS_SUMMARY.csv"
)

FROZEN_INV = (
    META
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)


print()
print("=" * 80)
print("STEP22A-7D: FREEZE FINAL OMIX005710 CONCLUSION")
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


def kv(
    rows,
    key_col,
    value_col
):

    return {
        clean(
            x.get(
                key_col
            )
        ):
        clean(
            x.get(
                value_col
            )
        )
        for x in rows
    }


def as_float(x):

    try:
        return float(
            clean(x)
        )

    except Exception:
        return float("nan")


def as_int(x):

    return int(
        float(
            clean(x)
        )
    )


def file_sha256(path):

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


# ============================================================
# 1. Verify primary gate remains closed
# ============================================================

print("===== 1. GOVERNANCE STATUS =====")
print()

gates = kv(
    read_csv(
        GATE_FILE
    ),
    "gate",
    "allowed"
)

primary_gate = gates.get(
    "paired_CORE2_test",
    "<missing>"
)

manuscript_gate = gates.get(
    "manuscript_update",
    "<missing>"
)

print(
    "Primary CORE2 gate :",
    primary_gate
)

print(
    "Manuscript gate    :",
    manuscript_gate
)

print()


if primary_gate == "YES":

    fail(
        "Primary CORE2 gate unexpectedly open."
    )


if manuscript_gate == "YES":

    fail(
        "Manuscript gate unexpectedly open."
    )


# ============================================================
# 2. Verify 6F governance
# ============================================================

print("===== 2. STEP22A-6F GOVERNANCE =====")
print()

gov = read_csv(
    GOVERNANCE_6F
)

gov_fail = [
    x
    for x in gov
    if clean(
        x.get(
            "status"
        )
    )
    != "PASS"
]

print(
    "Governance checks:",
    len(
        gov
    )
)

print(
    "Failures:",
    len(
        gov_fail
    )
)

print()


if gov_fail:

    fail(
        "Step22A-6F governance no longer all PASS."
    )


# ============================================================
# 3. Verify secondary authorization
# ============================================================

print("===== 3. SECONDARY-BRANCH AUTHORIZATION =====")
print()

auth = kv(
    read_csv(
        AUTH_7A
    ),
    "item",
    "status"
)

secondary_auth = auth.get(
    "SECONDARY_HISTOLOGY_SENSITIVITY",
    "<missing>"
)

scenario_selection = auth.get(
    "SCENARIO_SELECTION_ALLOWED",
    "<missing>"
)

all_scenarios_required = auth.get(
    "ALL_PRE_FROZEN_EASC_SCENARIOS_REQUIRED",
    "<missing>"
)

print(
    "Secondary histology sensitivity :",
    secondary_auth
)

print(
    "Scenario selection allowed      :",
    scenario_selection
)

print(
    "All frozen scenarios required   :",
    all_scenarios_required
)

print()


if secondary_auth != "AUTHORIZED":

    fail(
        "Secondary branch not authorized."
    )


if scenario_selection != "NO":

    fail(
        "Scenario-selection firewall broken."
    )


if all_scenarios_required != "YES_22_OF_22":

    fail(
        "22-scenario completeness lock missing."
    )


# ============================================================
# 4. Verify 7A CORE2 directional result
# ============================================================

print("===== 4. STEP22A-7A CORE2 DIRECTION =====")
print()

results_7a = read_csv(
    RESULTS_7A
)

hyp_core2 = [
    x
    for x in results_7a
    if (
        clean(
            x.get(
                "scenario_role"
            )
        )
        ==
        "SECONDARY_HISTOLOGY_SENSITIVITY"
        and
        clean(
            x.get(
                "endpoint"
            )
        )
        ==
        "CORE2_STRONG"
    )
]


if len(
    hyp_core2
) != 22:

    fail(
        "Expected exactly 22 hypothetical-EASC CORE2 results."
    )


mean_negative_7a = sum(
    as_float(
        x.get(
            "mean_delta"
        )
    )
    < 0
    for x in hyp_core2
)

median_negative_7a = sum(
    as_float(
        x.get(
            "median_delta"
        )
    )
    < 0
    for x in hyp_core2
)

majority_negative_7a = sum(
    as_int(
        x.get(
            "negative_n"
        )
    )
    >
    (
        as_int(
            x.get(
                "n_paired"
            )
        )
        / 2
    )
    for x in hyp_core2
)


mean_delta_min = min(
    as_float(
        x.get(
            "mean_delta"
        )
    )
    for x in hyp_core2
)

mean_delta_max = max(
    as_float(
        x.get(
            "mean_delta"
        )
    )
    for x in hyp_core2
)

median_delta_min = min(
    as_float(
        x.get(
            "median_delta"
        )
    )
    for x in hyp_core2
)

median_delta_max = max(
    as_float(
        x.get(
            "median_delta"
        )
    )
    for x in hyp_core2
)


print(
    "CORE2 mean delta negative   :",
    mean_negative_7a,
    "/22"
)

print(
    "CORE2 median delta negative :",
    median_negative_7a,
    "/22"
)

print(
    "Majority patients negative  :",
    majority_negative_7a,
    "/22"
)

print(
    "Mean delta range            :",
    mean_delta_min,
    "to",
    mean_delta_max
)

print(
    "Median delta range          :",
    median_delta_min,
    "to",
    median_delta_max
)

print()


if not (
    mean_negative_7a == 22
    and
    median_negative_7a == 22
    and
    majority_negative_7a == 22
):

    fail(
        "Step22A-7A is not uniformly negative as expected."
    )


# ============================================================
# 5. Verify Step22A-7B technical replication
# ============================================================

print("===== 5. STEP22A-7B TECHNICAL REPLICATION =====")
print()

audit_7b = read_csv(
    AUDIT_7B
)

scenario_numeric_rows = [
    x
    for x in audit_7b
    if clean(
        x.get(
            "scenario_id"
        )
    )
]


if len(
    scenario_numeric_rows
) < 22:

    fail(
        "Step22A-7B scenario audit incomplete."
    )


max_core2_score_diff_7b = max(
    as_float(
        x.get(
            "CORE2_score_max_abs_diff"
        )
    )
    for x in scenario_numeric_rows
)

p15_reference_violations = sum(
    clean(
        x.get(
            "P15_Before_in_reference"
        )
    ).upper()
    in {
        "TRUE",
        "T",
        "YES",
        "1"
    }
    for x in scenario_numeric_rows
)


print(
    "CORE2 max abs score difference:",
    max_core2_score_diff_7b
)

print(
    "P15 Before reference violations:",
    p15_reference_violations
)

print()


if max_core2_score_diff_7b > 1e-10:

    fail(
        "7B independent replication exceeds tolerance."
    )


if p15_reference_violations != 0:

    fail(
        "P15 Before entered a technical reference."
    )


# ============================================================
# 6. Verify Step22A-7C comparator / component pattern
# ============================================================

print("===== 6. STEP22A-7C MODULE PATTERN =====")
print()

sc7c = read_csv(
    SCENARIO_7C
)

if len(
    sc7c
) != 22:

    fail(
        "Expected exactly 22 Step22A-7C scenario rows."
    )


core2_mean_neg_7c = sum(
    as_float(
        x.get(
            "CORE2_mean_delta"
        )
    )
    < 0
    for x in sc7c
)

core4_mean_neg_7c = sum(
    as_float(
        x.get(
            "CORE4_mean_delta"
        )
    )
    < 0
    for x in sc7c
)

weak_mean_neg_7c = sum(
    as_float(
        x.get(
            "BGN_TIMP1_mean_delta"
        )
    )
    < 0
    for x in sc7c
)

weak_mean_pos_7c = sum(
    as_float(
        x.get(
            "BGN_TIMP1_mean_delta"
        )
    )
    > 0
    for x in sc7c
)

cdkn1b_contrib_neg = sum(
    as_float(
        x.get(
            "CDKN1B_mean_CORE2_contribution"
        )
    )
    < 0
    for x in sc7c
)

gsn_contrib_neg = sum(
    as_float(
        x.get(
            "GSN_mean_CORE2_contribution"
        )
    )
    < 0
    for x in sc7c
)


print(
    "CORE2 mean negative     :",
    core2_mean_neg_7c,
    "/22"
)

print(
    "CORE4 mean negative     :",
    core4_mean_neg_7c,
    "/22"
)

print(
    "BGN/TIMP1 mean negative :",
    weak_mean_neg_7c,
    "/22"
)

print(
    "BGN/TIMP1 mean positive :",
    weak_mean_pos_7c,
    "/22"
)

print(
    "CDKN1B contribution neg.:",
    cdkn1b_contrib_neg,
    "/22"
)

print(
    "GSN contribution neg.   :",
    gsn_contrib_neg,
    "/22"
)

print()


if core2_mean_neg_7c != 22:

    fail(
        "CORE2 Step22A-7C direction changed."
    )


if core4_mean_neg_7c != 22:

    fail(
        "CORE4 is not negative in all 22 scenarios."
    )


if cdkn1b_contrib_neg != 22:

    fail(
        "CDKN1B contribution not negative in all 22 scenarios."
    )


if gsn_contrib_neg != 22:

    fail(
        "GSN contribution not negative in all 22 scenarios."
    )


# ============================================================
# 7. Verify 105 LOPO runs
# ============================================================

print("===== 7. STEP22A-7C LOPO =====")
print()

lopo = read_csv(
    LOPO_7C
)

print(
    "LOPO runs:",
    len(
        lopo
    )
)

if len(
    lopo
) != 105:

    fail(
        "Expected exactly 105 LOPO runs."
    )


lopo_mean_neg = sum(
    as_float(
        x.get(
            "CORE2_mean_delta"
        )
    )
    < 0
    for x in lopo
)

lopo_median_neg = sum(
    as_float(
        x.get(
            "CORE2_median_delta"
        )
    )
    < 0
    for x in lopo
)

lopo_majority_neg = sum(
    as_int(
        x.get(
            "CORE2_negative_n"
        )
    )
    >
    (
        as_int(
            x.get(
                "remaining_pair_n"
            )
        )
        / 2
    )
    for x in lopo
)

lopo_all_neg = sum(
    clean(
        x.get(
            "CORE2_all_remaining_negative"
        )
    ).upper()
    in {
        "TRUE",
        "T",
        "YES",
        "1"
    }
    for x in lopo
)


print(
    "LOPO mean negative       :",
    lopo_mean_neg,
    "/105"
)

print(
    "LOPO median negative     :",
    lopo_median_neg,
    "/105"
)

print(
    "LOPO majority negative   :",
    lopo_majority_neg,
    "/105"
)

print(
    "All remaining pairs neg. :",
    lopo_all_neg,
    "/105"
)

print()


if not (
    lopo_mean_neg == 105
    and
    lopo_median_neg == 105
    and
    lopo_majority_neg == 105
    and
    lopo_all_neg == 105
):

    fail(
        "LOPO negative-direction lock not satisfied."
    )


# ============================================================
# 8. Verify robustness classifications
# ============================================================

print("===== 8. FROZEN CLASSIFICATIONS =====")
print()

rob7c = kv(
    read_csv(
        ROBUSTNESS_7C
    ),
    "metric",
    "value"
)


expected_lopo_status = (
    "CORE2_NEGATIVE_ROBUST_TO_ALL_PAIRED_LOPO_RUNS"
)

expected_component = (
    "CDKN1B_AND_GSN_BOTH_CONTRIBUTE_NEGATIVELY_IN_ALL_22_SCENARIOS"
)

expected_module = (
    "CORE4_NEGATIVE_ALL_22_WITH_MIXED_BGN_TIMP1_ARM"
)


observed_lopo_status = rob7c.get(
    "CORE2_LOPO_status",
    "<missing>"
)

observed_component = rob7c.get(
    "CORE2_component_pattern",
    "<missing>"
)

observed_module = rob7c.get(
    "four_gene_module_pattern",
    "<missing>"
)


print(
    "LOPO status      :",
    observed_lopo_status
)

print(
    "Component pattern:",
    observed_component
)

print(
    "Module pattern   :",
    observed_module
)

print()


if observed_lopo_status != expected_lopo_status:

    fail(
        "Unexpected LOPO classification."
    )


if observed_component != expected_component:

    fail(
        "Unexpected CORE2 component classification."
    )


if observed_module != expected_module:

    fail(
        "Unexpected four-gene module classification."
    )


# ============================================================
# 9. Frozen project integrity
# ============================================================

print("===== 9. FROZEN STEP19-21 INTEGRITY =====")
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

            expected_size = int(
                float(
                    expected
                )
            )

            actual_size = (
                path.stat().st_size
            )

            if (
                expected_size
                != actual_size
            ):

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

print()


if changes:

    for x in changes[:20]:
        print(x)

    fail(
        "FROZEN_STEP19_21_CHANGED"
    )


# ============================================================
# 10. Final scientific classification
# ============================================================

print("===== 10. FINAL SCIENTIFIC CLASSIFICATION =====")
print()


primary_status = (
    "PRIMARY_ESCC_ONLY_VALIDATION_UNAVAILABLE_HISTOLOGY_UNRESOLVED"
)

secondary_status = (
    "SECONDARY_HISTOLOGY_SENSITIVITY_DIRECTION_ROBUSTLY_OPPOSITE"
)

technical_status = (
    "TECHNICAL_AUDITS_PASS_NO_SIGN_PAIR_REFERENCE_OR_SINGLE_PATIENT_ARTIFACT"
)

core2_status = (
    "CORE2_STRONG_DIRECTION_NOT_REPLICATED_IN_OMIX005710_SECONDARY_ANALYSIS"
)

core4_status = (
    "CORE4_DIRECTION_ALSO_OPPOSITE_ACROSS_ALL_22_HISTOLOGY_SCENARIOS"
)

component_status = (
    "CDKN1B_AND_GSN_BOTH_SUPPORT_OPPOSITE_CORE2_DIRECTION"
)

weak_arm_status = (
    "BGN_TIMP1_COMPONENT_MIXED_18_NEGATIVE_4_POSITIVE_SCENARIOS"
)


print(
    "Primary status  :",
    primary_status
)

print(
    "Secondary status:",
    secondary_status
)

print(
    "Technical status:",
    technical_status
)

print()


# ============================================================
# 11. Freeze lock
# ============================================================

lock_rows = [
    {
        "item":
            "primary_external_validation",

        "value":
            primary_status
    },

    {
        "item":
            "secondary_histology_sensitivity",

        "value":
            secondary_status
    },

    {
        "item":
            "technical_audit_status",

        "value":
            technical_status
    },

    {
        "item":
            "CORE2_STRONG_status",

        "value":
            core2_status
    },

    {
        "item":
            "CORE4_status",

        "value":
            core4_status
    },

    {
        "item":
            "CORE2_component_status",

        "value":
            component_status
    },

    {
        "item":
            "BGN_TIMP1_component_status",

        "value":
            weak_arm_status
    },

    {
        "item":
            "hypothetical_EASC_scenarios_completed",

        "value":
            "22_OF_22"
    },

    {
        "item":
            "scenario_selection",

        "value":
            "NONE"
    },

    {
        "item":
            "CORE2_mean_negative_scenarios",

        "value":
            "22_OF_22"
    },

    {
        "item":
            "CORE2_median_negative_scenarios",

        "value":
            "22_OF_22"
    },

    {
        "item":
            "CORE4_mean_negative_scenarios",

        "value":
            "22_OF_22"
    },

    {
        "item":
            "CDKN1B_negative_contribution_scenarios",

        "value":
            "22_OF_22"
    },

    {
        "item":
            "GSN_negative_contribution_scenarios",

        "value":
            "22_OF_22"
    },

    {
        "item":
            "BGN_TIMP1_negative_scenarios",

        "value":
            "18_OF_22"
    },

    {
        "item":
            "BGN_TIMP1_positive_scenarios",

        "value":
            "4_OF_22"
    },

    {
        "item":
            "paired_LOPO_runs",

        "value":
            "105"
    },

    {
        "item":
            "paired_LOPO_mean_negative",

        "value":
            "105_OF_105"
    },

    {
        "item":
            "paired_LOPO_median_negative",

        "value":
            "105_OF_105"
    },

    {
        "item":
            "paired_LOPO_majority_negative",

        "value":
            "105_OF_105"
    },

    {
        "item":
            "paired_LOPO_all_remaining_pairs_negative",

        "value":
            "105_OF_105"
    },

    {
        "item":
            "posthoc_CAF_subtype_explanation_performed",

        "value":
            "NO"
    },

    {
        "item":
            "CORE4_used_as_rescue_endpoint",

        "value":
            "NO"
    },

    {
        "item":
            "manuscript_modified",

        "value":
            "NO"
    },

    {
        "item":
            "final_lock_status",

        "value":
            "STEP22A_OMIX005710_EXTERNAL_VALIDATION_CONCLUSION_FROZEN"
    }
]


write_csv(
    FINAL_LOCK_OUT,
    lock_rows,
    [
        "item",
        "value"
    ]
)


# ============================================================
# 12. Human-readable dossier
# ============================================================

with open(
    FINAL_DOSSIER_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-7D FINAL OMIX005710 EXTERNAL VALIDATION DOSSIER\n"
    )

    f.write(
        "=" * 76
        + "\n\n"
    )

    f.write(
        "PRIMARY ANALYSIS\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Status: UNAVAILABLE\n"
    )

    f.write(
        "Reason: the single EASC patient cannot be mapped to an\n"
    )

    f.write(
        "individual patient using the available public metadata.\n"
    )

    f.write(
        "Therefore a definitive ESCC-only primary reference cannot\n"
    )

    f.write(
        "be constructed without making an unsupported assumption.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "SECONDARY PRE-SPECIFIED HISTOLOGY SENSITIVITY\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "All 22 hypothetical EASC identities were evaluated.\n"
    )

    f.write(
        "No scenario was selected or promoted to primary.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CORE2_STRONG:\n"
    )

    f.write(
        "  mean After-Before delta negative: 22/22 scenarios\n"
    )

    f.write(
        "  median After-Before delta negative: 22/22 scenarios\n"
    )

    f.write(
        "  majority paired patients negative: 22/22 scenarios\n"
    )

    f.write(
        f"  mean delta range: {mean_delta_min:.6f} to {mean_delta_max:.6f}\n"
    )

    f.write(
        f"  median delta range: {median_delta_min:.6f} to {median_delta_max:.6f}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CORE4:\n"
    )

    f.write(
        "  mean After-Before delta negative: 22/22 scenarios\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CORE2 COMPONENTS:\n"
    )

    f.write(
        "  CDKN1B contribution negative: 22/22 scenarios\n"
    )

    f.write(
        "  GSN contribution negative: 22/22 scenarios\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "BGN/TIMP1 ARM:\n"
    )

    f.write(
        "  negative: 18/22 scenarios\n"
    )

    f.write(
        "  positive: 4/22 scenarios\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "PAIRED LOPO INFLUENCE AUDIT:\n"
    )

    f.write(
        "  total runs: 105\n"
    )

    f.write(
        "  mean delta negative: 105/105\n"
    )

    f.write(
        "  median delta negative: 105/105\n"
    )

    f.write(
        "  majority negative: 105/105\n"
    )

    f.write(
        "  all remaining paired patients negative: 105/105\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "TECHNICAL AUDIT\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Independent recalculation reproduced Step22A-7A within\n"
    )

    f.write(
        "floating-point precision. No sign, Before/After mapping,\n"
    )

    f.write(
        "reference-membership, P15 exclusion, or single-patient\n"
    )

    f.write(
        "influence artifact was detected.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "INTERPRETATION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "The pre-specified secondary analysis does not reproduce the\n"
    )

    f.write(
        "frozen discovery-direction increase in CORE2_STRONG.\n"
    )

    f.write(
        "Instead, the treatment-associated direction is consistently\n"
    )

    f.write(
        "opposite across all 22 histology scenarios.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CORE4 is also negative across all 22 scenarios, and both\n"
    )

    f.write(
        "CDKN1B and GSN components contribute in the same negative\n"
    )

    f.write(
        "CORE2 direction. The BGN/TIMP1 arm is mildly mixed.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "This should be reported as directional non-replication /\n"
    )

    f.write(
        "external discordance in the secondary sensitivity analysis.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "It should NOT be described as definitive ESCC-only external\n"
    )

    f.write(
        "validation because the EASC patient identity remains unresolved.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "No post-hoc CAF-subtype analysis is required to establish\n"
    )

    f.write(
        "this validation conclusion. Such analyses, if pursued later,\n"
    )

    f.write(
        "must be explicitly labeled exploratory and hypothesis-generating.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "FINAL STATUS:\n"
    )

    f.write(
        "PRIMARY_ESCC_ONLY_UNAVAILABLE__SECONDARY_DIRECTION_ROBUSTLY_OPPOSITE\n"
    )


# ============================================================
# 13. Reporting language
# ============================================================

with open(
    MANUSCRIPT_TEXT_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "SUGGESTED REPORTING LANGUAGE — OMIX005710\n"
    )

    f.write(
        "=" * 72
        + "\n\n"
    )

    f.write(
        "RESULTS VERSION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "In OMIX005710, definitive ESCC-only external validation was\n"
    )

    f.write(
        "not possible because the single patient with esophageal\n"
    )

    f.write(
        "adenosquamous carcinoma could not be identified from the\n"
    )

    f.write(
        "publicly available patient-level metadata. We therefore used\n"
    )

    f.write(
        "a pre-specified secondary sensitivity framework in which each\n"
    )

    f.write(
        "of the 22 patients was, in turn, assumed to represent the EASC\n"
    )

    f.write(
        "case and excluded from the ESCC reference.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Across all 22 histology scenarios, the paired treatment-related\n"
    )

    f.write(
        "change in CORE2_STRONG was consistently negative, opposite to\n"
    )

    f.write(
        "the frozen replication direction from the discovery cohorts.\n"
    )

    f.write(
        "This directional discordance was unchanged in leave-one-paired-\n"
    )

    f.write(
        "patient-out analyses. CORE4 was likewise negative across all\n"
    )

    f.write(
        "22 scenarios, while both CDKN1B and GSN contributed consistently\n"
    )

    f.write(
        "to the negative CORE2_STRONG direction. Thus, this secondary\n"
    )

    f.write(
        "analysis did not replicate the direction of the previously\n"
    )

    f.write(
        "observed CORE response in OMIX005710.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "DISCUSSION VERSION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "The OMIX005710 analysis provides an important negative external\n"
    )

    f.write(
        "finding. Although a definitive ESCC-only primary validation was\n"
    )

    f.write(
        "precluded by unresolved patient-level histology, the pre-specified\n"
    )

    f.write(
        "histology-sensitivity analysis yielded a highly consistent\n"
    )

    f.write(
        "direction opposite to that observed in the discovery cohorts.\n"
    )

    f.write(
        "The result was robust to all hypothetical EASC assignments,\n"
    )

    f.write(
        "independent score reconstruction, and paired leave-one-out\n"
    )

    f.write(
        "analysis, arguing against a simple computational, reference-set,\n"
    )

    f.write(
        "or single-patient explanation. These findings therefore indicate\n"
    )

    f.write(
        "that the treatment-associated CORE response is not uniformly\n"
    )

    f.write(
        "replicated across independent cohorts and may be context-dependent.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CAUTION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Do not state that OMIX005710 proves the opposite biological\n"
    )

    f.write(
        "mechanism. The result establishes directional discordance in\n"
    )

    f.write(
        "the frozen secondary sensitivity analysis, not its biological cause.\n"
    )


# ============================================================
# 14. Compact final summary
# ============================================================

with open(
    FINAL_SUMMARY_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-7D FINAL SUMMARY\n"
    )

    f.write(
        "=" * 60
        + "\n\n"
    )

    f.write(
        "Primary ESCC-only validation: UNAVAILABLE\n"
    )

    f.write(
        "Reason: EASC patient identity unresolved\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Secondary frozen histology sensitivity:\n"
    )

    f.write(
        "CORE2 mean negative: 22/22\n"
    )

    f.write(
        "CORE2 median negative: 22/22\n"
    )

    f.write(
        "CORE4 mean negative: 22/22\n"
    )

    f.write(
        "CDKN1B contribution negative: 22/22\n"
    )

    f.write(
        "GSN contribution negative: 22/22\n"
    )

    f.write(
        "BGN/TIMP1 negative: 18/22\n"
    )

    f.write(
        "BGN/TIMP1 positive: 4/22\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "LOPO mean negative: 105/105\n"
    )

    f.write(
        "LOPO median negative: 105/105\n"
    )

    f.write(
        "LOPO majority negative: 105/105\n"
    )

    f.write(
        "LOPO all remaining patients negative: 105/105\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Technical artifact detected: NO\n"
    )

    f.write(
        "Single-patient driver detected: NO\n"
    )

    f.write(
        "Scenario selection performed: NO\n"
    )

    f.write(
        "CORE4 used as rescue endpoint: NO\n"
    )

    f.write(
        "Post-hoc biological explanation performed: NO\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "FINAL INTERPRETATION:\n"
    )

    f.write(
        "PRIMARY_ESCC_ONLY_UNAVAILABLE__SECONDARY_DIRECTION_ROBUSTLY_OPPOSITE\n"
    )


# ============================================================
# FINAL
# ============================================================

print("=" * 80)
print("STEP22A-7D COMPLETE")
print("=" * 80)
print()

print(
    "PRIMARY ESCC-ONLY VALIDATION   : UNAVAILABLE"
)

print(
    "PRIMARY BLOCKER                : EASC IDENTITY UNRESOLVED"
)

print()

print(
    "SECONDARY HISTOLOGY SCENARIOS : 22 / 22"
)

print(
    "CORE2 MEAN NEGATIVE            : 22 / 22"
)

print(
    "CORE2 MEDIAN NEGATIVE          : 22 / 22"
)

print(
    "CORE4 MEAN NEGATIVE            : 22 / 22"
)

print(
    "CDKN1B CONTRIB. NEGATIVE       : 22 / 22"
)

print(
    "GSN CONTRIB. NEGATIVE          : 22 / 22"
)

print(
    "BGN/TIMP1 NEGATIVE             : 18 / 22"
)

print(
    "BGN/TIMP1 POSITIVE             : 4 / 22"
)

print()

print(
    "PAIRED LOPO                    : 105 RUNS"
)

print(
    "LOPO MEAN NEGATIVE             : 105 / 105"
)

print(
    "LOPO MEDIAN NEGATIVE           : 105 / 105"
)

print(
    "LOPO MAJORITY NEGATIVE         : 105 / 105"
)

print(
    "LOPO ALL REMAINING NEGATIVE    : 105 / 105"
)

print()

print(
    "TECHNICAL ARTIFACT DETECTED    : NO"
)

print(
    "SINGLE-PATIENT DRIVER          : NO"
)

print(
    "SCENARIO SELECTION             : NO"
)

print(
    "CORE4 USED AS RESCUE           : NO"
)

print(
    "POST-HOC BIOLOGY PERFORMED     : NO"
)

print()

print(
    "FINAL SCIENTIFIC STATUS:"
)

print(
    "PRIMARY_ESCC_ONLY_UNAVAILABLE__SECONDARY_DIRECTION_ROBUSTLY_OPPOSITE"
)

print()

print(
    "INTERPRETATION:"
)

print(
    "OMIX005710 does not replicate the frozen treatment-associated"
)

print(
    "CORE direction in the pre-specified secondary sensitivity analysis."
)

print()

print(
    "Do NOT claim a definitive biological mechanism for the reversal."
)

print()

print(
    "STEP22A EXTERNAL-VALIDATION PIPELINE: CLOSED."
)

print()
