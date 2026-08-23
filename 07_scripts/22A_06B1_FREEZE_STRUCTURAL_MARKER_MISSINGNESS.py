#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-6B1
#
# FREEZE STRUCTURAL MARKER MISSINGNESS
#
# Decision:
#
# - The biological 33-marker dictionary remains FROZEN.
# - KRT5 / KRT14 are NOT removed from that dictionary.
# - Feature-table absence is treated as NOT MEASURED,
#   never as zero expression / negative evidence.
# - Joint all-27 annotation may use only markers measurable
#   in the all-27 common ENSG universe.
# - No sample is excluded because of KRT5/KRT14 absence.
# - CORE genes must remain 4/4 available.
#
# NO expression matrices are read.
# NO cell annotation is performed.
# NO CORE values are read.
# NO scoring / inference is performed.
# ============================================================

from pathlib import Path
import csv
import sys


PROJECT = Path.cwd()

META_DIR = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

FEATURE_TABLE_COUNTS = (
    META_DIR
    / "STEP22A_06B_ALL27_FEATURE_TABLE_COUNTS.csv"
)

ALL27_COMMON = (
    META_DIR
    / "STEP22A_06B_ALL27_COMMON_ENSG_FEATURES.csv"
)

COMBINED_INV = (
    META_DIR
    / "STEP22A_06A_COMBINED_27_TUMOR_TECHNICAL_INVENTORY.csv"
)

GATE_FILE = (
    META_DIR
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

BY_SAMPLE_OUT = (
    META_DIR
    / "STEP22A_06B1_MARKER_AVAILABILITY_BY_SAMPLE.csv"
)

JOINT_LOCK_OUT = (
    META_DIR
    / "STEP22A_06B1_JOINT_ANNOTATION_MARKER_LOCK.csv"
)

DECISION_OUT = (
    META_DIR
    / "STEP22A_06B1_STRUCTURAL_MISSINGNESS_DECISION.txt"
)

STATUS_OUT = (
    META_DIR
    / "STEP22A_06B1_ALIGNMENT_FINAL_STATUS.csv"
)


print()
print("=" * 80)
print("STEP22A-6B1: FREEZE STRUCTURAL MARKER MISSINGNESS")
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


def read_symbols(path):

    symbols = set()

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for line in f:

            line = line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) >= 2:
                symbol = clean(
                    fields[1]
                )
            else:
                symbol = clean(
                    fields[0]
                )

            if symbol:
                symbols.add(
                    symbol.upper()
                )

    return symbols


def frozen_integrity():

    rows = read_csv(
        FROZEN_INV
    )

    changes = []

    for row in rows:

        raw = clean(
            row.get(
                "path"
            )
        )

        if not raw:
            continue

        path = Path(raw)

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
                        (
                            str(path),
                            "SIZE_CHANGED"
                        )
                    )

            except Exception:
                pass

    return changes


# ============================================================
# 1. Gate check
# ============================================================

print("===== 1. ANALYSIS GATE =====")
print()

gate_rows = read_csv(
    GATE_FILE
)

gates = {
    clean(
        x.get("gate")
    ):
    clean(
        x.get("allowed")
    )
    for x in gate_rows
}


for g in [
    "technical_preprocessing",
    "paired_CORE2_test",
    "manuscript_update"
]:

    if g not in gates:
        fail(
            "Missing gate: "
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


if gates[
    "technical_preprocessing"
] != "YES":

    fail(
        "Technical preprocessing blocked."
    )


if gates[
    "paired_CORE2_test"
] == "YES":

    fail(
        "CORE2 inference gate unexpectedly open."
    )


print(
    "✓ Technical marker-availability audit allowed."
)

print(
    "✓ CORE2 inference remains blocked."
)

print()


# ============================================================
# 2. Frozen project check
# ============================================================

print("===== 2. FROZEN PROJECT CHECK =====")
print()

changes = frozen_integrity()

print(
    "Frozen original files changed:",
    len(changes)
)

if changes:

    for x in changes[:20]:
        print(x)

    fail(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )


print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# 3. Frozen marker dictionary
# ============================================================

print("===== 3. FROZEN 33-MARKER DICTIONARY =====")
print()

marker_modules = {

    "Fibroblast_CAF": [
        "COL1A1",
        "COL1A2",
        "DCN",
        "LUM",
        "COL3A1",
        "COL6A1",
        "COL6A2",
        "COL6A3",
        "DPT",
        "CFD",
        "PDGFRA",
        "FAP"
    ],

    "Epithelial": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT5",
        "KRT14"
    ],

    "Immune": [
        "PTPRC",
        "CD3D",
        "CD3E",
        "NKG7",
        "LST1",
        "MS4A1",
        "CD79A"
    ],

    "Endothelial": [
        "PECAM1",
        "VWF",
        "EMCN",
        "KDR"
    ],

    "Mural_Pericyte": [
        "RGS5",
        "CSPG4",
        "MCAM",
        "NOTCH3"
    ]
}


core_genes = [
    "CDKN1B",
    "GSN",
    "BGN",
    "TIMP1"
]


all_broad_markers = []

for module, genes in marker_modules.items():

    for gene in genes:

        all_broad_markers.append(
            (
                module,
                gene
            )
        )


if len(
    all_broad_markers
) != 33:

    fail(
        "Frozen broad-marker count is not 33."
    )


print(
    "Frozen broad markers:",
    len(
        all_broad_markers
    )
)

print(
    "CORE genes:",
    len(
        core_genes
    )
)

print()

print(
    "IMPORTANT:"
)

print(
    "KRT5 and KRT14 remain in the frozen biological dictionary."
)

print(
    "They are NOT being deleted."
)

print()


# ============================================================
# 4. Load 27 gene-table paths
# ============================================================

print("===== 4. AUDIT MARKER AVAILABILITY IN ALL 27 FEATURE TABLES =====")
print()

feature_rows = read_csv(
    FEATURE_TABLE_COUNTS
)

inventory = read_csv(
    COMBINED_INV
)


if len(
    feature_rows
) != 27:

    fail(
        "Expected 27 feature-table rows."
    )


if len(
    inventory
) != 27:

    fail(
        "Expected 27 combined tumor rows."
    )


inv_lookup = {
    clean(
        x.get(
            "sample_name"
        )
    ):
    x
    for x in inventory
}


availability_rows = []

missing_events = []


for row in feature_rows:

    sample = clean(
        row.get(
            "sample_name"
        )
    )

    gene_file = (
        PROJECT
        / clean(
            row.get(
                "genes_file"
            )
        )
    )


    if not gene_file.exists():

        fail(
            "Missing gene table:\n"
            + str(
                gene_file
            )
        )


    symbols = read_symbols(
        gene_file
    )


    inv = inv_lookup.get(
        sample
    )

    if inv is None:

        fail(
            "Sample absent from combined inventory: "
            + sample
        )


    for module, gene in all_broad_markers:

        present = (
            gene.upper()
            in symbols
        )


        availability_rows.append({
            "technical_set":
                clean(
                    inv.get(
                        "technical_set"
                    )
                ),

            "patient_id":
                clean(
                    inv.get(
                        "patient_id"
                    )
                ),

            "timepoint":
                clean(
                    inv.get(
                        "timepoint"
                    )
                ),

            "sample_name":
                sample,

            "marker_class":
                "BROAD_ANNOTATION",

            "marker_module":
                module,

            "gene":
                gene,

            "feature_table_status":
                (
                    "MEASURED"
                    if present
                    else
                    "NOT_MEASURED"
                ),

            "absence_interpretation":
                (
                    "NA_NOT_ZERO"
                    if not present
                    else
                    "MEASURED"
                )
        })


        if not present:

            missing_events.append(
                (
                    sample,
                    gene
                )
            )


    for gene in core_genes:

        present = (
            gene.upper()
            in symbols
        )


        availability_rows.append({
            "technical_set":
                clean(
                    inv.get(
                        "technical_set"
                    )
                ),

            "patient_id":
                clean(
                    inv.get(
                        "patient_id"
                    )
                ),

            "timepoint":
                clean(
                    inv.get(
                        "timepoint"
                    )
                ),

            "sample_name":
                sample,

            "marker_class":
                "CORE_PRESENCE_ONLY",

            "marker_module":
                "CORE",

            "gene":
                gene,

            "feature_table_status":
                (
                    "MEASURED"
                    if present
                    else
                    "NOT_MEASURED"
                ),

            "absence_interpretation":
                (
                    "NA_NOT_ZERO"
                    if not present
                    else
                    "MEASURED"
                )
        })


write_csv(
    BY_SAMPLE_OUT,
    availability_rows,
    [
        "technical_set",
        "patient_id",
        "timepoint",
        "sample_name",
        "marker_class",
        "marker_module",
        "gene",
        "feature_table_status",
        "absence_interpretation"
    ]
)


missing_events = sorted(
    set(
        missing_events
    )
)


print(
    "Broad-marker structural missing events:",
    len(
        missing_events
    )
)


for sample, gene in missing_events:

    print(
        " ",
        sample,
        "->",
        gene
    )


print()


# ============================================================
# 5. Confirm exact observed missingness pattern
#
# This is the already-diagnosed technical issue.
# No biological interpretation is allowed.
# ============================================================

print("===== 5. STRUCTURAL-MISSINGNESS PATTERN CHECK =====")
print()


expected_missing = {
    (
        "P5_T_A",
        "KRT5"
    ),
    (
        "P5_T_A",
        "KRT14"
    ),
    (
        "P21_T_B",
        "KRT14"
    )
}


observed_missing = set(
    missing_events
)


print(
    "Expected missing events:",
    len(
        expected_missing
    )
)

print(
    "Observed missing events:",
    len(
        observed_missing
    )
)

print()


if observed_missing != expected_missing:

    print(
        "Unexpected missing events:"
    )

    for x in sorted(
        observed_missing
        -
        expected_missing
    ):

        print(
            " ",
            x
        )

    print(
        "Expected but not observed:"
    )

    for x in sorted(
        expected_missing
        -
        observed_missing
    ):

        print(
            " ",
            x
        )

    fail(
        "Marker missingness differs from diagnosed pattern."
    )


print(
    "✓ Exact diagnosed pattern confirmed:"
)

print(
    "  P5_T_A  : KRT5 + KRT14 NOT MEASURED"
)

print(
    "  P21_T_B : KRT14 NOT MEASURED"
)

print(
    "  Other 25 samples: KRT5/KRT14 measured"
)

print()


# ============================================================
# 6. Verify CORE genes measured in EVERY sample
# ============================================================

print("===== 6. CORE-GENE FEATURE AVAILABILITY =====")
print()


core_missing = [
    x
    for x in availability_rows
    if (
        x[
            "marker_class"
        ]
        == "CORE_PRESENCE_ONLY"
        and
        x[
            "feature_table_status"
        ]
        != "MEASURED"
    )
]


print(
    "CORE missing events:",
    len(
        core_missing
    )
)


if core_missing:

    for x in core_missing:

        print(
            " ",
            x[
                "sample_name"
            ],
            "->",
            x[
                "gene"
            ]
        )

    fail(
        "At least one CORE gene is structurally unavailable."
    )


print(
    "✓ CDKN1B / GSN / BGN / TIMP1 measured in all 27 samples."
)

print()


# ============================================================
# 7. All-27 common-feature availability
# ============================================================

print("===== 7. JOINT ALL-27 MARKER AVAILABILITY =====")
print()


common_rows = read_csv(
    ALL27_COMMON
)


common_symbols = {
    clean(
        x.get(
            "gene_symbol"
        )
    ).upper()
    for x in common_rows
}


joint_lock = []


for module, genes in marker_modules.items():

    for gene in genes:

        in_common = (
            gene.upper()
            in common_symbols
        )


        joint_lock.append({
            "marker_class":
                "BROAD_ANNOTATION",

            "marker_module":
                module,

            "gene":
                gene,

            "biological_dictionary_status":
                "FROZEN_RETAINED",

            "all27_common_feature_status":
                (
                    "MEASURABLE_IN_JOINT_OBJECT"
                    if in_common
                    else
                    "STRUCTURALLY_UNAVAILABLE_IN_JOINT_OBJECT"
                ),

            "joint_annotation_use":
                (
                    "YES"
                    if in_common
                    else
                    "NO_NA_NOT_ZERO"
                ),

            "reason":
                (
                    "PRESENT_IN_ALL27_COMMON_ENSG"
                    if in_common
                    else
                    "FEATURE_TABLE_ABSENT_IN_AT_LEAST_ONE_SAMPLE"
                )
        })


for gene in core_genes:

    in_common = (
        gene.upper()
        in common_symbols
    )


    joint_lock.append({
        "marker_class":
            "CORE_PRESENCE_ONLY",

        "marker_module":
            "CORE",

        "gene":
            gene,

        "biological_dictionary_status":
            "FROZEN_RETAINED",

        "all27_common_feature_status":
            (
                "MEASURABLE_IN_JOINT_OBJECT"
                if in_common
                else
                "STRUCTURALLY_UNAVAILABLE_IN_JOINT_OBJECT"
            ),

        "joint_annotation_use":
            "NOT_USED_FOR_CELL_ANNOTATION",

        "reason":
            (
                "PRESENT_IN_ALL27_COMMON_ENSG"
                if in_common
                else
                "FEATURE_TABLE_ABSENT_IN_AT_LEAST_ONE_SAMPLE"
            )
    })


write_csv(
    JOINT_LOCK_OUT,
    joint_lock,
    [
        "marker_class",
        "marker_module",
        "gene",
        "biological_dictionary_status",
        "all27_common_feature_status",
        "joint_annotation_use",
        "reason"
    ]
)


joint_broad = [
    x
    for x in joint_lock
    if (
        x[
            "marker_class"
        ]
        == "BROAD_ANNOTATION"
    )
]


usable_broad = [
    x
    for x in joint_broad
    if (
        x[
            "joint_annotation_use"
        ]
        == "YES"
    )
]


unavailable_joint = [
    x
    for x in joint_broad
    if (
        x[
            "joint_annotation_use"
        ]
        != "YES"
    )
]


print(
    "Frozen biological marker dictionary:",
    len(
        joint_broad
    ),
    "/ 33 retained"
)

print(
    "Joint all-27 measurable broad markers:",
    len(
        usable_broad
    ),
    "/ 33"
)

print(
    "Structurally unavailable in joint object:",
    len(
        unavailable_joint
    )
)


for x in unavailable_joint:

    print(
        " ",
        x[
            "marker_module"
        ],
        "->",
        x[
            "gene"
        ]
    )


print()


# ============================================================
# 8. Hard-code the technical adaptation
#
# Biological marker definition remains 33.
#
# Joint annotation:
#   Fibroblast  12/12
#   Epithelial   4/6
#   Immune       7/7
#   Endothelial  4/4
#   Mural        4/4
#
# KRT5/KRT14 absence is NEVER interpreted as epithelial-negative.
# ============================================================

print("===== 8. MODULE-COVERAGE LOCK =====")
print()


coverage_rows = []


for module, genes in marker_modules.items():

    measured = [
        g
        for g in genes
        if g.upper()
        in common_symbols
    ]

    unavailable = [
        g
        for g in genes
        if g.upper()
        not in common_symbols
    ]


    coverage_rows.append({
        "marker_module":
            module,

        "frozen_marker_count":
            len(
                genes
            ),

        "joint_measurable_count":
            len(
                measured
            ),

        "joint_measurable_markers":
            ";".join(
                measured
            ),

        "structurally_unavailable_markers":
            ";".join(
                unavailable
            ),

        "module_status":
            (
                "PASS_TECHNICAL_ADAPTATION"
                if (
                    module == "Epithelial"
                    and set(
                        unavailable
                    )
                    == {
                        "KRT5",
                        "KRT14"
                    }
                    and set(
                        measured
                    )
                    == {
                        "EPCAM",
                        "KRT8",
                        "KRT18",
                        "KRT19"
                    }
                )
                else
                (
                    "PASS_FULL_FROZEN_SET"
                    if not unavailable
                    else
                    "FAIL_UNEXPECTED_MARKER_LOSS"
                )
            )
    })


for x in coverage_rows:

    print(
        f"{x['marker_module']:<18s}"
        f"{x['joint_measurable_count']}/"
        f"{x['frozen_marker_count']}  "
        f"{x['module_status']}"
    )


print()


failed_modules = [
    x
    for x in coverage_rows
    if not x[
        "module_status"
    ].startswith(
        "PASS"
    )
]


if failed_modules:

    fail(
        "Unexpected marker loss remains in one or more modules."
    )


print(
    "✓ Fibroblast/Immune/Endothelial/Mural modules unchanged."
)

print(
    "✓ Epithelial module technically evaluable with:"
)

print(
    "  EPCAM / KRT8 / KRT18 / KRT19"
)

print(
    "✓ KRT5/KRT14 remain in frozen dictionary as NOT MEASURED."
)

print(
    "✓ They are never interpreted as zero or epithelial-negative."
)

print()


# ============================================================
# 9. Final scientific decision lock
# ============================================================

print("===== 9. DECISION LOCK =====")
print()


status_rows = [

    {
        "item":
            "frozen_broad_marker_dictionary",

        "value":
            "33_RETAINED"
    },

    {
        "item":
            "joint_measurable_broad_markers",

        "value":
            "31"
    },

    {
        "item":
            "KRT5_status",

        "value":
            "FROZEN_MARKER_STRUCTURALLY_UNAVAILABLE_IN_ALL27_COMMON"
    },

    {
        "item":
            "KRT14_status",

        "value":
            "FROZEN_MARKER_STRUCTURALLY_UNAVAILABLE_IN_ALL27_COMMON"
    },

    {
        "item":
            "structural_absence_interpretation",

        "value":
            "NOT_MEASURED_NA_NEVER_ZERO"
    },

    {
        "item":
            "P5_T_A",

        "value":
            "KRT5_AND_KRT14_NOT_MEASURED"
    },

    {
        "item":
            "P21_T_B",

        "value":
            "KRT14_NOT_MEASURED"
    },

    {
        "item":
            "sample_exclusion_due_to_KRT_missingness",

        "value":
            "NONE"
    },

    {
        "item":
            "epithelial_joint_markers",

        "value":
            "EPCAM;KRT8;KRT18;KRT19"
    },

    {
        "item":
            "CORE_genes_all27",

        "value":
            "4_OF_4_PRESENT"
    },

    {
        "item":
            "CORE_scoring_impact",

        "value":
            "NONE_TECHNICALLY"
    },

    {
        "item":
            "CORE_values_read",

        "value":
            "NO"
    },

    {
        "item":
            "histology_mapping",

        "value":
            "INCOMPLETE"
    },

    {
        "item":
            "primary_CORE2_inference",

        "value":
            "BLOCKED"
    },

    {
        "item":
            "STEP22A_06B_final_status",

        "value":
            "PASS_WITH_FROZEN_STRUCTURAL_MARKER_MISSINGNESS"
    }
]


write_csv(
    STATUS_OUT,
    status_rows,
    [
        "item",
        "value"
    ]
)


with open(
    DECISION_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-6B1 STRUCTURAL MARKER MISSINGNESS DECISION\n"
    )

    f.write(
        "=" * 70
        + "\n\n"
    )

    f.write(
        "BIOLOGICAL MARKER DICTIONARY\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "The frozen 33-marker biological dictionary is NOT changed.\n"
    )

    f.write(
        "KRT5 and KRT14 remain frozen epithelial-exclusion markers.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "STRUCTURAL FEATURE-TABLE MISSINGNESS\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "P5_T_A  : KRT5 and KRT14 are NOT MEASURED.\n"
    )

    f.write(
        "P21_T_B : KRT14 is NOT MEASURED.\n"
    )

    f.write(
        "Other 25 tumors: KRT5 and KRT14 are measured.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "These absences are feature-table absences.\n"
    )

    f.write(
        "They are NOT biological zero expression.\n"
    )

    f.write(
        "They are NOT evidence against epithelial identity.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "JOINT 27-TUMOR ANNOTATION\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "Use only frozen markers measurable in the true all-27\n"
    )

    f.write(
        "common ENSG universe.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Fibroblast_CAF : 12/12 markers\n"
    )

    f.write(
        "Epithelial     : 4/6 markers\n"
    )

    f.write(
        "                 EPCAM, KRT8, KRT18, KRT19\n"
    )

    f.write(
        "Immune         : 7/7 markers\n"
    )

    f.write(
        "Endothelial    : 4/4 markers\n"
    )

    f.write(
        "Mural_Pericyte : 4/4 markers\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "KRT5/KRT14 remain in provenance as structurally unavailable.\n"
    )

    f.write(
        "No sample is excluded for this reason.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "CORE SCORE TECHNICAL AVAILABILITY\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "CDKN1B / GSN / BGN / TIMP1 = 4/4 available in all 27.\n"
    )

    f.write(
        "No CORE expression values were read.\n"
    )

    f.write(
        "No CORE score was calculated.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "HISTOLOGY\n"
    )

    f.write(
        "------------------------------------------------------------\n"
    )

    f.write(
        "EASC identity remains unresolved.\n"
    )

    f.write(
        "Primary CORE2 inference remains blocked.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "FINAL STEP22A-6B STATUS:\n"
    )

    f.write(
        "PASS_WITH_FROZEN_STRUCTURAL_MARKER_MISSINGNESS\n"
    )


# ============================================================
# 10. Final frozen integrity
# ============================================================

print("===== 10. FINAL FROZEN PROJECT CHECK =====")
print()

changes_after = frozen_integrity()

print(
    "Frozen original files changed:",
    len(
        changes_after
    )
)

if changes_after:

    for x in changes_after[:20]:
        print(x)

    fail(
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
print("STEP22A-6B1 COMPLETE")
print("=" * 80)
print()

print(
    "FROZEN BIOLOGICAL MARKERS    : 33 / 33 RETAINED"
)

print(
    "JOINT MEASURABLE MARKERS     : 31 / 33"
)

print()

print(
    "STRUCTURAL MISSINGNESS:"
)

print(
    "  P5_T_A  : KRT5, KRT14"
)

print(
    "  P21_T_B : KRT14"
)

print()

print(
    "ABSENCE INTERPRETATION       : NOT_MEASURED / NA"
)

print(
    "ABSENCE TREATED AS ZERO      : NO"
)

print(
    "SAMPLES EXCLUDED             : 0"
)

print()

print(
    "JOINT EPITHELIAL MARKERS:"
)

print(
    "  EPCAM / KRT8 / KRT18 / KRT19"
)

print()

print(
    "FIBROBLAST MARKERS           : 12 / 12"
)

print(
    "IMMUNE MARKERS               : 7 / 7"
)

print(
    "ENDOTHELIAL MARKERS          : 4 / 4"
)

print(
    "MURAL/PERICYTE MARKERS       : 4 / 4"
)

print()

print(
    "CORE GENES                   : 4 / 4"
)

print(
    "CORE TECHNICAL IMPACT        : NONE"
)

print(
    "CORE GENE VALUES READ        : NO"
)

print(
    "CORE2 SCORING                : NOT PERFORMED"
)

print()

print(
    "STEP22A-6B FINAL STATUS:"
)

print(
    "PASS_WITH_FROZEN_STRUCTURAL_MARKER_MISSINGNESS"
)

print()

print(
    "HISTOLOGY MAPPING            : INCOMPLETE"
)

print(
    "PRIMARY CORE2 INFERENCE      : BLOCKED"
)

print()

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
    "STEP22A-6C = APPLY THE FROZEN QC +"
)

print(
    "31-MEASURABLE-MARKER JOINT BROAD ANNOTATION"
)

print(
    "TO ALL 27 TUMOR OBSERVATIONS."
)

print()

print(
    "DO NOT SCORE CORE2."
)

print()
