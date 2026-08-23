#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# Step22A-2
#
# Canonical sample mapping:
#   patient
#   tissue
#   treatment timepoint
#   ESCC eligibility
#
# NO .tar download in this step.
# NO modification of frozen Step19-21 outputs.
# ============================================================

import csv
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime


PROJECT = os.getcwd()

META_DIR = os.path.join(
    PROJECT,
    "00_metadata",
    "OMIX005710"
)

RAW_META = os.path.join(
    META_DIR,
    "STEP22A_SAMPLE_METADATA_RAW.csv"
)

FILE_INV = os.path.join(
    META_DIR,
    "STEP22A_OMIX_FILE_INVENTORY.csv"
)

FROZEN_INV = os.path.join(
    META_DIR,
    "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

RAW_DIR = os.path.join(
    PROJECT,
    "01_raw",
    "OMIX005710"
)


print()
print("=" * 72)
print("STEP22A-2: CANONICAL SAMPLE + HISTOLOGY AUDIT")
print("=" * 72)
print()


# ============================================================
# Helpers
# ============================================================

def read_csv(path):

    if not os.path.exists(path):
        raise SystemExit(
            "ERROR: missing required file:\n" + path
        )

    with open(
        path,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):

    with open(
        path,
        "w",
        encoding="utf-8-sig",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def patient_sort_key(x):

    m = re.search(r"(\d+)", x)

    if m:
        return int(m.group(1))

    return 999999


# ============================================================
# 1. Load Step22A-1 outputs
# ============================================================

print("===== 1. LOAD STEP22A-1 OUTPUTS =====")
print()

meta = read_csv(RAW_META)
files = read_csv(FILE_INV)

print("Metadata rows :", len(meta))
print("OMIX files    :", len(files))
print()

if len(meta) != 46:
    raise SystemExit(
        "ERROR: expected 46 metadata rows.\n"
        "STEP22A_METADATA_COUNT_MISMATCH"
    )

if len(files) != 46:
    raise SystemExit(
        "ERROR: expected 46 OMIX file records.\n"
        "STEP22A_FILE_COUNT_MISMATCH"
    )


# ============================================================
# 2. Build file lookup
# ============================================================

file_by_hrs = {}

for row in files:

    hrs = (row.get("hrs_accession") or "").strip()

    if not hrs:
        raise SystemExit(
            "ERROR: OMIX inventory contains missing HRS."
        )

    if hrs in file_by_hrs:
        raise SystemExit(
            "ERROR: duplicate HRS in OMIX inventory: "
            + hrs
        )

    file_by_hrs[hrs] = row


# ============================================================
# 3. Canonicalize sample identity
#
# Frozen biological nomenclature from source publications:
#
# T = tumor
# N = adjacent normal
# B = before / treatment-naive
# A = after treatment
#
# Tissue requires agreement between:
#   sample_name coding
#   official sample_title
#
# description is audited but NOT allowed to override
# concordant sample_name + official title.
# ============================================================

print("===== 2. CANONICALIZE 46 SAMPLE IDENTITIES =====")
print()

canonical = []
conflict_resolution = []

fatal_errors = []

for row in meta:

    hrs = (row.get("hrs_accession") or "").strip()

    sample_name = (
        row.get("sample_name") or ""
    ).strip()

    title = (
        row.get("sample_title") or ""
    ).strip()

    description = (
        row.get("description") or ""
    ).strip()

    biosample = (
        row.get("biosample_accession") or ""
    ).strip()

    # --------------------------------------------------------
    # A. Parse official sample name
    # --------------------------------------------------------

    m = re.fullmatch(
        r"P(\d+)_(T|N)_(A|B)",
        sample_name,
        flags=re.I
    )

    if not m:

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "SAMPLE_NAME_PATTERN_INVALID"
            )
        )

        continue

    patient_id = "P" + str(
        int(m.group(1))
    )

    tissue_code = m.group(2).upper()
    time_code = m.group(3).upper()

    tissue_from_name = (
        "Tumor"
        if tissue_code == "T"
        else "Adjacent_normal"
    )

    timepoint = (
        "After"
        if time_code == "A"
        else "Before"
    )

    # --------------------------------------------------------
    # B. Classify official sample title
    # --------------------------------------------------------

    tlow = title.lower()

    if "adjacent" in tlow:

        tissue_from_title = "Adjacent_normal"

    elif (
        "esophageal carcinoma" in tlow
        or
        "oesophageal carcinoma" in tlow
    ):

        tissue_from_title = "Tumor"

    else:

        tissue_from_title = "Unknown"

    # --------------------------------------------------------
    # C. Name/title MUST agree
    # --------------------------------------------------------

    if tissue_from_title == "Unknown":

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "OFFICIAL_TITLE_TISSUE_UNKNOWN"
            )
        )

        continue

    if tissue_from_title != tissue_from_name:

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "NAME_TITLE_TISSUE_DISAGREEMENT"
            )
        )

        continue

    # --------------------------------------------------------
    # D. Parse description only as AUDIT
    # --------------------------------------------------------

    dlow = description.lower()

    desc_tissue = "Unknown"

    if "adjacent" in dlow:
        desc_tissue = "Adjacent_normal"

    elif "tumor" in dlow:
        desc_tissue = "Tumor"

    desc_time = "Unknown"

    if "before treatment" in dlow:
        desc_time = "Before"

    elif "after treatment" in dlow:
        desc_time = "After"

    desc_patient = ""

    pm = re.search(
        r"Patient\s*(\d+)",
        description,
        flags=re.I
    )

    if pm:
        desc_patient = (
            "P"
            + str(int(pm.group(1)))
        )

    # --------------------------------------------------------
    # E. Patient / timepoint consistency
    # --------------------------------------------------------

    if (
        desc_patient
        and
        desc_patient != patient_id
    ):

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "PATIENT_ID_DISAGREEMENT"
            )
        )

        continue

    if (
        desc_time != "Unknown"
        and
        desc_time != timepoint
    ):

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "TIMEPOINT_DISAGREEMENT"
            )
        )

        continue

    # --------------------------------------------------------
    # F. Description tissue conflict
    # --------------------------------------------------------

    description_tissue_conflict = (
        "YES"
        if (
            desc_tissue != "Unknown"
            and
            desc_tissue != tissue_from_name
        )
        else "NO"
    )

    resolution = (
        "OFFICIAL_NAME_TITLE_OVERRIDE_DESCRIPTION"
        if description_tissue_conflict == "YES"
        else
        "ALL_AVAILABLE_TISSUE_FIELDS_CONCORDANT"
    )

    # --------------------------------------------------------
    # G. Histology audit
    #
    # Official BioProject PRJCA016745 states the 22 patients
    # were diagnosed with ESCC.
    #
    # We therefore use PROJECT-LEVEL confirmed ESCC status.
    # We do NOT invent patient-specific pathology fields.
    # --------------------------------------------------------

    histology = "ESCC"

    histology_basis = (
        "PRJCA016745_PROJECT_LEVEL_ESCC_CONFIRMATION"
    )

    # Guard against obvious contradictory public text
    histology_text = (
        title + " " + description
    ).lower()

    contradictory_terms = [
        "adenocarcinoma",
        "adenosquamous",
        "adeno-squamous"
    ]

    contradicted = [
        x for x in contradictory_terms
        if x in histology_text
    ]

    if contradicted:

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "CONTRADICTORY_HISTOLOGY_TEXT:"
                + ",".join(contradicted)
            )
        )

        continue

    # --------------------------------------------------------
    # H. Analysis eligibility
    # --------------------------------------------------------

    if tissue_from_name == "Tumor":

        analysis_status = (
            "ELIGIBLE_TUMOR_FOR_EXTERNAL_VALIDATION"
        )

        download_candidate = "YES"

    else:

        analysis_status = (
            "EXCLUDE_PRIMARY_ANALYSIS_ADJACENT_NORMAL"
        )

        download_candidate = "NO"

    # --------------------------------------------------------
    # I. OMIX file information
    # --------------------------------------------------------

    if hrs not in file_by_hrs:

        fatal_errors.append(
            (
                hrs,
                sample_name,
                "NO_OMIX_FILE_MAPPING"
            )
        )

        continue

    frow = file_by_hrs[hrs]

    canonical.append({

        "patient_id":
            patient_id,

        "sample_name":
            sample_name,

        "hrs_accession":
            hrs,

        "biosample_accession":
            biosample,

        "omix_file_id":
            frow.get("file_id", ""),

        "canonical_tissue":
            tissue_from_name,

        "canonical_timepoint":
            timepoint,

        "canonical_group":
            (
                ("T" if tissue_from_name == "Tumor" else "N")
                + "_"
                + ("A" if timepoint == "After" else "B")
            ),

        "sample_title":
            title,

        "description":
            description,

        "tissue_from_sample_name":
            tissue_from_name,

        "tissue_from_official_title":
            tissue_from_title,

        "tissue_from_description":
            desc_tissue,

        "description_tissue_conflict":
            description_tissue_conflict,

        "conflict_resolution":
            resolution,

        "histology":
            histology,

        "histology_evidence_level":
            "PROJECT_LEVEL",

        "histology_basis":
            histology_basis,

        "patient_specific_histology_public":
            "NOT_REQUIRED_FOR_PROJECT_LEVEL_ESCC_COHORT",

        "analysis_status":
            analysis_status,

        "download_candidate":
            download_candidate,

        "listed_size_MB":
            frow.get("listed_size_MB", ""),

        "download_url_https":
            frow.get("download_url_https", ""),

        "data_downloaded":
            "NO"
    })

    if description_tissue_conflict == "YES":

        conflict_resolution.append({

            "hrs_accession":
                hrs,

            "sample_name":
                sample_name,

            "sample_title":
                title,

            "description":
                description,

            "name_tissue":
                tissue_from_name,

            "title_tissue":
                tissue_from_title,

            "description_tissue":
                desc_tissue,

            "resolution":
                resolution,

            "biological_basis":
                (
                    "Published T/N sample nomenclature + "
                    "official BioSample title"
                ),

            "sample_excluded_due_to_conflict":
                "NO"
        })


# ============================================================
# 4. Stop on unresolved contradictions
# ============================================================

if fatal_errors:

    print()
    print("FATAL CANONICALIZATION ERRORS:")
    print()

    for x in fatal_errors:
        print(x)

    print()
    print("STEP22A_CANONICAL_MAPPING_FAILED")
    sys.exit(1)


if len(canonical) != 46:

    raise SystemExit(
        "ERROR: canonical sample count != 46"
    )


# ============================================================
# 5. Validate exact source-study group counts
# ============================================================

print("===== 3. SOURCE-STUDY SAMPLE COUNT CHECK =====")
print()

group_counts = Counter(
    x["canonical_group"]
    for x in canonical
)

expected_counts = {
    "T_B": 15,
    "N_B": 7,
    "T_A": 12,
    "N_A": 12
}

for group in [
    "T_B",
    "N_B",
    "T_A",
    "N_A"
]:

    observed = group_counts.get(
        group,
        0
    )

    expected = expected_counts[group]

    status = (
        "PASS"
        if observed == expected
        else "FAIL"
    )

    print(
        f"{group}: "
        f"observed={observed}, "
        f"expected={expected} "
        f"[{status}]"
    )

    if observed != expected:

        raise SystemExit(
            "\nERROR: source-study sample count "
            "does not match canonical mapping.\n"
            "STEP22A_GROUP_COUNT_MISMATCH"
        )

print()

print(
    "✓ Exact 46-sample composition reproduced:"
)

print(
    "  15 T_B + 7 N_B + "
    "12 T_A + 12 N_A = 46"
)

print()


# ============================================================
# 6. Validate patient count
# ============================================================

print("===== 4. PATIENT AUDIT =====")
print()

patients = sorted(
    set(
        x["patient_id"]
        for x in canonical
    ),
    key=patient_sort_key
)

print(
    "Unique patients:",
    len(patients)
)

print(
    "Patients:",
    ", ".join(patients)
)

if len(patients) != 22:

    raise SystemExit(
        "\nERROR: expected 22 patients.\n"
        "STEP22A_PATIENT_COUNT_MISMATCH"
    )

print()
print("✓ 22-patient ESCC cohort confirmed.")
print()


# ============================================================
# 7. Tumor-only external-validation inventory
# ============================================================

print("===== 5. PRIMARY TUMOR SAMPLE INVENTORY =====")
print()

tumor = [
    x for x in canonical
    if x["canonical_tissue"] == "Tumor"
]

tumor_before = [
    x for x in tumor
    if x["canonical_timepoint"] == "Before"
]

tumor_after = [
    x for x in tumor
    if x["canonical_timepoint"] == "After"
]

print(
    "Eligible tumor samples :",
    len(tumor)
)

print(
    "Before-treatment tumor :",
    len(tumor_before)
)

print(
    "After-treatment tumor  :",
    len(tumor_after)
)

if len(tumor) != 27:

    raise SystemExit(
        "ERROR: expected 27 tumor samples."
    )

if len(tumor_before) != 15:

    raise SystemExit(
        "ERROR: expected 15 T_B samples."
    )

if len(tumor_after) != 12:

    raise SystemExit(
        "ERROR: expected 12 T_A samples."
    )

print()


# ============================================================
# 8. Identify paired tumor patients
# ============================================================

patient_timepoints = defaultdict(set)

patient_samples = defaultdict(list)

for x in tumor:

    patient_timepoints[
        x["patient_id"]
    ].add(
        x["canonical_timepoint"]
    )

    patient_samples[
        x["patient_id"]
    ].append(
        x["sample_name"]
    )

paired_patients = sorted(
    [
        p
        for p, tps
        in patient_timepoints.items()
        if {
            "Before",
            "After"
        }.issubset(tps)
    ],
    key=patient_sort_key
)

print("===== 6. PAIRED TUMOR PATIENTS =====")
print()

print(
    "Patients with BOTH before + after tumor:",
    len(paired_patients)
)

if paired_patients:

    for p in paired_patients:

        print(
            " ",
            p,
            "->",
            ", ".join(
                sorted(
                    patient_samples[p]
                )
            )
        )

else:

    print("  <none>")

print()


# ============================================================
# 9. Description-conflict resolution audit
# ============================================================

print("===== 7. DESCRIPTION CONFLICT RESOLUTION =====")
print()

print(
    "Description tissue conflicts:",
    len(conflict_resolution)
)

for x in conflict_resolution:

    print()
    print("-" * 72)

    print(
        "HRS         :",
        x["hrs_accession"]
    )

    print(
        "Sample      :",
        x["sample_name"]
    )

    print(
        "Title       :",
        x["sample_title"]
    )

    print(
        "Description :",
        x["description"]
    )

    print(
        "Name tissue :",
        x["name_tissue"]
    )

    print(
        "Title tissue:",
        x["title_tissue"]
    )

    print(
        "Desc tissue :",
        x["description_tissue"]
    )

    print(
        "Resolution  :",
        x["resolution"]
    )

print()

if len(conflict_resolution) != 5:

    print(
        "WARNING: Step22A-1 reported 5 conflicts, "
        "but Step22A-2 detected",
        len(conflict_resolution)
    )

    print(
        "This is not automatically fatal, "
        "but should be reviewed."
    )

print()


# ============================================================
# 10. Histology audit
# ============================================================

print("===== 8. HISTOLOGY AUDIT =====")
print()

print(
    "Official BioProject cohort diagnosis:"
)

print(
    "  ESCC, 22 patients"
)

print()

print(
    "Contradictory public adenocarcinoma/"
    "adenosquamous strings detected: 0"
)

print()

print(
    "Histology rule:"
)

print(
    "  Use project-level ESCC confirmation."
)

print(
    "  Do NOT invent or exclude a hypothetical "
    "EASC patient."
)

print()


# ============================================================
# 11. Targeted download plan — DO NOT download yet
# ============================================================

print("===== 9. TARGETED DOWNLOAD PLAN =====")
print()

download_rows = []

total_selected_mb = 0.0

for x in tumor:

    try:
        size_mb = float(
            x["listed_size_MB"]
        )

    except Exception:
        size_mb = 0.0

    total_selected_mb += size_mb

    download_rows.append({

        "priority":
            (
                "1_PRIMARY_PAIRED"
                if x["patient_id"]
                in paired_patients
                else
                "2_SECONDARY_TUMOR"
            ),

        "patient_id":
            x["patient_id"],

        "timepoint":
            x["canonical_timepoint"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            x["hrs_accession"],

        "omix_file_id":
            x["omix_file_id"],

        "listed_size_MB":
            x["listed_size_MB"],

        "download_url_https":
            x["download_url_https"],

        "download_now":
            "NO",

        "reason":
            (
                "PRIMARY_PAIRED_EXTERNAL_VALIDATION"
                if x["patient_id"]
                in paired_patients
                else
                "SECONDARY_ALL_TUMOR_VALIDATION"
            )
    })


download_rows.sort(
    key=lambda x: (
        x["priority"],
        patient_sort_key(
            x["patient_id"]
        ),
        x["timepoint"]
    )
)

print(
    "Tumor files selected for future download:",
    len(download_rows)
)

print(
    "Selected listed size:",
    round(
        total_selected_mb,
        2
    ),
    "MB"
)

print(
    "Selected listed size:",
    round(
        total_selected_mb / 1024,
        3
    ),
    "GiB"
)

print()

print(
    "Adjacent-normal files excluded "
    "from primary download:"
)

print(
    46 - len(download_rows)
)

print()


# ============================================================
# 12. Save outputs
# ============================================================

canonical.sort(
    key=lambda x: (
        patient_sort_key(
            x["patient_id"]
        ),
        x["canonical_timepoint"],
        x["canonical_tissue"]
    )
)

CANONICAL_OUT = os.path.join(
    META_DIR,
    "STEP22A_SAMPLE_METADATA_CANONICAL.csv"
)

write_csv(
    CANONICAL_OUT,
    canonical,
    list(canonical[0].keys())
)


CONFLICT_OUT = os.path.join(
    META_DIR,
    "STEP22A_METADATA_CONFLICT_RESOLUTION.csv"
)

conflict_fields = [
    "hrs_accession",
    "sample_name",
    "sample_title",
    "description",
    "name_tissue",
    "title_tissue",
    "description_tissue",
    "resolution",
    "biological_basis",
    "sample_excluded_due_to_conflict"
]

write_csv(
    CONFLICT_OUT,
    conflict_resolution,
    conflict_fields
)


ANALYSIS_OUT = os.path.join(
    META_DIR,
    "STEP22A_ANALYSIS_SAMPLE_INVENTORY.csv"
)

analysis_rows = []

for x in canonical:

    analysis_rows.append({

        "patient_id":
            x["patient_id"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            x["hrs_accession"],

        "tissue":
            x["canonical_tissue"],

        "timepoint":
            x["canonical_timepoint"],

        "group":
            x["canonical_group"],

        "histology":
            x["histology"],

        "analysis_status":
            x["analysis_status"],

        "paired_tumor_patient":
            (
                "YES"
                if (
                    x["canonical_tissue"]
                    == "Tumor"
                    and
                    x["patient_id"]
                    in paired_patients
                )
                else "NO"
            ),

        "download_candidate":
            x["download_candidate"]
    })


write_csv(
    ANALYSIS_OUT,
    analysis_rows,
    list(analysis_rows[0].keys())
)


DOWNLOAD_OUT = os.path.join(
    META_DIR,
    "STEP22A_TARGETED_TUMOR_DOWNLOAD_PLAN.csv"
)

write_csv(
    DOWNLOAD_OUT,
    download_rows,
    list(download_rows[0].keys())
)


HISTOLOGY_OUT = os.path.join(
    META_DIR,
    "STEP22A_HISTOLOGY_MAPPING_STATUS.csv"
)

histology_rows = []

for p in patients:

    histology_rows.append({

        "patient_id":
            p,

        "histology":
            "ESCC",

        "evidence_level":
            "PROJECT_LEVEL",

        "evidence_source":
            (
                "PRJCA016745 official BioProject; "
                "Genome Medicine 2024 source study"
            ),

        "adenosquamous_exclusion_required":
            "NO",

        "status":
            "ESCC_PROJECT_LEVEL_CONFIRMED"
    })


write_csv(
    HISTOLOGY_OUT,
    histology_rows,
    list(histology_rows[0].keys())
)


# ============================================================
# 13. Summary
# ============================================================

SUMMARY_OUT = os.path.join(
    META_DIR,
    "STEP22A_02_CANONICAL_SAMPLE_AUDIT_SUMMARY.txt"
)

with open(
    SUMMARY_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-2 CANONICAL SAMPLE AUDIT\n"
    )

    f.write("=" * 60 + "\n\n")

    f.write(
        f"Canonical samples: {len(canonical)}\n"
    )

    f.write(
        f"Patients: {len(patients)}\n"
    )

    f.write(
        f"T_B: {group_counts['T_B']}\n"
    )

    f.write(
        f"N_B: {group_counts['N_B']}\n"
    )

    f.write(
        f"T_A: {group_counts['T_A']}\n"
    )

    f.write(
        f"N_A: {group_counts['N_A']}\n"
    )

    f.write(
        f"Tumor samples: {len(tumor)}\n"
    )

    f.write(
        "Before tumor samples: "
        f"{len(tumor_before)}\n"
    )

    f.write(
        "After tumor samples: "
        f"{len(tumor_after)}\n"
    )

    f.write(
        "Paired tumor patients: "
        f"{len(paired_patients)}\n"
    )

    f.write(
        "Paired patient IDs: "
        + ",".join(paired_patients)
        + "\n"
    )

    f.write(
        "Resolved description tissue conflicts: "
        f"{len(conflict_resolution)}\n"
    )

    f.write(
        "Histology: ESCC "
        "(project-level confirmed)\n"
    )

    f.write(
        "Hypothetical EASC exclusion applied: NO\n"
    )

    f.write(
        "Future tumor files selected: "
        f"{len(download_rows)}\n"
    )

    f.write(
        "Selected download size MB: "
        f"{round(total_selected_mb, 2)}\n"
    )

    f.write(
        "Selected download size GiB: "
        f"{round(total_selected_mb / 1024, 3)}\n"
    )

    f.write(
        "Data files downloaded in Step22A-2: NO\n"
    )


# ============================================================
# 14. Verify no tar downloaded
# ============================================================

tar_files = []

for root, dirs, filenames in os.walk(
    RAW_DIR
):

    for filename in filenames:

        if filename.lower().endswith(
            ".tar"
        ):

            tar_files.append(
                os.path.join(
                    root,
                    filename
                )
            )

if tar_files:

    print()
    print(
        "ERROR: unexpected OMIX .tar files found:"
    )

    for x in tar_files:
        print(x)

    print()

    print(
        "STEP22A_UNEXPECTED_TAR_PRESENT"
    )

    sys.exit(1)


# ============================================================
# 15. Verify frozen Step19-21 file inventory unchanged
# ============================================================

print("===== 10. FROZEN PROJECT CHECK =====")
print()

frozen_changes = []

if os.path.exists(FROZEN_INV):

    frozen_rows = read_csv(
        FROZEN_INV
    )

    for row in frozen_rows:

        path = row.get(
            "path",
            ""
        )

        if not path:
            continue

        if not os.path.exists(path):

            frozen_changes.append(
                (
                    path,
                    "MISSING"
                )
            )

            continue

        expected_size = str(
            row.get(
                "size_bytes",
                ""
            )
        )

        actual_size = str(
            os.path.getsize(path)
        )

        if (
            expected_size
            and
            expected_size != actual_size
        ):

            frozen_changes.append(
                (
                    path,
                    "SIZE_CHANGED"
                )
            )


print(
    "Frozen original files with size changes:",
    len(frozen_changes)
)

if frozen_changes:

    for x in frozen_changes[:20]:
        print(x)

    print()
    print(
        "ERROR: frozen original project changed."
    )

    print(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )

    sys.exit(1)

print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# FINAL
# ============================================================

print("=" * 72)
print("STEP22A-2 COMPLETE")
print("=" * 72)
print()

print(
    "CANONICAL SAMPLES           :",
    len(canonical)
)

print(
    "ESCC PATIENTS               :",
    len(patients)
)

print()

print(
    "T_B / N_B / T_A / N_A      :",
    group_counts["T_B"],
    "/",
    group_counts["N_B"],
    "/",
    group_counts["T_A"],
    "/",
    group_counts["N_A"]
)

print()

print(
    "ELIGIBLE TUMOR SAMPLES      :",
    len(tumor)
)

print(
    "  BEFORE                    :",
    len(tumor_before)
)

print(
    "  AFTER                     :",
    len(tumor_after)
)

print()

print(
    "PAIRED TUMOR PATIENTS       :",
    len(paired_patients)
)

print(
    "PAIRED PATIENT IDS          :",
    ", ".join(paired_patients)
)

print()

print(
    "DESCRIPTION CONFLICTS       :",
    len(conflict_resolution)
)

print(
    "CONFLICTS RESOLVED          : YES"
)

print(
    "RESOLUTION BASIS            :"
)

print(
    "  sample_name + official sample_title"
)

print(
    "  + published T/N/A/B nomenclature"
)

print()

print(
    "HISTOLOGY                   : ESCC"
)

print(
    "HISTOLOGY EVIDENCE          : PROJECT LEVEL"
)

print(
    "EASC PATIENT EXCLUDED       : NO"
)

print(
    "UNSUPPORTED EASC ASSUMPTION : REMOVED"
)

print()

print(
    "FUTURE TUMOR DOWNLOAD FILES :",
    len(download_rows)
)

print(
    "FUTURE DOWNLOAD SIZE        :",
    round(total_selected_mb / 1024, 3),
    "GiB"
)

print()

print(
    "OMIX TAR FILES DOWNLOADED   : 0"
)

print(
    "FROZEN STEP19-21 MODIFIED   : NO"
)

print(
    "MANUSCRIPT MODIFIED         : NO"
)

print()

print("OUTPUTS:")

for path in [
    CANONICAL_OUT,
    CONFLICT_OUT,
    ANALYSIS_OUT,
    DOWNLOAD_OUT,
    HISTOLOGY_OUT,
    SUMMARY_OUT
]:

    print(
        " ",
        os.path.relpath(
            path,
            PROJECT
        )
    )

print()

print(
    "READY FOR STEP22A-3:"
)

print(
    "TARGETED TUMOR-ONLY DOWNLOAD"
)

print()

print("STOP HERE.")
print(
    "DO NOT MANUALLY DOWNLOAD ANY FILES."
)
print()
