#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-4B
#
# PUBLIC GSA-Human HRI HISTOLOGY AUDIT
#
# Purpose:
# - Map sample -> patient -> HRI individual
# - Fetch PUBLIC individualDetail pages only
# - Search for explicit EASC / adenosquamous evidence
#
# IMPORTANT:
# - PUBLIC pages only
# - NO controlled-access circumvention
# - NO OMIX TAR download
# - NO CAF annotation
# - NO CORE2 scoring
# - NO exclusion unless a UNIQUE public mapping is found
# - NO changes to frozen Step19-21
# ============================================================

import csv
import html
import os
import re
import subprocess
import sys
from collections import defaultdict
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

FROZEN_INV = os.path.join(
    META_DIR,
    "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

SOURCE_DIR = os.path.join(
    META_DIR,
    "STEP22A_04B_public_HRI_pages"
)

os.makedirs(
    SOURCE_DIR,
    exist_ok=True
)


print()
print("=" * 78)
print("STEP22A-4B: PUBLIC HRI HISTOLOGY AUDIT")
print("=" * 78)
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


def strip_html(raw):

    raw = re.sub(
        r"<script\b.*?</script>",
        " ",
        raw,
        flags=re.I | re.S
    )

    raw = re.sub(
        r"<style\b.*?</style>",
        " ",
        raw,
        flags=re.I | re.S
    )

    raw = re.sub(
        r"<[^>]+>",
        " ",
        raw
    )

    raw = html.unescape(raw)

    raw = raw.replace(
        "\xa0",
        " "
    )

    raw = re.sub(
        r"\s+",
        " ",
        raw
    )

    return raw.strip()


def fetch_public_page(url, outfile):

    tmp = outfile + ".part"

    if os.path.exists(tmp):
        os.remove(tmp)

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--silent",
        "--show-error",
        "--compressed",
        "--retry", "4",
        "--retry-delay", "2",
        "--connect-timeout", "20",
        "--max-time", "90",
        "-A",
        "Mozilla/5.0 Step22A-public-HRI-audit",
        url,
        "-o",
        tmp
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:

        if os.path.exists(tmp):
            os.remove(tmp)

        return (
            False,
            result.stderr.strip()
        )

    if (
        not os.path.exists(tmp)
        or
        os.path.getsize(tmp) < 100
    ):

        if os.path.exists(tmp):
            os.remove(tmp)

        return (
            False,
            "PAGE_TOO_SMALL"
        )

    os.replace(
        tmp,
        outfile
    )

    with open(
        outfile,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        return (
            True,
            f.read()
        )


def patient_from_sample_name(sample_name):

    m = re.fullmatch(
        r"P(\d+)_[TN]_[AB]",
        sample_name,
        flags=re.I
    )

    if not m:
        return ""

    return "P" + str(
        int(m.group(1))
    )


def patient_key(x):

    m = re.search(
        r"\d+",
        x
    )

    if m:
        return int(
            m.group()
        )

    return 999999


# ============================================================
# 1. Read Step22A metadata
# ============================================================

print("===== 1. LOAD EXISTING SAMPLE METADATA =====")
print()

rows = read_csv(
    RAW_META
)

print(
    "Sample metadata rows:",
    len(rows)
)

if len(rows) != 46:

    raise SystemExit(
        "ERROR: expected 46 sample metadata rows."
    )

print()


# ============================================================
# 2. Build patient ↔ HRI mapping
# ============================================================

print("===== 2. BUILD PATIENT ↔ HRI MAPPING =====")
print()

mapping_rows = []

patient_to_hri = defaultdict(set)
hri_to_patient = defaultdict(set)

for row in rows:

    sample_name = (
        row.get("sample_name")
        or ""
    ).strip()

    hrs = (
        row.get("hrs_accession")
        or ""
    ).strip()

    hri = (
        row.get("individual_accession")
        or ""
    ).strip()

    patient = patient_from_sample_name(
        sample_name
    )

    if not patient:

        raise SystemExit(
            "ERROR: cannot parse patient from "
            + sample_name
        )

    if not re.fullmatch(
        r"HRI\d+",
        hri
    ):

        raise SystemExit(
            "ERROR: invalid/missing HRI for "
            + sample_name
            + ": "
            + hri
        )

    patient_to_hri[
        patient
    ].add(
        hri
    )

    hri_to_patient[
        hri
    ].add(
        patient
    )

    mapping_rows.append({

        "patient_id":
            patient,

        "sample_name":
            sample_name,

        "hrs_accession":
            hrs,

        "hri_accession":
            hri
    })


patients = sorted(
    patient_to_hri,
    key=patient_key
)

unique_hri = sorted(
    hri_to_patient
)


print(
    "Unique patients:",
    len(patients)
)

print(
    "Unique HRI individuals:",
    len(unique_hri)
)

print()


# HRI must never map across multiple patient IDs.

cross_patient_hri = {

    hri: vals

    for hri, vals
    in hri_to_patient.items()

    if len(vals) != 1
}

if cross_patient_hri:

    print(
        "ERROR: one HRI maps to multiple P IDs:"
    )

    for hri, vals in cross_patient_hri.items():

        print(
            hri,
            sorted(vals)
        )

    print()

    print(
        "STEP22A_HRI_PATIENT_MAPPING_INVALID"
    )

    sys.exit(1)


mapping_out = os.path.join(
    META_DIR,
    "STEP22A_04B_PATIENT_HRI_MAPPING.csv"
)

write_csv(
    mapping_out,
    mapping_rows,
    [
        "patient_id",
        "sample_name",
        "hrs_accession",
        "hri_accession"
    ]
)


# ============================================================
# 3. Fetch PUBLIC individualDetail pages
# ============================================================

print("===== 3. FETCH PUBLIC HRI INDIVIDUAL PAGES =====")
print()

hri_audit = []
fetch_failures = []

for i, hri in enumerate(
    unique_hri,
    1
):

    patient = next(
        iter(
            hri_to_patient[hri]
        )
    )

    url = (
        "https://ngdc.cncb.ac.cn/"
        "gsa-human/browse/"
        "individualDetail/"
        + hri
    )

    outfile = os.path.join(
        SOURCE_DIR,
        hri + ".html"
    )

    print(
        f"[{i:02d}/{len(unique_hri)}] "
        f"{patient:<4s} {hri}",
        end=" ... ",
        flush=True
    )

    ok, result = fetch_public_page(
        url,
        outfile
    )

    if not ok:

        print("FAILED")

        fetch_failures.append({

            "patient_id":
                patient,

            "hri_accession":
                hri,

            "url":
                url,

            "error":
                result
        })

        continue

    text = strip_html(
        result
    )

    low = text.lower()

    # --------------------------------------------------------
    # Explicit pathology terms only.
    #
    # "esophageal carcinoma" alone is NOT sufficient to assign
    # ESCC/EASC.
    # --------------------------------------------------------

    easc_hit = bool(
        re.search(
            r"\bEASC\b",
            text,
            flags=re.I
        )
        or
        re.search(
            r"adenosquamous",
            text,
            flags=re.I
        )
        or
        re.search(
            r"adeno[\s-]*squamous",
            text,
            flags=re.I
        )
    )

    adenocarcinoma_hit = bool(
        re.search(
            r"adenocarcinoma",
            text,
            flags=re.I
        )
    )

    escc_explicit_hit = bool(
        re.search(
            r"\bESCC\b",
            text,
            flags=re.I
        )
        or
        re.search(
            r"esophageal\s+squamous"
            r"\s+cell\s+carcinoma",
            text,
            flags=re.I
        )
    )

    histology_word_hit = bool(
        re.search(
            r"histolog|"
            r"patholog|"
            r"diagnos|"
            r"disease",
            text,
            flags=re.I
        )
    )

    # --------------------------------------------------------
    # Context snippets
    # --------------------------------------------------------

    snippets = []

    terms = [
        "EASC",
        "adenosquamous",
        "adeno-squamous",
        "adenocarcinoma",
        "ESCC",
        "squamous cell carcinoma",
        "histology",
        "pathology",
        "diagnosis",
        "disease"
    ]

    for term in terms:

        for match in re.finditer(
            re.escape(term),
            text,
            flags=re.I
        ):

            start = max(
                0,
                match.start() - 180
            )

            end = min(
                len(text),
                match.end() + 260
            )

            snippet = text[
                start:end
            ].strip()

            if snippet not in snippets:

                snippets.append(
                    snippet
                )

            if len(snippets) >= 10:
                break

        if len(snippets) >= 10:
            break

    hri_audit.append({

        "patient_id":
            patient,

        "hri_accession":
            hri,

        "public_page_url":
            url,

        "page_bytes":
            os.path.getsize(
                outfile
            ),

        "explicit_EASC_or_adenosquamous":
            "YES"
            if easc_hit
            else "NO",

        "explicit_adenocarcinoma":
            "YES"
            if adenocarcinoma_hit
            else "NO",

        "explicit_ESCC":
            "YES"
            if escc_explicit_hit
            else "NO",

        "histology_pathology_diagnosis_words":
            "YES"
            if histology_word_hit
            else "NO",

        "context_snippets":
            " || ".join(
                snippets
            ),

        "fetch_status":
            "SUCCESS"
    })

    print(
        "OK",
        (
            "  <-- EASC TERM"
            if easc_hit
            else ""
        )
    )


print()


# ============================================================
# 4. Require complete public-page retrieval
# ============================================================

if fetch_failures:

    failures_out = os.path.join(
        META_DIR,
        "STEP22A_04B_HRI_FETCH_FAILURES.csv"
    )

    write_csv(
        failures_out,
        fetch_failures,
        [
            "patient_id",
            "hri_accession",
            "url",
            "error"
        ]
    )

    print(
        "ERROR:",
        len(fetch_failures),
        "public HRI pages failed to fetch."
    )

    print()

    print(
        "Saved:"
    )

    print(
        failures_out
    )

    print()

    print(
        "PUBLIC_HRI_AUDIT_INCOMPLETE"
    )

    sys.exit(1)


audit_out = os.path.join(
    META_DIR,
    "STEP22A_04B_HRI_PUBLIC_METADATA_AUDIT.csv"
)

write_csv(
    audit_out,
    hri_audit,
    [
        "patient_id",
        "hri_accession",
        "public_page_url",
        "page_bytes",
        "explicit_EASC_or_adenosquamous",
        "explicit_adenocarcinoma",
        "explicit_ESCC",
        "histology_pathology_diagnosis_words",
        "context_snippets",
        "fetch_status"
    ]
)


# ============================================================
# 5. Identify explicit EASC candidates
# ============================================================

print("===== 4. EXPLICIT PUBLIC HISTOLOGY RESULT =====")
print()

easc_rows = [

    x for x in hri_audit

    if (
        x[
            "explicit_EASC_or_adenosquamous"
        ]
        == "YES"
    )
]

easc_patients = sorted(
    set(
        x["patient_id"]
        for x in easc_rows
    ),
    key=patient_key
)

print(
    "HRI pages with explicit EASC/"
    "adenosquamous:",
    len(easc_rows)
)

print(
    "Patient candidates:",
    (
        ", ".join(
            easc_patients
        )
        if easc_patients
        else "<none>"
    )
)

print()


# ============================================================
# 6. Conservative decision rule
# ============================================================

if len(easc_patients) == 1:

    decision = (
        "UNIQUE_PUBLIC_EASC_PATIENT_IDENTIFIED"
    )

    easc_patient = (
        easc_patients[0]
    )

elif len(easc_patients) == 0:

    decision = (
        "EASC_NOT_RESOLVED_ON_PUBLIC_HRI_PAGES"
    )

    easc_patient = ""

else:

    decision = (
        "EASC_PUBLIC_HRI_MAPPING_AMBIGUOUS"
    )

    easc_patient = ""


# ============================================================
# 7. Patient-level status table
# ============================================================

patient_status_rows = []

for patient in patients:

    relevant = [

        x for x in hri_audit

        if x["patient_id"]
        == patient
    ]

    explicit_easc = any(

        x[
            "explicit_EASC_or_adenosquamous"
        ]
        == "YES"

        for x in relevant
    )

    if (
        decision
        ==
        "UNIQUE_PUBLIC_EASC_PATIENT_IDENTIFIED"
        and
        patient == easc_patient
    ):

        status = (
            "PUBLICLY_IDENTIFIED_EASC"
        )

    elif (
        decision
        ==
        "UNIQUE_PUBLIC_EASC_PATIENT_IDENTIFIED"
    ):

        status = (
            "NOT_THE_PUBLICLY_IDENTIFIED_EASC_PATIENT"
        )

    else:

        status = (
            "HISTOLOGY_NOT_INDIVIDUALLY_RESOLVED"
        )

    patient_status_rows.append({

        "patient_id":
            patient,

        "HRI_accessions":
            ";".join(
                sorted(
                    patient_to_hri[
                        patient
                    ]
                )
            ),

        "explicit_EASC_public_HRI":
            "YES"
            if explicit_easc
            else "NO",

        "histology_mapping_status":
            status,

        "excluded_now":
            "NO"
    })


patient_status_out = os.path.join(
    META_DIR,
    "STEP22A_04B_PATIENT_HISTOLOGY_PUBLIC_STATUS.csv"
)

write_csv(
    patient_status_out,
    patient_status_rows,
    [
        "patient_id",
        "HRI_accessions",
        "explicit_EASC_public_HRI",
        "histology_mapping_status",
        "excluded_now"
    ]
)


# ============================================================
# 8. Make explicit decision file
# ============================================================

decision_out = os.path.join(
    META_DIR,
    "STEP22A_04B_HISTOLOGY_MAPPING_DECISION.txt"
)

with open(
    decision_out,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-4B PUBLIC HRI HISTOLOGY AUDIT\n"
    )

    f.write(
        "=" * 65
        + "\n\n"
    )

    f.write(
        "Source cohort paper reports one EASC patient: YES\n"
    )

    f.write(
        f"Patients audited: {len(patients)}\n"
    )

    f.write(
        f"Unique HRI pages audited: {len(unique_hri)}\n"
    )

    f.write(
        "Public HRI fetch failures: 0\n"
    )

    f.write(
        "Public HRI pages containing explicit "
        f"EASC/adenosquamous: {len(easc_rows)}\n"
    )

    f.write(
        "Candidate patients: "
        + (
            ",".join(
                easc_patients
            )
            if easc_patients
            else "NONE"
        )
        + "\n"
    )

    f.write(
        f"Decision: {decision}\n"
    )

    f.write(
        "EASC patient: "
        + (
            easc_patient
            if easc_patient
            else "UNRESOLVED"
        )
        + "\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Samples excluded in Step22A-4B: 0\n"
    )

    f.write(
        "CAF annotation performed: NO\n"
    )

    f.write(
        "CORE2 validation performed: NO\n"
    )

    f.write(
        "Frozen Step19-21 modified: NO\n"
    )

    f.write(
        "Manuscript modified: NO\n"
    )


# ============================================================
# 9. Confirm no extra OMIX TAR download
# ============================================================

print(
    "===== 5. OMIX DOWNLOAD SAFETY CHECK ====="
)
print()

raw_root = os.path.join(
    PROJECT,
    "01_raw",
    "OMIX005710"
)

tar_files = []

for root, dirs, filenames in os.walk(
    raw_root
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


print(
    "Current OMIX TAR files:",
    len(tar_files)
)

print(
    "Expected from Step22A-3A:",
    12
)

if len(tar_files) != 12:

    print()

    print(
        "WARNING: TAR count differs from "
        "the 12 paired files expected."
    )

else:

    print(
        "✓ No additional OMIX TAR files detected."
    )

print()


# ============================================================
# 10. Frozen project check
# ============================================================

print(
    "===== 6. FROZEN PROJECT CHECK ====="
)
print()

frozen_changes = []

if os.path.exists(
    FROZEN_INV
):

    frozen_rows = read_csv(
        FROZEN_INV
    )

    for row in frozen_rows:

        path = (
            row.get("path")
            or ""
        )

        if not path:
            continue

        if not os.path.exists(
            path
        ):

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
            os.path.getsize(
                path
            )
        )

        if (
            expected_size
            and
            expected_size
            != actual_size
        ):

            frozen_changes.append(
                (
                    path,
                    "SIZE_CHANGED"
                )
            )


print(
    "Frozen original files changed:",
    len(frozen_changes)
)

if frozen_changes:

    for x in frozen_changes[:20]:
        print(x)

    print()

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

print("=" * 78)
print("STEP22A-4B COMPLETE")
print("=" * 78)
print()

print(
    "PATIENTS AUDITED             :",
    len(patients)
)

print(
    "UNIQUE PUBLIC HRI PAGES      :",
    len(unique_hri)
)

print(
    "HRI FETCH FAILURES           : 0"
)

print()

print(
    "PAPER REPORTS ONE EASC       : YES"
)

print(
    "EXPLICIT EASC HRI PAGES      :",
    len(easc_rows)
)

print(
    "EASC PATIENT CANDIDATES      :",
    (
        ", ".join(
            easc_patients
        )
        if easc_patients
        else "<none>"
    )
)

print()

print(
    "PUBLIC HISTOLOGY DECISION    :",
    decision
)

print(
    "EASC PATIENT                 :",
    (
        easc_patient
        if easc_patient
        else "UNRESOLVED"
    )
)

print()

print(
    "SAMPLES EXCLUDED NOW         : 0"
)

print(
    "CAF ANNOTATION PERFORMED     : NO"
)

print(
    "CORE2 VALIDATION PERFORMED   : NO"
)

print(
    "NEW OMIX TAR DOWNLOADED      : 0"
)

print()

print(
    "FROZEN STEP19-21 MODIFIED    : NO"
)

print(
    "MANUSCRIPT MODIFIED          : NO"
)

print()

print("OUTPUTS:")

for path in [
    mapping_out,
    audit_out,
    patient_status_out,
    decision_out
]:

    print(
        " ",
        os.path.relpath(
            path,
            PROJECT
        )
    )

print()

if (
    decision
    ==
    "UNIQUE_PUBLIC_EASC_PATIENT_IDENTIFIED"
):

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "STEP22A-4C = EXCLUDE THE IDENTIFIED EASC "
        "PATIENT + REBUILD ESCC-ONLY PAIRED SET"
    )

else:

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "PUBLIC HISTOLOGY MAPPING EXHAUSTED"
    )

    print(
        "DO NOT START PRIMARY CORE2 VALIDATION."
    )

    print(
        "Do not guess the EASC patient."
    )

print()

print("STOP HERE.")
print()
