#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# Step22A-1
#
# Official OMIX005710 / PRJCA016745 metadata inventory
#
# IMPORTANT:
# - Fetch metadata pages ONLY
# - DO NOT download any .tar data files
# - DO NOT modify frozen Step19/20/21 outputs
# - DO NOT modify manuscript
# ============================================================

import csv
import hashlib
import html as html_lib
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from urllib.parse import urljoin


# ============================================================
# 0. Paths / sources
# ============================================================

PROJECT = os.getcwd()

META_DIR = os.path.join(
    PROJECT,
    "00_metadata",
    "OMIX005710"
)

SOURCE_DIR = os.path.join(
    META_DIR,
    "STEP22A_01_source_html"
)

LOG_DIR = os.path.join(
    PROJECT,
    "08_logs"
)

os.makedirs(META_DIR, exist_ok=True)
os.makedirs(SOURCE_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

OMIX_URL = (
    "https://ngdc.cncb.ac.cn/"
    "omix/release/OMIX005710"
)

BIOPROJECT_URL = (
    "https://ngdc.cncb.ac.cn/"
    "bioproject/browse/PRJCA016745"
)

SAMPLE_URL_PREFIX = (
    "https://ngdc.cncb.ac.cn/"
    "gsa-human/browse/sampleDetail/"
)

EXPECTED_FILES = 46
EXPECTED_BIOSAMPLES = 46


print()
print("=" * 72)
print("STEP22A-1: OFFICIAL OMIX005710 METADATA INVENTORY")
print("=" * 72)
print()

print("PROJECT:")
print(PROJECT)
print()

print("IMPORTANT:")
print("  Metadata pages only.")
print("  No .tar data files will be downloaded.")
print()


# ============================================================
# 1. Helpers
# ============================================================

def fetch_page(url, outfile, retries=4):
    """
    Download a small HTML metadata page using curl.
    Large-data host download.cncb.ac.cn is explicitly blocked.
    """

    if "download.cncb.ac.cn" in url:
        raise RuntimeError(
            "GUARDRAIL: attempted large-data download URL: "
            + url
        )

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
        "--retry", str(retries),
        "--retry-delay", "2",
        "--connect-timeout", "20",
        "--max-time", "90",
        "-A",
        "Mozilla/5.0 Step22A-metadata-audit",
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

        raise RuntimeError(
            "curl failed for "
            + url
            + "\n"
            + result.stderr.strip()
        )

    if (
        not os.path.exists(tmp)
        or os.path.getsize(tmp) < 100
    ):
        raise RuntimeError(
            "Downloaded metadata page is unexpectedly small: "
            + url
        )

    os.replace(tmp, outfile)

    with open(
        outfile,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:
        return f.read()


def strip_html(x):
    x = re.sub(
        r"<script\b.*?</script>",
        " ",
        x,
        flags=re.I | re.S
    )

    x = re.sub(
        r"<style\b.*?</style>",
        " ",
        x,
        flags=re.I | re.S
    )

    x = re.sub(
        r"<br\s*/?>",
        " ",
        x,
        flags=re.I
    )

    x = re.sub(
        r"<[^>]+>",
        " ",
        x
    )

    x = html_lib.unescape(x)
    x = x.replace("\xa0", " ")

    return re.sub(
        r"\s+",
        " ",
        x
    ).strip()


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


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


def rex(pattern, text, default=""):
    m = re.search(
        pattern,
        text,
        flags=re.I | re.S
    )

    if not m:
        return default

    return m.group(1).strip()


def size_to_mb(value, unit):
    value = float(value)
    unit = unit.upper()

    if unit == "KB":
        return value / 1024

    if unit == "MB":
        return value

    if unit == "GB":
        return value * 1024

    if unit == "TB":
        return value * 1024 * 1024

    return None


# ============================================================
# 2. Fetch official OMIX page
# ============================================================

print("===== 1. FETCH OFFICIAL OMIX PAGE =====")
print()

omix_html_path = os.path.join(
    SOURCE_DIR,
    "OMIX005710.html"
)

try:
    omix_html = fetch_page(
        OMIX_URL,
        omix_html_path
    )

except Exception as e:
    print("ERROR:")
    print(e)
    print()
    print("STEP22A_OMIX_PAGE_FETCH_FAILED")
    sys.exit(1)

print("✓ OMIX page saved:")
print(omix_html_path)
print()


# ============================================================
# 3. Parse OMIX file inventory
# ============================================================

print("===== 2. PARSE OMIX FILE INVENTORY =====")
print()

rows_html = re.findall(
    r"<tr\b[^>]*>(.*?)</tr>",
    omix_html,
    flags=re.I | re.S
)

file_rows = []

for row_html in rows_html:

    row_text = strip_html(row_html)

    file_match = re.search(
        r"\b(OMIX005710-\d+)\b",
        row_text
    )

    hrs_match = re.search(
        r"\b(HRS\d+)\b",
        row_text
    )

    size_match = re.search(
        r"\b(\d+(?:\.\d+)?)\s*"
        r"(KB|MB|GB|TB)\b",
        row_text,
        flags=re.I
    )

    if not (
        file_match
        and hrs_match
        and size_match
    ):
        continue

    file_id = file_match.group(1)
    hrs = hrs_match.group(1)

    size_value = float(
        size_match.group(1)
    )

    size_unit = size_match.group(2).upper()

    size_mb = size_to_mb(
        size_value,
        size_unit
    )

    # Extract official href if exposed in HTML.
    hrefs = re.findall(
        r'href=["\']([^"\']+)["\']',
        row_html,
        flags=re.I
    )

    official_https = ""

    for href in hrefs:

        full = urljoin(
            OMIX_URL,
            html_lib.unescape(href)
        )

        if (
            "download.cncb.ac.cn" in full
            and file_id in full
        ):
            official_https = full
            break

    # Official OMIX download path.
    # This is RECORDED ONLY.
    # It is NOT accessed in Step22A-1.
    if not official_https:
        official_https = (
            "https://download.cncb.ac.cn/"
            "OMIX/OMIX005710/"
            + file_id
            + ".tar"
        )

    file_rows.append({
        "omix_accession": "OMIX005710",
        "bioproject": "PRJCA016745",
        "file_id": file_id,
        "hrs_accession": hrs,
        "file_suffix": "tar",
        "listed_size_value": size_value,
        "listed_size_unit": size_unit,
        "listed_size_MB": round(size_mb, 3),
        "download_url_https": official_https,
        "downloaded_now": "NO",
        "source_page": OMIX_URL
    })


# Remove possible duplicate parsing
dedup = {}

for row in file_rows:
    dedup[row["file_id"]] = row

file_rows = list(dedup.values())


def file_number(x):
    m = re.search(
        r"-(\d+)$",
        x["file_id"]
    )

    if m:
        return int(m.group(1))

    return 999999


file_rows.sort(
    key=file_number
)

print("OMIX file records found:", len(file_rows))

if len(file_rows) != EXPECTED_FILES:
    print()
    print(
        "ERROR: expected",
        EXPECTED_FILES,
        "official OMIX files but parsed",
        len(file_rows)
    )
    print()
    print("OMIX_FILE_INVENTORY_INCOMPLETE")
    sys.exit(1)

hrs_list = [
    x["hrs_accession"]
    for x in file_rows
]

if len(set(hrs_list)) != EXPECTED_FILES:
    print()
    print(
        "ERROR: HRS accessions are not uniquely mapped "
        "1:1 to OMIX files."
    )
    print("OMIX_HRS_MAPPING_NOT_UNIQUE")
    sys.exit(1)

total_mb = sum(
    x["listed_size_MB"]
    for x in file_rows
)

print("Unique HRS accessions :", len(set(hrs_list)))
print("Listed total size MB :", round(total_mb, 2))
print("Listed total size GiB:", round(total_mb / 1024, 3))
print()

file_inventory_csv = os.path.join(
    META_DIR,
    "STEP22A_OMIX_FILE_INVENTORY.csv"
)

write_csv(
    file_inventory_csv,
    file_rows,
    [
        "omix_accession",
        "bioproject",
        "file_id",
        "hrs_accession",
        "file_suffix",
        "listed_size_value",
        "listed_size_unit",
        "listed_size_MB",
        "download_url_https",
        "downloaded_now",
        "source_page"
    ]
)

print("✓ Saved:")
print(file_inventory_csv)
print()


# ============================================================
# 4. Fetch official BioProject page
# ============================================================

print("===== 3. FETCH OFFICIAL BIOPROJECT PAGE =====")
print()

bioproject_html_path = os.path.join(
    SOURCE_DIR,
    "PRJCA016745.html"
)

try:
    bioproject_html = fetch_page(
        BIOPROJECT_URL,
        bioproject_html_path
    )

except Exception as e:
    print("ERROR:")
    print(e)
    print()
    print("STEP22A_BIOPROJECT_FETCH_FAILED")
    sys.exit(1)

print("✓ BioProject page saved:")
print(bioproject_html_path)
print()

official_samc = sorted(
    set(
        re.findall(
            r"\bSAMC\d+\b",
            strip_html(bioproject_html)
        )
    )
)

print(
    "BioSample accessions detected on BioProject page:",
    len(official_samc)
)

if len(official_samc) != EXPECTED_BIOSAMPLES:
    print()
    print(
        "WARNING: expected 46 BioSample accessions, "
        "but detected",
        len(official_samc)
    )
    print(
        "The HRS sample-detail audit will still continue."
    )

print()


# ============================================================
# 5. Fetch official HRS sample-detail metadata
# ============================================================

print("===== 4. FETCH 46 OFFICIAL HRS SAMPLE METADATA PAGES =====")
print()

sample_rows = []
failed_hrs = []

for i, file_row in enumerate(file_rows, 1):

    hrs = file_row["hrs_accession"]

    url = SAMPLE_URL_PREFIX + hrs

    outfile = os.path.join(
        SOURCE_DIR,
        hrs + ".html"
    )

    print(
        f"[{i:02d}/{EXPECTED_FILES}] {hrs}",
        end=" ... ",
        flush=True
    )

    success = False
    raw = ""
    error_message = ""

    for attempt in range(1, 4):

        try:
            raw = fetch_page(
                url,
                outfile
            )

            success = True
            break

        except Exception as e:
            error_message = str(e)

            if attempt < 3:
                time.sleep(2)

    if not success:
        print("FAILED")

        failed_hrs.append(
            {
                "hrs_accession": hrs,
                "url": url,
                "error": error_message
            }
        )

        continue

    text = strip_html(raw)

    page_hrs = rex(
        r"\bAccession\s+(HRS\d+)\b",
        text
    )

    biosample = rex(
        r"BioSample\s+Accession\s+"
        r"(SAMC\d+)",
        text
    )

    individual = rex(
        r"Individual\s+accession\s+"
        r"(HRI\d+)",
        text
    )

    study = rex(
        r"Study\s+accession\s+"
        r"(HRA\d+)",
        text
    )

    sample_name = rex(
        r"Sample\s+name\s+"
        r"([A-Za-z0-9_.-]+)",
        text
    )

    sample_title = rex(
        r"Sample\s+title\s+(.+?)"
        r"\s+Collection\s+date\b",
        text
    )

    collection_date = rex(
        r"Collection\s+date\s+"
        r"([0-9]{4}-[0-9]{2}-[0-9]{2})",
        text
    )

    sample_type = rex(
        r"Sample\s+Type\s+(.+?)"
        r"\s+Description\b",
        text
    )

    description = rex(
        r"Description\s+(.+?)"
        r"\s+Tissue\b",
        text
    )

    tissue = rex(
        r"Tissue\s+(.+?)"
        r"\s+Age\b",
        text
    )

    age = rex(
        r"\bAge\s+([0-9]+(?:\.[0-9]+)?)\b",
        text
    )

    release_date = rex(
        r"Release\s+date\s+"
        r"([0-9]{4}-[0-9]{2}-[0-9]{2})",
        text
    )

    # --------------------------------------------------------
    # Sample-name coding audit
    # DO NOT treat this as final biological truth yet.
    # --------------------------------------------------------

    name_patient = ""
    name_tissue = ""
    name_timepoint = ""

    m = re.fullmatch(
        r"P(\d+)_(T|N)_(A|B)",
        sample_name,
        flags=re.I
    )

    if m:
        name_patient = "P" + str(int(m.group(1)))

        name_tissue = (
            "Tumor"
            if m.group(2).upper() == "T"
            else "Adjacent_normal"
        )

        name_timepoint = (
            "After"
            if m.group(3).upper() == "A"
            else "Before"
        )

    # --------------------------------------------------------
    # Description-derived coding
    # --------------------------------------------------------

    desc_patient = ""
    desc_tissue = ""
    desc_timepoint = ""

    m_patient = re.search(
        r"Patient\s*(\d+)",
        description,
        flags=re.I
    )

    if m_patient:
        desc_patient = (
            "P"
            + str(int(m_patient.group(1)))
        )

    dlow = description.lower()

    if (
        "normal adjacent" in dlow
        or "adjacent" in dlow
    ):
        desc_tissue = "Adjacent_normal"

    elif "tumor" in dlow:
        desc_tissue = "Tumor"

    if "before treatment" in dlow:
        desc_timepoint = "Before"

    elif "after treatment" in dlow:
        desc_timepoint = "After"

    # --------------------------------------------------------
    # Sample-title-derived tissue
    # --------------------------------------------------------

    title_tissue = ""

    tlow = sample_title.lower()

    if "adjacent" in tlow:
        title_tissue = "Adjacent_normal"

    elif (
        "esophageal carcinoma" in tlow
        or "oesophageal carcinoma" in tlow
    ):
        title_tissue = "Tumor"

    # --------------------------------------------------------
    # Conflict audit
    # --------------------------------------------------------

    patient_values = set(
        x for x in [
            name_patient,
            desc_patient
        ]
        if x
    )

    tissue_values = set(
        x for x in [
            name_tissue,
            desc_tissue,
            title_tissue
        ]
        if x
    )

    timepoint_values = set(
        x for x in [
            name_timepoint,
            desc_timepoint
        ]
        if x
    )

    patient_conflict = (
        "YES"
        if len(patient_values) > 1
        else "NO"
    )

    tissue_conflict = (
        "YES"
        if len(tissue_values) > 1
        else "NO"
    )

    timepoint_conflict = (
        "YES"
        if len(timepoint_values) > 1
        else "NO"
    )

    sample_rows.append({
        "omix_file_id": file_row["file_id"],
        "hrs_accession": hrs,
        "hrs_accession_page": page_hrs,
        "biosample_accession": biosample,
        "individual_accession": individual,
        "study_accession": study,
        "sample_name": sample_name,
        "sample_title": sample_title,
        "collection_date": collection_date,
        "sample_type": sample_type,
        "description": description,
        "tissue_field": tissue,
        "age": age,
        "release_date": release_date,

        "patient_from_name": name_patient,
        "patient_from_description": desc_patient,

        "tissue_from_name": name_tissue,
        "tissue_from_description": desc_tissue,
        "tissue_from_title": title_tissue,

        "timepoint_from_name": name_timepoint,
        "timepoint_from_description": desc_timepoint,

        "patient_conflict": patient_conflict,
        "tissue_conflict": tissue_conflict,
        "timepoint_conflict": timepoint_conflict,

        "histology_status": "NOT_AUDITED_YET",
        "canonical_inclusion_status": "NOT_ASSIGNED_YET",

        "metadata_source_url": url
    })

    print("OK")

    # Gentle request pacing
    time.sleep(0.15)


print()

if failed_hrs:

    failed_csv = os.path.join(
        META_DIR,
        "STEP22A_SAMPLE_METADATA_FETCH_FAILURES.csv"
    )

    write_csv(
        failed_csv,
        failed_hrs,
        [
            "hrs_accession",
            "url",
            "error"
        ]
    )

    print(
        "ERROR:",
        len(failed_hrs),
        "HRS metadata pages failed."
    )
    print()
    print("Saved failures:")
    print(failed_csv)
    print()
    print("STEP22A_SAMPLE_METADATA_INCOMPLETE")
    sys.exit(1)


# ============================================================
# 6. Validate sample metadata
# ============================================================

print("===== 5. SAMPLE METADATA VALIDATION =====")
print()

print("HRS metadata pages parsed :", len(sample_rows))

if len(sample_rows) != EXPECTED_FILES:
    print(
        "ERROR: expected 46 parsed HRS metadata records."
    )
    print("STEP22A_SAMPLE_METADATA_INCOMPLETE")
    sys.exit(1)

page_hrs_values = [
    x["hrs_accession_page"]
    for x in sample_rows
]

bad_hrs_page = [
    x
    for x in sample_rows
    if x["hrs_accession_page"]
    != x["hrs_accession"]
]

if bad_hrs_page:
    print(
        "ERROR: HRS accession mismatch detected "
        "between requested and returned pages."
    )
    print("STEP22A_HRS_PAGE_MISMATCH")
    sys.exit(1)

biosamples_mapped = [
    x["biosample_accession"]
    for x in sample_rows
    if x["biosample_accession"]
]

print(
    "BioSample accessions mapped:",
    len(set(biosamples_mapped))
)

if len(set(biosamples_mapped)) != EXPECTED_BIOSAMPLES:
    print()
    print(
        "ERROR: expected 46 unique BioSample mappings."
    )
    print("STEP22A_BIOSAMPLE_MAPPING_INCOMPLETE")
    sys.exit(1)

if official_samc:

    mapped_set = set(biosamples_mapped)
    official_set = set(official_samc)

    only_project = sorted(
        official_set - mapped_set
    )

    only_hrs = sorted(
        mapped_set - official_set
    )

    print(
        "BioProject ↔ HRS BioSample set match:",
        "YES"
        if not only_project and not only_hrs
        else "NO"
    )

    if only_project or only_hrs:
        print()
        print(
            "ERROR: BioSample sets do not match."
        )
        print(
            "Only in BioProject:",
            only_project
        )
        print(
            "Only in HRS mapping:",
            only_hrs
        )
        print("STEP22A_BIOSAMPLE_SET_MISMATCH")
        sys.exit(1)

print()


# ============================================================
# 7. Save RAW official metadata
# ============================================================

raw_metadata_csv = os.path.join(
    META_DIR,
    "STEP22A_SAMPLE_METADATA_RAW.csv"
)

sample_fields = [
    "omix_file_id",
    "hrs_accession",
    "hrs_accession_page",
    "biosample_accession",
    "individual_accession",
    "study_accession",
    "sample_name",
    "sample_title",
    "collection_date",
    "sample_type",
    "description",
    "tissue_field",
    "age",
    "release_date",

    "patient_from_name",
    "patient_from_description",

    "tissue_from_name",
    "tissue_from_description",
    "tissue_from_title",

    "timepoint_from_name",
    "timepoint_from_description",

    "patient_conflict",
    "tissue_conflict",
    "timepoint_conflict",

    "histology_status",
    "canonical_inclusion_status",

    "metadata_source_url"
]

write_csv(
    raw_metadata_csv,
    sample_rows,
    sample_fields
)

print("✓ Saved:")
print(raw_metadata_csv)
print()


# ============================================================
# 8. Create coding-conflict audit
# ============================================================

audit_rows = []

for x in sample_rows:

    any_conflict = (
        "YES"
        if (
            x["patient_conflict"] == "YES"
            or x["tissue_conflict"] == "YES"
            or x["timepoint_conflict"] == "YES"
        )
        else "NO"
    )

    audit_rows.append({
        "hrs_accession": x["hrs_accession"],
        "biosample_accession": x["biosample_accession"],
        "sample_name": x["sample_name"],
        "sample_title": x["sample_title"],
        "description": x["description"],

        "patient_from_name":
            x["patient_from_name"],

        "patient_from_description":
            x["patient_from_description"],

        "tissue_from_name":
            x["tissue_from_name"],

        "tissue_from_description":
            x["tissue_from_description"],

        "tissue_from_title":
            x["tissue_from_title"],

        "timepoint_from_name":
            x["timepoint_from_name"],

        "timepoint_from_description":
            x["timepoint_from_description"],

        "patient_conflict":
            x["patient_conflict"],

        "tissue_conflict":
            x["tissue_conflict"],

        "timepoint_conflict":
            x["timepoint_conflict"],

        "any_conflict": any_conflict,

        "review_status":
            "MANUAL_BIOLOGICAL_REVIEW_REQUIRED"
            if any_conflict == "YES"
            else "NO_STRING_CONFLICT"
    })


audit_csv = os.path.join(
    META_DIR,
    "STEP22A_SAMPLE_NAME_DESCRIPTION_AUDIT.csv"
)

write_csv(
    audit_csv,
    audit_rows,
    [
        "hrs_accession",
        "biosample_accession",
        "sample_name",
        "sample_title",
        "description",

        "patient_from_name",
        "patient_from_description",

        "tissue_from_name",
        "tissue_from_description",
        "tissue_from_title",

        "timepoint_from_name",
        "timepoint_from_description",

        "patient_conflict",
        "tissue_conflict",
        "timepoint_conflict",

        "any_conflict",
        "review_status"
    ]
)

print("✓ Saved:")
print(audit_csv)
print()


# ============================================================
# 9. Source provenance manifest
# ============================================================

source_rows = []

for filename in sorted(
    os.listdir(SOURCE_DIR)
):

    path = os.path.join(
        SOURCE_DIR,
        filename
    )

    if not os.path.isfile(path):
        continue

    source_rows.append({
        "file": filename,
        "size_bytes": os.path.getsize(path),
        "sha256": sha256_file(path),
        "captured_at": datetime.now().isoformat(
            timespec="seconds"
        )
    })


source_manifest_csv = os.path.join(
    META_DIR,
    "STEP22A_METADATA_SOURCE_MANIFEST.csv"
)

write_csv(
    source_manifest_csv,
    source_rows,
    [
        "file",
        "size_bytes",
        "sha256",
        "captured_at"
    ]
)

print("✓ Saved:")
print(source_manifest_csv)
print()


# ============================================================
# 10. Summary audit
# ============================================================

patients_name = sorted(
    set(
        x["patient_from_name"]
        for x in sample_rows
        if x["patient_from_name"]
    ),
    key=lambda z: int(
        re.sub(r"\D", "", z)
    )
)

tissue_name_counts = {
    "Tumor": 0,
    "Adjacent_normal": 0,
    "Unknown": 0
}

timepoint_name_counts = {
    "Before": 0,
    "After": 0,
    "Unknown": 0
}

for x in sample_rows:

    tissue_value = x["tissue_from_name"]

    if tissue_value in tissue_name_counts:
        tissue_name_counts[tissue_value] += 1
    else:
        tissue_name_counts["Unknown"] += 1

    tp = x["timepoint_from_name"]

    if tp in timepoint_name_counts:
        timepoint_name_counts[tp] += 1
    else:
        timepoint_name_counts["Unknown"] += 1


conflict_rows = [
    x
    for x in audit_rows
    if x["any_conflict"] == "YES"
]

summary_txt = os.path.join(
    META_DIR,
    "STEP22A_01_OFFICIAL_METADATA_SUMMARY.txt"
)

with open(
    summary_txt,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-1 OFFICIAL METADATA SUMMARY\n"
    )
    f.write("=" * 60 + "\n\n")

    f.write(
        f"OMIX accession: OMIX005710\n"
    )

    f.write(
        f"BioProject: PRJCA016745\n"
    )

    f.write(
        f"OMIX files: {len(file_rows)}\n"
    )

    f.write(
        f"Unique HRS: {len(set(hrs_list))}\n"
    )

    f.write(
        "Unique BioSamples: "
        f"{len(set(biosamples_mapped))}\n"
    )

    f.write(
        "Patients derived from SAMPLE NAME ONLY: "
        f"{len(patients_name)}\n"
    )

    f.write(
        "Tumor-coded sample names: "
        f"{tissue_name_counts['Tumor']}\n"
    )

    f.write(
        "Adjacent-normal-coded sample names: "
        f"{tissue_name_counts['Adjacent_normal']}\n"
    )

    f.write(
        "Before-coded sample names: "
        f"{timepoint_name_counts['Before']}\n"
    )

    f.write(
        "After-coded sample names: "
        f"{timepoint_name_counts['After']}\n"
    )

    f.write(
        f"String conflicts requiring review: "
        f"{len(conflict_rows)}\n"
    )

    f.write(
        f"Listed total size MB: "
        f"{round(total_mb, 2)}\n"
    )

    f.write(
        f"Listed total size GiB: "
        f"{round(total_mb / 1024, 3)}\n"
    )

    f.write("\n")
    f.write(
        "HISTOLOGY STATUS: NOT AUDITED YET\n"
    )

    f.write(
        "CANONICAL SAMPLE INCLUSION: "
        "NOT ASSIGNED YET\n"
    )

    f.write(
        "DATA FILES DOWNLOADED: NO\n"
    )


print("✓ Saved:")
print(summary_txt)
print()


# ============================================================
# 11. Print conflicts, if any
# ============================================================

print("===== 6. STRING-CONFLICT AUDIT =====")
print()

print(
    "Rows with patient/tissue/timepoint "
    "metadata conflict:",
    len(conflict_rows)
)

if conflict_rows:

    print()
    print(
        "These rows are NOT automatically corrected "
        "or excluded."
    )
    print()

    for x in conflict_rows:

        print("-" * 72)
        print("HRS        :", x["hrs_accession"])
        print("Sample name:", x["sample_name"])
        print("Title      :", x["sample_title"])
        print("Description:", x["description"])

        print(
            "Name tissue:",
            x["tissue_from_name"]
        )

        print(
            "Desc tissue:",
            x["tissue_from_description"]
        )

        print(
            "Title tissue:",
            x["tissue_from_title"]
        )

        print(
            "Name time  :",
            x["timepoint_from_name"]
        )

        print(
            "Desc time  :",
            x["timepoint_from_description"]
        )

        print(
            "Patient conflict  :",
            x["patient_conflict"]
        )

        print(
            "Tissue conflict   :",
            x["tissue_conflict"]
        )

        print(
            "Timepoint conflict:",
            x["timepoint_conflict"]
        )

    print("-" * 72)

print()


# ============================================================
# 12. Final safety checks
# ============================================================

print("===== 7. FINAL SAFETY CHECK =====")
print()

tar_files_local = []

raw_root = os.path.join(
    PROJECT,
    "01_raw",
    "OMIX005710"
)

for root, dirs, files in os.walk(raw_root):

    for filename in files:

        if filename.lower().endswith(".tar"):
            tar_files_local.append(
                os.path.join(
                    root,
                    filename
                )
            )

print(
    "OMIX .tar files downloaded in Step22A:",
    len(tar_files_local)
)

if tar_files_local:

    print()
    print(
        "ERROR: unexpected .tar files found."
    )

    for x in tar_files_local:
        print(x)

    print()
    print("STEP22A_UNEXPECTED_DATA_DOWNLOAD")
    sys.exit(1)


# ============================================================
# Final
# ============================================================

print()
print("=" * 72)
print("STEP22A-1 COMPLETE")
print("=" * 72)
print()

print("OFFICIAL OMIX FILE RECORDS :", len(file_rows))
print("OFFICIAL HRS RECORDS       :", len(set(hrs_list)))
print(
    "OFFICIAL BIOSAMPLE RECORDS :",
    len(set(biosamples_mapped))
)
print(
    "PATIENT IDs FROM NAMES     :",
    len(patients_name)
)

print()
print(
    "LISTED TOTAL DATA SIZE     : "
    f"{round(total_mb, 2)} MB "
    f"({round(total_mb / 1024, 3)} GiB)"
)

print()
print(
    "STRING CONFLICTS TO REVIEW :",
    len(conflict_rows)
)

print()
print("HISTOLOGY AUDITED          : NO")
print("CANONICAL INCLUSION SET    : NOT YET")
print("OMIX TAR FILES DOWNLOADED  : 0")
print()
print("FROZEN STEP19-21 MODIFIED  : NO")
print("MANUSCRIPT MODIFIED        : NO")

print()
print("OUTPUTS:")
print(
    "  00_metadata/OMIX005710/"
    "STEP22A_OMIX_FILE_INVENTORY.csv"
)

print(
    "  00_metadata/OMIX005710/"
    "STEP22A_SAMPLE_METADATA_RAW.csv"
)

print(
    "  00_metadata/OMIX005710/"
    "STEP22A_SAMPLE_NAME_DESCRIPTION_AUDIT.csv"
)

print(
    "  00_metadata/OMIX005710/"
    "STEP22A_METADATA_SOURCE_MANIFEST.csv"
)

print(
    "  00_metadata/OMIX005710/"
    "STEP22A_01_OFFICIAL_METADATA_SUMMARY.txt"
)

print()

if len(conflict_rows) == 0:
    print(
        "READY FOR STEP22A-2: "
        "HISTOLOGY + CANONICAL SAMPLE AUDIT"
    )
else:
    print(
        "READY FOR STEP22A-2: "
        "RESOLVE METADATA CONFLICTS + "
        "HISTOLOGY AUDIT"
    )

print()
print("STOP HERE.")
print("DO NOT DOWNLOAD OMIX TAR FILES YET.")
print()
