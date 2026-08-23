#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-4A
#
# AUTHOR SOURCE + HISTOLOGY AUDIT
#
# Purposes:
# 1. Audit official supplementary materials from the two
#    primary publications using PRJCA016745 / OMIX005710.
# 2. Identify the reported EASC patient if publicly resolvable.
# 3. Search for author-provided cell-level / cluster annotation.
# 4. Test exact barcode compatibility against our 12 paired
#    tumor samples when possible.
#
# IMPORTANT:
# - NO additional OMIX TAR download
# - NO Seurat
# - NO clustering
# - NO cell annotation
# - NO CORE2 scoring
# - NO changes to frozen Step19-21
# - NO manuscript changes
# ============================================================

import csv
import hashlib
import os
import re
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime


PROJECT = os.getcwd()

META_DIR = os.path.join(
    PROJECT,
    "00_metadata",
    "OMIX005710"
)

SOURCE_DIR = os.path.join(
    META_DIR,
    "STEP22A_04_author_sources"
)

MATRIX_AUDIT = os.path.join(
    META_DIR,
    "STEP22A_03B_MATRIX_AUDIT.csv"
)

FROZEN_INV = os.path.join(
    META_DIR,
    "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

os.makedirs(
    SOURCE_DIR,
    exist_ok=True
)


print()
print("=" * 78)
print("STEP22A-4A: AUTHOR SOURCE + HISTOLOGY AUDIT")
print("=" * 78)
print()


# ============================================================
# 1. Helpers
# ============================================================

def read_csv(path, delimiter=","):

    if not os.path.exists(path):
        raise SystemExit(
            "ERROR: missing required file:\n"
            + path
        )

    with open(
        path,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as f:

        return list(
            csv.DictReader(
                f,
                delimiter=delimiter
            )
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


def sha256_file(path):

    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            chunk = f.read(
                1024 * 1024
            )

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


def download_small_file(url, outfile):

    # Only Springer supplementary material is allowed here.
    allowed = (
        url.startswith(
            "https://media.springernature.com/"
        )
    )

    if not allowed:
        raise RuntimeError(
            "Blocked unexpected download host:\n"
            + url
        )

    if os.path.exists(outfile):

        size = os.path.getsize(outfile)

        if size > 100:
            return "EXISTING"

    tmp = outfile + ".part"

    if os.path.exists(tmp):
        os.remove(tmp)

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--silent",
        "--show-error",
        "--retry", "5",
        "--retry-delay", "3",
        "--connect-timeout", "20",
        "--max-time", "180",
        "-A",
        "Mozilla/5.0 Step22A-author-source-audit",
        url,
        "-o",
        tmp
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:

        if os.path.exists(tmp):
            os.remove(tmp)

        raise RuntimeError(
            result.stderr.strip()
        )

    if (
        not os.path.exists(tmp)
        or os.path.getsize(tmp) < 100
    ):
        raise RuntimeError(
            "Downloaded file unexpectedly small."
        )

    os.replace(
        tmp,
        outfile
    )

    return "DOWNLOADED"


def normalize_space(x):

    x = str(x)

    x = x.replace(
        "\u00a0",
        " "
    )

    x = re.sub(
        r"\s+",
        " ",
        x
    )

    return x.strip()


# ============================================================
# 2. Official supplementary files
# ============================================================

print("===== 1. OFFICIAL SUPPLEMENTARY SOURCE LIST =====")
print()

sources = []


# ------------------------------------------------------------
# Genome Medicine 2024
#
# DOI:
# 10.1186/s13073-024-01320-9
# ------------------------------------------------------------

gm_base = (
    "https://media.springernature.com/original/"
    "springer-static/esm/"
    "art%3A10.1186%2Fs13073-024-01320-9/"
    "MediaObjects/"
)

gm_files = [
    (
        1,
        "docx",
        "All supplementary figures"
    ),
    (
        2,
        "xlsx",
        "scRNA-seq quality control data"
    ),
    (
        3,
        "docx",
        "Table S1 clinical characteristics"
    ),
    (
        4,
        "xlsx",
        "Cell count data"
    ),
    (
        5,
        "xlsx",
        "Cancer-cell differential expression"
    ),
    (
        6,
        "xlsx",
        "CD8 effector T-cell differential expression"
    ),
    (
        7,
        "xlsx",
        "NK-cell differential expression"
    ),
    (
        8,
        "xlsx",
        "Treg differential expression"
    ),
    (
        9,
        "xlsx",
        "CellPhoneDB raw data"
    )
]

for n, ext, description in gm_files:

    filename = (
        f"13073_2024_1320_MOESM{n}_ESM.{ext}"
    )

    sources.append({

        "paper":
            "Genome_Medicine_2024",

        "doi":
            "10.1186/s13073-024-01320-9",

        "source_number":
            n,

        "description":
            description,

        "filename":
            filename,

        "url":
            gm_base + filename
    })


# ------------------------------------------------------------
# Biomarker Research 2024 CAF paper
#
# DOI:
# 10.1186/s40364-024-00656-z
# ------------------------------------------------------------

br_base = (
    "https://media.springernature.com/original/"
    "springer-static/esm/"
    "art%3A10.1186%2Fs40364-024-00656-z/"
    "MediaObjects/"
)

br_files = [
    (1, "xlsx"),
    (2, "xlsx"),
    (3, "xlsx"),
    (4, "csv"),
    (5, "xlsx"),
    (6, "xlsx"),
    (7, "xlsx"),
    (8, "docx")
]

for n, ext in br_files:

    filename = (
        f"40364_2024_656_MOESM{n}_ESM.{ext}"
    )

    sources.append({

        "paper":
            "Biomarker_Research_CAF_2024",

        "doi":
            "10.1186/s40364-024-00656-z",

        "source_number":
            n,

        "description":
            "Official electronic supplementary material",

        "filename":
            filename,

        "url":
            br_base + filename
    })


print(
    "Official supplementary files scheduled:",
    len(sources)
)

print(
    "  Genome Medicine :",
    len(gm_files)
)

print(
    "  Biomarker Res.  :",
    len(br_files)
)

print()


# ============================================================
# 3. Download supplementary materials
# ============================================================

print("===== 2. DOWNLOAD SMALL OFFICIAL SUPPLEMENTS =====")
print()

manifest = []
download_failures = []

for i, src in enumerate(
    sources,
    1
):

    outfile = os.path.join(
        SOURCE_DIR,
        src["filename"]
    )

    print(
        f"[{i:02d}/{len(sources)}] "
        f"{src['filename']}",
        end=" ... ",
        flush=True
    )

    try:

        status = download_small_file(
            src["url"],
            outfile
        )

        size = os.path.getsize(
            outfile
        )

        print(
            f"{status} "
            f"({size / 1024:.1f} KB)"
        )

        manifest.append({

            "paper":
                src["paper"],

            "doi":
                src["doi"],

            "source_number":
                src["source_number"],

            "description":
                src["description"],

            "filename":
                src["filename"],

            "url":
                src["url"],

            "size_bytes":
                size,

            "sha256":
                sha256_file(outfile),

            "status":
                status
        })

    except Exception as e:

        print("FAILED")

        download_failures.append({

            "paper":
                src["paper"],

            "filename":
                src["filename"],

            "url":
                src["url"],

            "error":
                str(e)
        })


print()

if download_failures:

    fail_out = os.path.join(
        META_DIR,
        "STEP22A_04_SUPPLEMENT_DOWNLOAD_FAILURES.csv"
    )

    write_csv(
        fail_out,
        download_failures,
        [
            "paper",
            "filename",
            "url",
            "error"
        ]
    )

    print(
        "ERROR:",
        len(download_failures),
        "official supplementary files failed."
    )

    print(
        "Saved:"
    )

    print(
        fail_out
    )

    print()
    print(
        "STEP22A_OFFICIAL_SUPPLEMENT_DOWNLOAD_INCOMPLETE"
    )

    sys.exit(1)


manifest_out = os.path.join(
    META_DIR,
    "STEP22A_04_SOURCE_MANIFEST.csv"
)

write_csv(
    manifest_out,
    manifest,
    list(
        manifest[0].keys()
    )
)


# ============================================================
# 4. XLSX parser using Python standard library only
# ============================================================

def extract_xlsx_records(path):

    records = []

    with zipfile.ZipFile(
        path,
        "r"
    ) as z:

        names = set(
            z.namelist()
        )

        shared = []

        if (
            "xl/sharedStrings.xml"
            in names
        ):

            root = ET.fromstring(
                z.read(
                    "xl/sharedStrings.xml"
                )
            )

            for si in root.iter():

                if si.tag.endswith(
                    "}si"
                ):

                    text_parts = []

                    for node in si.iter():

                        if node.tag.endswith(
                            "}t"
                        ):

                            text_parts.append(
                                node.text or ""
                            )

                    shared.append(
                        "".join(
                            text_parts
                        )
                    )

        # Workbook sheet names
        sheet_names = []

        if (
            "xl/workbook.xml"
            in names
        ):

            wb = ET.fromstring(
                z.read(
                    "xl/workbook.xml"
                )
            )

            for node in wb.iter():

                if node.tag.endswith(
                    "}sheet"
                ):

                    sheet_names.append(
                        node.attrib.get(
                            "name",
                            ""
                        )
                    )

        worksheets = sorted(
            [
                x
                for x in names
                if re.fullmatch(
                    r"xl/worksheets/sheet\d+\.xml",
                    x
                )
            ]
        )

        for sheet_idx, sheet_path in enumerate(
            worksheets
        ):

            sheet_name = (
                sheet_names[sheet_idx]
                if sheet_idx < len(sheet_names)
                else os.path.basename(
                    sheet_path
                )
            )

            root = ET.fromstring(
                z.read(
                    sheet_path
                )
            )

            row_number = 0

            for row in root.iter():

                if not row.tag.endswith(
                    "}row"
                ):
                    continue

                row_number += 1

                values = []

                for cell in row:

                    if not cell.tag.endswith(
                        "}c"
                    ):
                        continue

                    ctype = cell.attrib.get(
                        "t",
                        ""
                    )

                    value = ""

                    # inline string
                    if ctype == "inlineStr":

                        parts = []

                        for node in cell.iter():

                            if node.tag.endswith(
                                "}t"
                            ):
                                parts.append(
                                    node.text or ""
                                )

                        value = "".join(
                            parts
                        )

                    else:

                        v_node = None

                        for node in cell:

                            if node.tag.endswith(
                                "}v"
                            ):
                                v_node = node
                                break

                        if v_node is not None:

                            raw = (
                                v_node.text
                                or ""
                            )

                            if ctype == "s":

                                try:
                                    value = shared[
                                        int(raw)
                                    ]

                                except Exception:
                                    value = raw

                            else:
                                value = raw

                    values.append(
                        normalize_space(
                            value
                        )
                    )

                row_text = normalize_space(
                    " | ".join(
                        x
                        for x in values
                        if x
                    )
                )

                if row_text:

                    records.append({

                        "section":
                            sheet_name,

                        "record_number":
                            row_number,

                        "text":
                            row_text
                    })

    return records


# ============================================================
# 5. DOCX parser
# ============================================================

def extract_docx_records(path):

    records = []

    with zipfile.ZipFile(
        path,
        "r"
    ) as z:

        xml_files = [
            x
            for x in z.namelist()
            if (
                x == "word/document.xml"
                or
                re.fullmatch(
                    r"word/header\d+\.xml",
                    x
                )
                or
                re.fullmatch(
                    r"word/footer\d+\.xml",
                    x
                )
            )
        ]

        n = 0

        for xmlfile in xml_files:

            root = ET.fromstring(
                z.read(
                    xmlfile
                )
            )

            # Table rows first.
            for tr in root.iter():

                if not tr.tag.endswith(
                    "}tr"
                ):
                    continue

                cells = []

                for tc in tr:

                    if not tc.tag.endswith(
                        "}tc"
                    ):
                        continue

                    text_parts = []

                    for node in tc.iter():

                        if node.tag.endswith(
                            "}t"
                        ):
                            text_parts.append(
                                node.text or ""
                            )

                    cells.append(
                        normalize_space(
                            " ".join(
                                text_parts
                            )
                        )
                    )

                text = normalize_space(
                    " | ".join(
                        x
                        for x in cells
                        if x
                    )
                )

                if text:

                    n += 1

                    records.append({

                        "section":
                            xmlfile + ":TABLE",

                        "record_number":
                            n,

                        "text":
                            text
                    })

            # Paragraphs.
            for p in root.iter():

                if not p.tag.endswith(
                    "}p"
                ):
                    continue

                parts = []

                for node in p.iter():

                    if node.tag.endswith(
                        "}t"
                    ):
                        parts.append(
                            node.text or ""
                        )

                text = normalize_space(
                    " ".join(
                        parts
                    )
                )

                if text:

                    n += 1

                    records.append({

                        "section":
                            xmlfile + ":PARAGRAPH",

                        "record_number":
                            n,

                        "text":
                            text
                    })

    return records


# ============================================================
# 6. CSV parser
# ============================================================

def extract_csv_records(path):

    records = []

    with open(
        path,
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as f:

        for i, line in enumerate(
            f,
            1
        ):

            text = normalize_space(
                line
            )

            if text:

                records.append({

                    "section":
                        "CSV",

                    "record_number":
                        i,

                    "text":
                        text
                })

    return records


def extract_records(path):

    low = path.lower()

    if low.endswith(
        ".xlsx"
    ):
        return extract_xlsx_records(
            path
        )

    if low.endswith(
        ".docx"
    ):
        return extract_docx_records(
            path
        )

    if low.endswith(
        ".csv"
    ):
        return extract_csv_records(
            path
        )

    return []


# ============================================================
# 7. Build exact barcode probes from our paired samples
# ============================================================

print()
print("===== 3. BUILD LOCAL BARCODE PROBES =====")
print()

matrix_rows = read_csv(
    MATRIX_AUDIT
)

if len(matrix_rows) != 12:

    raise SystemExit(
        "ERROR: STEP22A-03B matrix audit "
        "does not contain 12 paired samples."
    )


barcode_probes = []

for row in matrix_rows:

    barcode_path = os.path.join(
        PROJECT,
        row["barcodes_file"]
    )

    if not os.path.exists(
        barcode_path
    ):
        raise SystemExit(
            "ERROR: missing barcode file:\n"
            + barcode_path
        )

    sample_barcodes = []

    with open(
        barcode_path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for line in f:

            bc = line.strip()

            if not bc:
                continue

            sample_barcodes.append(
                bc
            )

            if (
                len(sample_barcodes)
                >= 30
            ):
                break

    for bc in sample_barcodes:

        barcode_probes.append({

            "patient_id":
                row["patient_id"],

            "timepoint":
                row["timepoint"],

            "sample_name":
                row["sample_name"],

            "hrs_accession":
                row["hrs_accession"],

            "barcode":
                bc
        })


print(
    "Exact barcode probes:",
    len(barcode_probes)
)

print(
    "Samples represented   :",
    len(
        set(
            x["sample_name"]
            for x in barcode_probes
        )
    )
)

print()


# ============================================================
# 8. Parse every supplementary file
# ============================================================

print("===== 4. PARSE OFFICIAL SUPPLEMENTARY CONTENT =====")
print()

all_records = []

records_by_file = {}

for m in manifest:

    path = os.path.join(
        SOURCE_DIR,
        m["filename"]
    )

    try:

        records = extract_records(
            path
        )

    except Exception as e:

        print(
            "ERROR parsing:",
            m["filename"]
        )

        print(
            str(e)
        )

        sys.exit(1)

    records_by_file[
        m["filename"]
    ] = records

    print(
        f"{m['filename']:<42s}"
        f"{len(records):>8d} records"
    )

    for r in records:

        all_records.append({

            "paper":
                m["paper"],

            "doi":
                m["doi"],

            "filename":
                m["filename"],

            "source_number":
                m["source_number"],

            "section":
                r["section"],

            "record_number":
                r["record_number"],

            "text":
                r["text"]
        })


print()
print(
    "Total parsed records:",
    len(all_records)
)
print()


# ============================================================
# 9. Keyword audit
# ============================================================

print("===== 5. KEYWORD AUDIT =====")
print()

patterns = {

    "EASC":
        r"\bEASC\b",

    "ADENOSQUAMOUS":
        r"adenosquamous",

    "ESCC":
        r"\bESCC\b",

    "SQUAMOUS":
        r"squamous",

    "FIBROBLAST":
        r"fibroblast",

    "CAF":
        r"\bCAF\b",

    "CELLTYPE":
        r"cell[\s_-]*type",

    "ANNOTATION":
        r"annotat",

    "BARCODE":
        r"barcode",

    "CLUSTER":
        r"cluster",

    "PATIENT_ID":
        r"\bP\d+\b",

    "SAMPLE_ID":
        r"\bP\d+_[TN]_[AB]\b"
}


keyword_hits = []

for rec in all_records:

    for label, pattern in patterns.items():

        if re.search(
            pattern,
            rec["text"],
            flags=re.I
        ):

            keyword_hits.append({

                "keyword":
                    label,

                "paper":
                    rec["paper"],

                "filename":
                    rec["filename"],

                "section":
                    rec["section"],

                "record_number":
                    rec["record_number"],

                "text":
                    rec["text"][:2000]
            })


keyword_out = os.path.join(
    META_DIR,
    "STEP22A_04_KEYWORD_HITS.csv"
)

write_csv(
    keyword_out,
    keyword_hits,
    [
        "keyword",
        "paper",
        "filename",
        "section",
        "record_number",
        "text"
    ]
)


counts = defaultdict(int)

for x in keyword_hits:
    counts[x["keyword"]] += 1


for key in [
    "EASC",
    "ADENOSQUAMOUS",
    "ESCC",
    "FIBROBLAST",
    "CAF",
    "CELLTYPE",
    "ANNOTATION",
    "BARCODE",
    "CLUSTER"
]:

    print(
        f"{key:<18s}: "
        f"{counts.get(key, 0)}"
    )

print()


# ============================================================
# 10. Histology evidence audit
# ============================================================

print("===== 6. EASC / HISTOLOGY PATIENT AUDIT =====")
print()

histology_evidence = []

patient_pattern = re.compile(
    r"\bP(?:atient\s*)?(\d+)\b",
    flags=re.I
)

for idx, rec in enumerate(
    all_records
):

    if not re.search(
        r"\bEASC\b|adenosquamous",
        rec["text"],
        flags=re.I
    ):
        continue

    # Same record.
    candidate_ids = set()

    for m in patient_pattern.finditer(
        rec["text"]
    ):

        candidate_ids.add(
            "P"
            + str(
                int(
                    m.group(1)
                )
            )
        )

    # Nearby records from the SAME file,
    # useful when table extraction separates text.
    nearby_texts = []

    for j in range(
        max(
            0,
            idx - 2
        ),
        min(
            len(all_records),
            idx + 3
        )
    ):

        near = all_records[j]

        if (
            near["filename"]
            ==
            rec["filename"]
        ):
            nearby_texts.append(
                near["text"]
            )

    nearby_text = normalize_space(
        " || ".join(
            nearby_texts
        )
    )

    nearby_ids = set()

    for m in patient_pattern.finditer(
        nearby_text
    ):

        nearby_ids.add(
            "P"
            + str(
                int(
                    m.group(1)
                )
            )
        )

    histology_evidence.append({

        "paper":
            rec["paper"],

        "filename":
            rec["filename"],

        "section":
            rec["section"],

        "record_number":
            rec["record_number"],

        "histology_keyword":
            (
                "EASC"
                if re.search(
                    r"\bEASC\b",
                    rec["text"],
                    flags=re.I
                )
                else "adenosquamous"
            ),

        "same_record_patient_ids":
            ";".join(
                sorted(
                    candidate_ids,
                    key=lambda x: int(
                        x[1:]
                    )
                )
            ),

        "nearby_patient_ids":
            ";".join(
                sorted(
                    nearby_ids,
                    key=lambda x: int(
                        x[1:]
                    )
                )
            ),

        "source_text":
            rec["text"][:3000],

        "nearby_context":
            nearby_text[:5000]
    })


histology_out = os.path.join(
    META_DIR,
    "STEP22A_04_HISTOLOGY_EVIDENCE.csv"
)

write_csv(
    histology_out,
    histology_evidence,
    [
        "paper",
        "filename",
        "section",
        "record_number",
        "histology_keyword",
        "same_record_patient_ids",
        "nearby_patient_ids",
        "source_text",
        "nearby_context"
    ]
)


same_record_candidates = set()

for x in histology_evidence:

    for p in (
        x["same_record_patient_ids"]
        .split(";")
    ):

        p = p.strip()

        if p:
            same_record_candidates.add(
                p
            )


print(
    "Records containing EASC/adenosquamous:",
    len(histology_evidence)
)

print(
    "Patient IDs on SAME evidence record:",
    (
        ", ".join(
            sorted(
                same_record_candidates,
                key=lambda x: int(
                    x[1:]
                )
            )
        )
        if same_record_candidates
        else "<none>"
    )
)

print()


# ============================================================
# 11. Exact barcode compatibility audit
# ============================================================

print("===== 7. EXACT BARCODE COMPATIBILITY AUDIT =====")
print()

barcode_hits = []

for m in manifest:

    filename = m["filename"]

    records = records_by_file[
        filename
    ]

    # Joining is safe here because these are small supplements.
    full_text = "\n".join(
        r["text"]
        for r in records
    )

    hit_barcodes = set()
    hit_samples = set()

    for probe in barcode_probes:

        bc = probe["barcode"]

        if (
            bc
            and
            bc in full_text
        ):

            hit_barcodes.add(
                bc
            )

            hit_samples.add(
                probe["sample_name"]
            )

            barcode_hits.append({

                "paper":
                    m["paper"],

                "filename":
                    filename,

                "sample_name":
                    probe["sample_name"],

                "hrs_accession":
                    probe["hrs_accession"],

                "barcode":
                    bc
            })

    print(
        f"{filename:<42s}"
        f" exact_barcodes={len(hit_barcodes):>4d}"
        f" samples={len(hit_samples):>2d}"
    )


barcode_hits_out = os.path.join(
    META_DIR,
    "STEP22A_04_EXACT_BARCODE_HITS.csv"
)

write_csv(
    barcode_hits_out,
    barcode_hits,
    [
        "paper",
        "filename",
        "sample_name",
        "hrs_accession",
        "barcode"
    ]
)

print()


# ============================================================
# 12. Annotation-source classification per file
# ============================================================

print("===== 8. AUTHOR ANNOTATION SOURCE CLASSIFICATION =====")
print()

annotation_audit = []

for m in manifest:

    filename = m["filename"]

    records = records_by_file[
        filename
    ]

    text = "\n".join(
        x["text"]
        for x in records
    )

    low = text.lower()

    fibroblast_hits = len(
        re.findall(
            r"fibroblast|\bcaf\b",
            text,
            flags=re.I
        )
    )

    annotation_word_hits = len(
        re.findall(
            r"cell[\s_-]*type|"
            r"annotat|"
            r"cluster|"
            r"barcode",
            text,
            flags=re.I
        )
    )

    sample_id_hits = len(
        re.findall(
            r"\bP\d+_[TN]_[AB]\b",
            text,
            flags=re.I
        )
    )

    exact_hits_for_file = [
        x
        for x in barcode_hits
        if x["filename"]
        == filename
    ]

    exact_barcode_count = len(
        set(
            x["barcode"]
            for x in exact_hits_for_file
        )
    )

    exact_sample_count = len(
        set(
            x["sample_name"]
            for x in exact_hits_for_file
        )
    )

    if (
        exact_barcode_count >= 5
        and
        annotation_word_hits > 0
    ):

        classification = (
            "POSSIBLE_AUTHOR_CELL_LEVEL_ANNOTATION"
        )

    elif (
        exact_barcode_count > 0
    ):

        classification = (
            "BARCODE_COMPATIBLE_DATA_REVIEW"
        )

    elif (
        fibroblast_hits > 0
        and
        annotation_word_hits > 0
    ):

        classification = (
            "AUTHOR_SUMMARY_OR_CLUSTER_LEVEL_DATA"
        )

    elif (
        fibroblast_hits > 0
    ):

        classification = (
            "FIBROBLAST_RELATED_SUMMARY_DATA"
        )

    else:

        classification = (
            "NO_CELL_ANNOTATION_EVIDENCE"
        )

    annotation_audit.append({

        "paper":
            m["paper"],

        "filename":
            filename,

        "description":
            m["description"],

        "parsed_records":
            len(records),

        "fibroblast_or_CAF_keyword_hits":
            fibroblast_hits,

        "annotation_cluster_barcode_keyword_hits":
            annotation_word_hits,

        "sample_id_hits":
            sample_id_hits,

        "exact_local_barcode_hits":
            exact_barcode_count,

        "paired_samples_with_exact_barcode_hits":
            exact_sample_count,

        "classification":
            classification
    })


annotation_out = os.path.join(
    META_DIR,
    "STEP22A_04_AUTHOR_ANNOTATION_SOURCE_AUDIT.csv"
)

write_csv(
    annotation_out,
    annotation_audit,
    list(
        annotation_audit[0].keys()
    )
)


for x in annotation_audit:

    print(
        f"{x['filename']:<42s} "
        f"{x['classification']}"
    )

print()


# ============================================================
# 13. Determine overall annotation status
# ============================================================

cell_level_candidates = [
    x
    for x in annotation_audit
    if x["classification"]
    ==
    "POSSIBLE_AUTHOR_CELL_LEVEL_ANNOTATION"
]

barcode_review_candidates = [
    x
    for x in annotation_audit
    if x["classification"]
    ==
    "BARCODE_COMPATIBLE_DATA_REVIEW"
]

summary_candidates = [
    x
    for x in annotation_audit
    if x["classification"]
    in {
        "AUTHOR_SUMMARY_OR_CLUSTER_LEVEL_DATA",
        "FIBROBLAST_RELATED_SUMMARY_DATA"
    }
]


if cell_level_candidates:

    author_annotation_status = (
        "POSSIBLE_PUBLIC_CELL_LEVEL_ANNOTATION_FOUND"
    )

elif barcode_review_candidates:

    author_annotation_status = (
        "BARCODE_COMPATIBLE_PUBLIC_DATA_FOUND_REVIEW_REQUIRED"
    )

elif summary_candidates:

    author_annotation_status = (
        "PUBLIC_SUMMARY_DATA_ONLY_NO_BARCODE_ANNOTATION"
    )

else:

    author_annotation_status = (
        "NO_PUBLIC_AUTHOR_CELL_ANNOTATION_FOUND"
    )


# ============================================================
# 14. Histology decision
# ============================================================

# We require a UNIQUE patient ID on the SAME record
# containing EASC / adenosquamous before automatic exclusion.

if len(
    same_record_candidates
) == 1:

    easc_mapping_status = (
        "UNIQUE_PUBLIC_EASC_PATIENT_IDENTIFIED"
    )

    easc_patient = next(
        iter(
            same_record_candidates
        )
    )

elif len(
    same_record_candidates
) == 0:

    easc_mapping_status = (
        "EASC_PATIENT_NOT_RESOLVED_FROM_PUBLIC_SUPPLEMENTS"
    )

    easc_patient = ""

else:

    easc_mapping_status = (
        "EASC_PATIENT_MAPPING_AMBIGUOUS"
    )

    easc_patient = ""


# ============================================================
# 15. Frozen original project check
# ============================================================

print("===== 9. FROZEN PROJECT CHECK =====")
print()

frozen_changes = []

if os.path.exists(
    FROZEN_INV
):

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
    "Frozen files changed:",
    len(frozen_changes)
)

if frozen_changes:

    for item in frozen_changes[:20]:
        print(item)

    print()
    print(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )

    sys.exit(1)

print(
    "✓ Frozen Step19-21 files unchanged."
)
print()


# ============================================================
# 16. Summary
# ============================================================

summary_out = os.path.join(
    META_DIR,
    "STEP22A_04_AUTHOR_SOURCE_HISTOLOGY_SUMMARY.txt"
)

with open(
    summary_out,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-4A AUTHOR SOURCE + HISTOLOGY AUDIT\n"
    )

    f.write(
        "=" * 65
        + "\n\n"
    )

    f.write(
        "Official supplementary files downloaded: "
        f"{len(manifest)}\n"
    )

    f.write(
        "Genome Medicine supplements: 9\n"
    )

    f.write(
        "Biomarker Research supplements: 8\n"
    )

    f.write(
        "Exact local barcode probes tested: "
        f"{len(barcode_probes)}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "PRIMARY PAPER HISTOLOGY FACT:\n"
    )

    f.write(
        "One EASC patient is reported in the "
        "22-patient source cohort.\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "EASC mapping status: "
        f"{easc_mapping_status}\n"
    )

    f.write(
        "EASC patient: "
        f"{easc_patient if easc_patient else 'UNRESOLVED'}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Author annotation status: "
        f"{author_annotation_status}\n"
    )

    f.write(
        "Possible public cell-level files: "
        f"{len(cell_level_candidates)}\n"
    )

    f.write(
        "Barcode-compatible review files: "
        f"{len(barcode_review_candidates)}\n"
    )

    f.write(
        "Summary/cluster-level candidate files: "
        f"{len(summary_candidates)}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Cell annotation performed: NO\n"
    )

    f.write(
        "Samples excluded in Step22A-4A: 0\n"
    )

    f.write(
        "Frozen Step19-21 changed: NO\n"
    )

    f.write(
        "Manuscript changed: NO\n"
    )


# ============================================================
# FINAL PRINT
# ============================================================

print("=" * 78)
print("STEP22A-4A COMPLETE")
print("=" * 78)
print()

print(
    "OFFICIAL SUPPLEMENTS        :",
    len(manifest)
)

print(
    "  GENOME MEDICINE           : 9"
)

print(
    "  BIOMARKER RESEARCH        : 8"
)

print()

print(
    "PRIMARY PAPER REPORTS EASC  : YES"
)

print(
    "EASC MAPPING STATUS         :",
    easc_mapping_status
)

print(
    "EASC PATIENT                :",
    (
        easc_patient
        if easc_patient
        else "UNRESOLVED"
    )
)

print()

print(
    "AUTHOR ANNOTATION STATUS    :",
    author_annotation_status
)

print(
    "CELL-LEVEL CANDIDATES       :",
    len(cell_level_candidates)
)

print(
    "BARCODE REVIEW CANDIDATES   :",
    len(barcode_review_candidates)
)

print(
    "SUMMARY/CLUSTER CANDIDATES  :",
    len(summary_candidates)
)

print()

print(
    "EXACT BARCODE HITS TOTAL    :",
    len(
        set(
            (
                x["filename"],
                x["barcode"]
            )
            for x in barcode_hits
        )
    )
)

print()

print(
    "CELL ANNOTATION PERFORMED   : NO"
)

print(
    "SAMPLES EXCLUDED NOW        : 0"
)

print(
    "NEW OMIX TAR DOWNLOADED     : 0"
)

print()

print(
    "FROZEN STEP19-21 MODIFIED   : NO"
)

print(
    "MANUSCRIPT MODIFIED         : NO"
)

print()

print("OUTPUTS:")

for path in [
    manifest_out,
    histology_out,
    annotation_out,
    barcode_hits_out,
    keyword_out,
    summary_out
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
    easc_mapping_status
    ==
    "UNIQUE_PUBLIC_EASC_PATIENT_IDENTIFIED"
):

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "STEP22A-4B = CORRECT ESCC ELIGIBILITY "
        "+ VALIDATE AUTHOR ANNOTATION ROUTE"
    )

elif (
    easc_mapping_status
    ==
    "EASC_PATIENT_NOT_RESOLVED_FROM_PUBLIC_SUPPLEMENTS"
):

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "HISTOLOGY_MAPPING_INCOMPLETE"
    )

    print(
        "Do NOT start primary CORE2 validation yet."
    )

else:

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "HISTOLOGY_MAPPING_AMBIGUOUS"
    )

    print(
        "Do NOT start primary CORE2 validation yet."
    )

print()

print("STOP HERE.")
print(
    "DO NOT START DE NOVO CAF ANNOTATION YET."
)
print()
