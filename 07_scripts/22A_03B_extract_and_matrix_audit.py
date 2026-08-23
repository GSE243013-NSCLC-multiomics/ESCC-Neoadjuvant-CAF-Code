#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-3B
#
# Paired tumor TAR extraction + matrix integrity audit
#
# IMPORTANT:
# - Extract ONLY the 12 already-downloaded paired tumor TARs
# - DO NOT download additional files
# - DO NOT run Seurat / clustering / QC filtering
# - DO NOT annotate fibroblasts yet
# - DO NOT modify Step19-21 or manuscript
# ============================================================

import csv
import os
import re
import shutil
import sys
import tarfile
from collections import Counter
from datetime import datetime


PROJECT = os.getcwd()

META_DIR = os.path.join(
    PROJECT,
    "00_metadata",
    "OMIX005710"
)

PAIRED_DIR = os.path.join(
    PROJECT,
    "01_raw",
    "OMIX005710",
    "paired_tumor"
)

EXTRACT_ROOT = os.path.join(
    PROJECT,
    "01_raw",
    "OMIX005710",
    "extracted"
)

TARGETS_FILE = os.path.join(
    META_DIR,
    "STEP22A_03A_PAIRED_DOWNLOAD_TARGETS.tsv"
)

FROZEN_INV = os.path.join(
    META_DIR,
    "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

os.makedirs(EXTRACT_ROOT, exist_ok=True)


print()
print("=" * 76)
print("STEP22A-3B: PAIRED TAR EXTRACTION + MATRIX AUDIT")
print("=" * 76)
print()


# ============================================================
# Helpers
# ============================================================

def read_csv(path, delimiter=","):

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
            csv.DictReader(
                f,
                delimiter=delimiter
            )
        )


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


def count_lines(path):

    n = 0

    with open(
        path,
        "rb"
    ) as f:

        for _ in f:
            n += 1

    return n


def read_matrix_market_dims(path):

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        first = f.readline().strip()

        if not first.lower().startswith(
            "%%matrixmarket"
        ):
            raise RuntimeError(
                "Not a Matrix Market file: "
                + path
            )

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith("%"):
                continue

            parts = line.split()

            if len(parts) != 3:
                raise RuntimeError(
                    "Cannot parse Matrix Market "
                    "dimension line: "
                    + line
                )

            nrow = int(parts[0])
            ncol = int(parts[1])
            nnz = int(parts[2])

            return nrow, ncol, nnz

    raise RuntimeError(
        "No dimension line found: "
        + path
    )


def safe_extract(tar_path, destination):
    """
    Safe extraction:
    - no absolute paths
    - no ../ traversal
    - no symlinks / hardlinks
    - regular files + directories only
    """

    root = os.path.realpath(destination)

    with tarfile.open(
        tar_path,
        "r:*"
    ) as tar:

        members = tar.getmembers()

        for member in members:

            name = member.name

            if os.path.isabs(name):
                raise RuntimeError(
                    "Unsafe absolute TAR path: "
                    + name
                )

            target = os.path.realpath(
                os.path.join(
                    destination,
                    name
                )
            )

            if not (
                target == root
                or
                target.startswith(
                    root + os.sep
                )
            ):
                raise RuntimeError(
                    "Unsafe TAR path traversal: "
                    + name
                )

            if (
                member.issym()
                or
                member.islnk()
            ):
                raise RuntimeError(
                    "Symlink/hardlink rejected: "
                    + name
                )

            if not (
                member.isfile()
                or
                member.isdir()
            ):
                raise RuntimeError(
                    "Unsupported TAR member: "
                    + name
                )

        for member in members:

            if member.isdir():

                out = os.path.join(
                    destination,
                    member.name
                )

                os.makedirs(
                    out,
                    exist_ok=True
                )

            elif member.isfile():

                out = os.path.join(
                    destination,
                    member.name
                )

                os.makedirs(
                    os.path.dirname(out),
                    exist_ok=True
                )

                src = tar.extractfile(member)

                if src is None:
                    raise RuntimeError(
                        "Cannot extract: "
                        + member.name
                    )

                with src:
                    with open(
                        out,
                        "wb"
                    ) as dst:
                        shutil.copyfileobj(
                            src,
                            dst
                        )


def recursive_files(root):

    out = []

    for dirpath, dirs, files in os.walk(root):

        for filename in files:

            if filename.startswith(
                ".STEP22A_"
            ):
                continue

            out.append(
                os.path.join(
                    dirpath,
                    filename
                )
            )

    return sorted(out)


def find_component(files, component):

    matches = []

    for path in files:

        b = os.path.basename(
            path
        ).lower()

        if component == "matrix":

            if (
                b == "matrix.mtx"
                or
                b.endswith(
                    "_matrix.mtx"
                )
            ):
                matches.append(path)

        elif component == "barcodes":

            if (
                b == "barcodes.tsv"
                or
                b.endswith(
                    "_barcodes.tsv"
                )
            ):
                matches.append(path)

        elif component == "genes":

            if (
                b == "genes.tsv"
                or
                b == "features.tsv"
                or
                b.endswith(
                    "_genes.tsv"
                )
                or
                b.endswith(
                    "_features.tsv"
                )
            ):
                matches.append(path)

    return matches


def read_gene_table(path):

    rows = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for line in f:

            line = line.rstrip("\r\n")

            if not line:
                continue

            parts = line.split("\t")

            gene_id = (
                parts[0].strip()
                if len(parts) >= 1
                else ""
            )

            gene_symbol = (
                parts[1].strip()
                if len(parts) >= 2
                else gene_id
            )

            rows.append(
                (
                    gene_id,
                    gene_symbol
                )
            )

    return rows


# ============================================================
# 1. Load paired target list
# ============================================================

print("===== 1. LOAD PAIRED TARGET LIST =====")
print()

targets = read_csv(
    TARGETS_FILE,
    delimiter="\t"
)

print(
    "Target records:",
    len(targets)
)

if len(targets) != 12:

    raise SystemExit(
        "ERROR: expected exactly 12 paired targets."
    )


expected_patients = {
    "P6",
    "P9",
    "P15",
    "P16",
    "P18",
    "P19"
}

observed_patients = {
    x["patient_id"]
    for x in targets
}

if observed_patients != expected_patients:

    raise SystemExit(
        "ERROR: paired patient set mismatch."
    )

print(
    "Paired patients:",
    ", ".join(
        sorted(
            observed_patients,
            key=lambda x: int(x[1:])
        )
    )
)

print()


# ============================================================
# 2. Confirm 12 TARs
# ============================================================

print("===== 2. CONFIRM LOCAL TAR FILES =====")
print()

for x in targets:

    tar_path = os.path.join(
        PAIRED_DIR,
        x["filename"]
    )

    if not os.path.exists(
        tar_path
    ):
        raise SystemExit(
            "ERROR: missing TAR:\n"
            + tar_path
        )

print("✓ All 12 TAR files found.")
print()


# ============================================================
# 3. Estimate extraction size before extracting
# ============================================================

print("===== 3. EXTRACTION SIZE PREFLIGHT =====")
print()

total_archive_bytes = 0
total_uncompressed_bytes = 0
archive_rows = []

for x in targets:

    tar_path = os.path.join(
        PAIRED_DIR,
        x["filename"]
    )

    archive_size = os.path.getsize(
        tar_path
    )

    total_archive_bytes += archive_size

    with tarfile.open(
        tar_path,
        "r:*"
    ) as tar:

        members = [
            m
            for m in tar.getmembers()
            if m.isfile()
        ]

        expanded = sum(
            m.size
            for m in members
        )

    total_uncompressed_bytes += expanded

    archive_rows.append({

        "patient_id":
            x["patient_id"],

        "timepoint":
            x["timepoint"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            x["hrs_accession"],

        "filename":
            x["filename"],

        "archive_size_MB":
            round(
                archive_size / 1024**2,
                3
            ),

        "estimated_extracted_MB":
            round(
                expanded / 1024**2,
                3
            ),

        "internal_file_count":
            len(members)
    })


disk = shutil.disk_usage(
    PROJECT
)

free_gib = disk.free / 1024**3
expanded_gib = (
    total_uncompressed_bytes
    / 1024**3
)

print(
    "Compressed TAR total :",
    round(
        total_archive_bytes / 1024**3,
        3
    ),
    "GiB"
)

print(
    "Estimated extraction :",
    round(
        expanded_gib,
        3
    ),
    "GiB"
)

print(
    "Current free space   :",
    round(
        free_gib,
        2
    ),
    "GiB"
)

projected_free = (
    free_gib - expanded_gib
)

print(
    "Projected free space :",
    round(
        projected_free,
        2
    ),
    "GiB"
)

print()

if projected_free < 20:

    print(
        "ERROR: extraction would leave "
        "less than 20 GiB free."
    )

    print(
        "STEP22A_EXTRACTION_DISK_GUARD"
    )

    sys.exit(1)

print(
    "✓ Disk-space guard passed."
)
print()


archive_preflight_out = os.path.join(
    META_DIR,
    "STEP22A_03B_EXTRACTION_PREFLIGHT.csv"
)

write_csv(
    archive_preflight_out,
    archive_rows,
    list(
        archive_rows[0].keys()
    )
)


# ============================================================
# 4. Extract each TAR safely
# ============================================================

print("===== 4. SAFE EXTRACTION =====")
print()

extraction_manifest = []

for i, x in enumerate(
    targets,
    1
):

    hrs = x["hrs_accession"]

    tar_path = os.path.join(
        PAIRED_DIR,
        x["filename"]
    )

    sample_dir = os.path.join(
        EXTRACT_ROOT,
        hrs
    )

    marker = os.path.join(
        sample_dir,
        ".STEP22A_EXTRACT_COMPLETE"
    )

    print(
        f"[{i:02d}/12] "
        f"{x['patient_id']} "
        f"{x['timepoint']} "
        f"{x['sample_name']} "
        f"({hrs})"
    )

    extraction_status = ""

    if os.path.exists(marker):

        print(
            "  ✓ Already extracted; "
            "marker present."
        )

        extraction_status = (
            "ALREADY_EXTRACTED_VERIFIED"
        )

    else:

        if os.path.exists(
            sample_dir
        ):

            print(
                "  Partial extraction directory "
                "detected; rebuilding safely."
            )

            shutil.rmtree(
                sample_dir
            )

        os.makedirs(
            sample_dir,
            exist_ok=True
        )

        safe_extract(
            tar_path,
            sample_dir
        )

        with open(
            marker,
            "w",
            encoding="utf-8"
        ) as f:

            f.write(
                datetime.now().isoformat(
                    timespec="seconds"
                )
                + "\n"
            )

        print("  ✓ Extraction complete.")

        extraction_status = (
            "EXTRACTED_NOW"
        )

    files_now = recursive_files(
        sample_dir
    )

    extraction_manifest.append({

        "patient_id":
            x["patient_id"],

        "timepoint":
            x["timepoint"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            hrs,

        "source_tar":
            x["filename"],

        "extract_directory":
            os.path.relpath(
                sample_dir,
                PROJECT
            ),

        "extracted_file_count":
            len(files_now),

        "status":
            extraction_status
    })


print()

extraction_manifest_out = os.path.join(
    META_DIR,
    "STEP22A_03B_EXTRACTION_MANIFEST.csv"
)

write_csv(
    extraction_manifest_out,
    extraction_manifest,
    list(
        extraction_manifest[0].keys()
    )
)


# ============================================================
# 5. Detect matrix / barcode / gene files
# ============================================================

print(
    "===== 5. MATRIX COMPONENT AUDIT ====="
)
print()

matrix_audit = []
fatal = []

for x in targets:

    hrs = x["hrs_accession"]

    sample_dir = os.path.join(
        EXTRACT_ROOT,
        hrs
    )

    files_now = recursive_files(
        sample_dir
    )

    matrices = find_component(
        files_now,
        "matrix"
    )

    barcodes = find_component(
        files_now,
        "barcodes"
    )

    genes = find_component(
        files_now,
        "genes"
    )

    if (
        len(matrices) != 1
        or
        len(barcodes) != 1
        or
        len(genes) != 1
    ):

        fatal.append({

            "sample_name":
                x["sample_name"],

            "hrs":
                hrs,

            "matrix_count":
                len(matrices),

            "barcode_count":
                len(barcodes),

            "gene_count":
                len(genes)
        })

        continue

    matrix_path = matrices[0]
    barcode_path = barcodes[0]
    gene_path = genes[0]

    nrow, ncol, nnz = (
        read_matrix_market_dims(
            matrix_path
        )
    )

    n_barcodes = count_lines(
        barcode_path
    )

    n_genes = count_lines(
        gene_path
    )

    dimension_status = (
        "PASS"
        if (
            nrow == n_genes
            and
            ncol == n_barcodes
        )
        else "FAIL"
    )

    density = (
        nnz / (nrow * ncol)
        if nrow > 0 and ncol > 0
        else 0
    )

    matrix_audit.append({

        "patient_id":
            x["patient_id"],

        "timepoint":
            x["timepoint"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            hrs,

        "matrix_file":
            os.path.relpath(
                matrix_path,
                PROJECT
            ),

        "barcodes_file":
            os.path.relpath(
                barcode_path,
                PROJECT
            ),

        "genes_file":
            os.path.relpath(
                gene_path,
                PROJECT
            ),

        "matrix_rows":
            nrow,

        "matrix_columns":
            ncol,

        "matrix_nonzero":
            nnz,

        "barcode_count":
            n_barcodes,

        "gene_count":
            n_genes,

        "matrix_density":
            round(
                density,
                8
            ),

        "dimension_status":
            dimension_status
    })


if fatal:

    print(
        "ERROR: unexpected component count."
    )

    for x in fatal:
        print(x)

    print()
    print(
        "STEP22A_MATRIX_COMPONENT_AUDIT_FAILED"
    )

    sys.exit(1)


print(
    f"{'PATIENT':<8}"
    f"{'TIME':<8}"
    f"{'SAMPLE':<12}"
    f"{'GENES':>8}"
    f"{'CELLS':>10}"
    f"{'STATUS':>10}"
)

print("-" * 62)

for x in matrix_audit:

    print(
        f"{x['patient_id']:<8}"
        f"{x['timepoint']:<8}"
        f"{x['sample_name']:<12}"
        f"{x['gene_count']:>8}"
        f"{x['barcode_count']:>10}"
        f"{x['dimension_status']:>10}"
    )


failed_dimensions = [
    x
    for x in matrix_audit
    if x["dimension_status"]
    != "PASS"
]

print()

print(
    "Dimension checks passed:",
    len(matrix_audit)
    -
    len(failed_dimensions),
    "/",
    len(matrix_audit)
)

if failed_dimensions:

    print()
    print(
        "ERROR: matrix/barcode/gene "
        "dimensions do not match."
    )

    print(
        "STEP22A_MATRIX_DIMENSION_MISMATCH"
    )

    sys.exit(1)


matrix_audit_out = os.path.join(
    META_DIR,
    "STEP22A_03B_MATRIX_AUDIT.csv"
)

write_csv(
    matrix_audit_out,
    matrix_audit,
    list(
        matrix_audit[0].keys()
    )
)


# ============================================================
# 6. Barcode uniqueness audit
# ============================================================

print()
print(
    "===== 6. BARCODE UNIQUENESS AUDIT ====="
)
print()

barcode_audit = []

for x in matrix_audit:

    path = os.path.join(
        PROJECT,
        x["barcodes_file"]
    )

    barcodes = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for line in f:

            value = line.strip()

            if value:
                barcodes.append(value)

    unique_count = len(
        set(barcodes)
    )

    duplicates = (
        len(barcodes)
        -
        unique_count
    )

    status = (
        "PASS"
        if duplicates == 0
        else "FAIL"
    )

    barcode_audit.append({

        "patient_id":
            x["patient_id"],

        "timepoint":
            x["timepoint"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            x["hrs_accession"],

        "barcode_count":
            len(barcodes),

        "unique_barcodes":
            unique_count,

        "duplicate_barcodes":
            duplicates,

        "status":
            status
    })


bad_barcode = [
    x for x in barcode_audit
    if x["status"] != "PASS"
]

print(
    "Samples with unique barcodes:",
    len(barcode_audit)
    -
    len(bad_barcode),
    "/",
    len(barcode_audit)
)

if bad_barcode:

    print()
    print(
        "ERROR: duplicate barcodes detected."
    )

    print(
        "STEP22A_DUPLICATE_BARCODES"
    )

    sys.exit(1)


barcode_out = os.path.join(
    META_DIR,
    "STEP22A_03B_BARCODE_AUDIT.csv"
)

write_csv(
    barcode_out,
    barcode_audit,
    list(
        barcode_audit[0].keys()
    )
)


# ============================================================
# 7. Gene table + marker availability audit
#
# IMPORTANT:
# Marker presence only.
# NO classification and NO selection here.
# ============================================================

print()
print(
    "===== 7. GENE / MARKER AVAILABILITY AUDIT ====="
)
print()

marker_groups = {

    "fibroblast_positive": [
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

    "epithelial_exclusion": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT5",
        "KRT14"
    ],

    "immune_exclusion": [
        "PTPRC",
        "CD3D",
        "CD3E",
        "NKG7",
        "LST1",
        "MS4A1",
        "CD79A"
    ],

    "endothelial_exclusion": [
        "PECAM1",
        "VWF",
        "EMCN",
        "KDR"
    ],

    "mural_pericyte": [
        "RGS5",
        "CSPG4",
        "MCAM",
        "NOTCH3"
    ]
}


marker_rows = []
gene_audit = []

for x in matrix_audit:

    gene_path = os.path.join(
        PROJECT,
        x["genes_file"]
    )

    genes = read_gene_table(
        gene_path
    )

    ids = [
        g[0]
        for g in genes
    ]

    symbols = [
        g[1]
        for g in genes
    ]

    symbols_upper = {
        g.upper()
        for g in symbols
        if g
    }

    duplicate_symbols = (
        len(symbols)
        -
        len(set(symbols))
    )

    gene_audit.append({

        "patient_id":
            x["patient_id"],

        "timepoint":
            x["timepoint"],

        "sample_name":
            x["sample_name"],

        "hrs_accession":
            x["hrs_accession"],

        "gene_rows":
            len(genes),

        "unique_gene_ids":
            len(set(ids)),

        "unique_gene_symbols":
            len(set(symbols)),

        "duplicate_gene_symbol_rows":
            duplicate_symbols,

        "first_gene_id":
            ids[0] if ids else "",

        "first_gene_symbol":
            symbols[0] if symbols else ""
    })


    for group, markers in (
        marker_groups.items()
    ):

        for marker in markers:

            marker_rows.append({

                "patient_id":
                    x["patient_id"],

                "timepoint":
                    x["timepoint"],

                "sample_name":
                    x["sample_name"],

                "hrs_accession":
                    x["hrs_accession"],

                "marker_group":
                    group,

                "gene":
                    marker,

                "present_in_gene_table":
                    (
                        "YES"
                        if marker.upper()
                        in symbols_upper
                        else "NO"
                    )
            })


gene_audit_out = os.path.join(
    META_DIR,
    "STEP22A_03B_GENE_TABLE_AUDIT.csv"
)

write_csv(
    gene_audit_out,
    gene_audit,
    list(
        gene_audit[0].keys()
    )
)


marker_out = os.path.join(
    META_DIR,
    "STEP22A_03B_MARKER_GENE_COVERAGE.csv"
)

write_csv(
    marker_out,
    marker_rows,
    list(
        marker_rows[0].keys()
    )
)


for group in marker_groups:

    subset = [
        x
        for x in marker_rows
        if x["marker_group"] == group
    ]

    present = sum(
        1
        for x in subset
        if x["present_in_gene_table"]
        == "YES"
    )

    total = len(subset)

    print(
        f"{group:<25s} "
        f"{present}/{total} "
        "sample-marker checks present"
    )


print()
print(
    "NOTE: marker presence was audited only."
)

print(
    "No cell was annotated or selected."
)


# ============================================================
# 8. Cell-count summary
# ============================================================

print()
print(
    "===== 8. RAW CELL COUNT SUMMARY ====="
)
print()

total_cells = sum(
    x["barcode_count"]
    for x in matrix_audit
)

before_cells = sum(
    x["barcode_count"]
    for x in matrix_audit
    if x["timepoint"] == "Before"
)

after_cells = sum(
    x["barcode_count"]
    for x in matrix_audit
    if x["timepoint"] == "After"
)

print(
    "Raw paired tumor barcodes:",
    total_cells
)

print(
    "Before-treatment barcodes:",
    before_cells
)

print(
    "After-treatment barcodes :",
    after_cells
)

print()

print(
    "IMPORTANT: barcode count is NOT "
    "a biological replicate count."
)


# ============================================================
# 9. Confirm no author annotation appeared after extraction
# ============================================================

print()
print(
    "===== 9. AUTHOR ANNOTATION FILE AUDIT ====="
)
print()

annotation_hits = []

for x in targets:

    sample_dir = os.path.join(
        EXTRACT_ROOT,
        x["hrs_accession"]
    )

    for path in recursive_files(
        sample_dir
    ):

        b = os.path.basename(
            path
        ).lower()

        if any(
            key in b
            for key in [
                "annotation",
                "celltype",
                "cell_type",
                "metadata",
                "cluster"
            ]
        ):

            annotation_hits.append({

                "patient_id":
                    x["patient_id"],

                "timepoint":
                    x["timepoint"],

                "sample_name":
                    x["sample_name"],

                "hrs_accession":
                    x["hrs_accession"],

                "file":
                    os.path.relpath(
                        path,
                        PROJECT
                    )
            })


annotation_out = os.path.join(
    META_DIR,
    "STEP22A_03B_AUTHOR_ANNOTATION_FILE_AUDIT.csv"
)

annotation_fields = [
    "patient_id",
    "timepoint",
    "sample_name",
    "hrs_accession",
    "file"
]

write_csv(
    annotation_out,
    annotation_hits,
    annotation_fields
)

print(
    "Annotation-like extracted files:",
    len(annotation_hits)
)

if annotation_hits:

    for x in annotation_hits:
        print(
            " ",
            x["sample_name"],
            "->",
            x["file"]
        )

else:

    print(
        "  <none detected>"
    )

print()


# ============================================================
# 10. Frozen Step19-21 file check
# ============================================================

print(
    "===== 10. FROZEN PROJECT CHECK ====="
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


# ============================================================
# 11. Disk space after extraction
# ============================================================

print()
print(
    "===== 11. DISK SPACE AFTER EXTRACTION ====="
)
print()

disk_after = shutil.disk_usage(
    PROJECT
)

print(
    "Free space:",
    round(
        disk_after.free / 1024**3,
        2
    ),
    "GiB"
)


# ============================================================
# 12. Summary
# ============================================================

summary_out = os.path.join(
    META_DIR,
    "STEP22A_03B_EXTRACTION_MATRIX_AUDIT_SUMMARY.txt"
)

with open(
    summary_out,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-3B EXTRACTION + MATRIX AUDIT\n"
    )

    f.write("=" * 60 + "\n\n")

    f.write(
        "Paired patients: 6\n"
    )

    f.write(
        "Paired tumor samples: 12\n"
    )

    f.write(
        "TARs extracted: 12\n"
    )

    f.write(
        "Matrix dimension PASS: "
        f"{len(matrix_audit)} / 12\n"
    )

    f.write(
        "Barcode uniqueness PASS: "
        f"{len(barcode_audit)} / 12\n"
    )

    f.write(
        "Raw paired tumor barcodes: "
        f"{total_cells}\n"
    )

    f.write(
        "Before raw barcodes: "
        f"{before_cells}\n"
    )

    f.write(
        "After raw barcodes: "
        f"{after_cells}\n"
    )

    f.write(
        "Annotation-like files after extraction: "
        f"{len(annotation_hits)}\n"
    )

    f.write(
        "Cell annotation performed: NO\n"
    )

    f.write(
        "QC filtering performed: NO\n"
    )

    f.write(
        "Clustering performed: NO\n"
    )

    f.write(
        "Seurat object created: NO\n"
    )

    f.write(
        "Frozen Step19-21 changed: NO\n"
    )

    f.write(
        "Manuscript changed: NO\n"
    )


print()
print("=" * 76)
print("STEP22A-3B COMPLETE")
print("=" * 76)
print()

print(
    "PAIRED PATIENTS              : 6"
)

print(
    "PAIRED TUMOR SAMPLES         : 12"
)

print(
    "TAR EXTRACTION               : COMPLETE"
)

print()

print(
    "MATRIX DIMENSION PASS        :",
    len(matrix_audit),
    "/ 12"
)

print(
    "BARCODE UNIQUENESS PASS      :",
    len(barcode_audit),
    "/ 12"
)

print()

print(
    "RAW PAIRED TUMOR BARCODES    :",
    total_cells
)

print(
    "  BEFORE                     :",
    before_cells
)

print(
    "  AFTER                      :",
    after_cells
)

print()

print(
    "AUTHOR ANNOTATION FILES      :",
    len(annotation_hits)
)

print()

print(
    "CELL ANNOTATION PERFORMED    : NO"
)

print(
    "QC FILTERING PERFORMED       : NO"
)

print(
    "CLUSTERING PERFORMED         : NO"
)

print(
    "SEURAT OBJECT CREATED        : NO"
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
    archive_preflight_out,
    extraction_manifest_out,
    matrix_audit_out,
    barcode_out,
    gene_audit_out,
    marker_out,
    annotation_out,
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

if len(annotation_hits) == 0:

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "STEP22A-4 = AUTHOR ANNOTATION SOURCE AUDIT"
    )

    print(
        "Do NOT perform de novo CAF annotation yet."
    )

else:

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "STEP22A-4 = VALIDATE AUTHOR ANNOTATION "
        "BARCODE COMPATIBILITY"
    )

print()

print("STOP HERE.")
print()
