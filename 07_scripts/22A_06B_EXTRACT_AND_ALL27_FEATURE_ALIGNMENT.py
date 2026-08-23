#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-6B
#
# SAFE EXTRACT REMAINING 15 TUMORS
# + STRUCTURAL MATRIX AUDIT
# + TRUE ALL-27 ENSG INTERSECTION
#
# IMPORTANT:
#
# We DO NOT force the remaining 15 samples onto the old
# 18,354-feature set.
#
# Instead:
#   - 18,354 = intersection of the original paired 12
#   - new final technical feature universe =
#       intersection across ALL 27 tumor observations
#
# The new universe must be a SUBSET of the frozen 18,354 set.
#
# NO union + zero filling.
#
# This step does NOT:
#   - read matrix expression values
#   - perform cell QC
#   - perform cell annotation
#   - calculate CPM
#   - calculate z-scores
#   - calculate CORE2 / CORE4
#   - perform Before/After inference
#
# Histology remains unresolved.
# ============================================================

from pathlib import Path, PurePosixPath
from datetime import datetime
import csv
import hashlib
import os
import shutil
import sys
import tarfile


PROJECT = Path.cwd()

META_DIR = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

RAW_ROOT = (
    PROJECT
    / "01_raw"
    / "OMIX005710"
)

REMAINING_TAR_DIR = (
    RAW_ROOT
    / "remaining_tumor"
)

# Use same extracted root already used by paired 12.
EXTRACT_ROOT = (
    RAW_ROOT
    / "extracted"
)

COMBINED_INV = (
    META_DIR
    / "STEP22A_06A_COMBINED_27_TUMOR_TECHNICAL_INVENTORY.csv"
)

PAIRED_MATRIX_AUDIT = (
    META_DIR
    / "STEP22A_03B_MATRIX_AUDIT.csv"
)

OLD_COMMON_FEATURES = (
    META_DIR
    / "STEP22A_05A_COMMON_ENSG_FEATURES.csv"
)

GATE_FILE = (
    META_DIR
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

PREFLIGHT_OUT = (
    META_DIR
    / "STEP22A_06B_EXTRACTION_PREFLIGHT.csv"
)

REMAINING_MATRIX_OUT = (
    META_DIR
    / "STEP22A_06B_REMAINING_15_MATRIX_AUDIT.csv"
)

ALL27_FEATURE_COUNT_OUT = (
    META_DIR
    / "STEP22A_06B_ALL27_FEATURE_TABLE_COUNTS.csv"
)

ALL27_COMMON_OUT = (
    META_DIR
    / "STEP22A_06B_ALL27_COMMON_ENSG_FEATURES.csv"
)

ALIGNMENT_OUT = (
    META_DIR
    / "STEP22A_06B_ALL27_FEATURE_ALIGNMENT_MANIFEST.csv"
)

MARKER_OUT = (
    META_DIR
    / "STEP22A_06B_ALL27_FROZEN_MARKER_COVERAGE.csv"
)

SUMMARY_OUT = (
    META_DIR
    / "STEP22A_06B_EXTRACTION_ALIGNMENT_SUMMARY.txt"
)


EXTRACT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)


print()
print("=" * 80)
print("STEP22A-6B: SAFE EXTRACTION + TRUE ALL-27 ENSG ALIGNMENT")
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


def read_csv(
    path,
    delimiter=","
):

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
            csv.DictReader(
                f,
                delimiter=delimiter
            )
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


def sha256_file(path):

    h = hashlib.sha256()

    with open(
        path,
        "rb"
    ) as f:

        while True:

            chunk = f.read(
                8 * 1024 * 1024
            )

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


def count_lines(path):

    n = 0

    with open(
        path,
        "rb"
    ) as f:

        for _ in f:
            n += 1

    return n


def read_matrix_dims(path):

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
                "Not Matrix Market: "
                + str(path)
            )

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith("%"):
                continue

            fields = line.split()

            if len(fields) != 3:

                raise RuntimeError(
                    "Invalid Matrix Market dimension line:\n"
                    + line
                )

            return (
                int(fields[0]),
                int(fields[1]),
                int(fields[2])
            )

    raise RuntimeError(
        "Matrix dimension line not found:\n"
        + str(path)
    )


def read_gene_table(path):

    ids = []
    symbols = []

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

            fields = line.split(
                "\t"
            )

            gid = clean(
                fields[0]
                if len(fields) >= 1
                else ""
            )

            symbol = clean(
                fields[1]
                if len(fields) >= 2
                else gid
            )

            if not gid:

                raise RuntimeError(
                    "Empty gene ID in:\n"
                    + str(path)
                )

            ids.append(
                gid
            )

            symbols.append(
                symbol
                if symbol
                else gid
            )

    if len(ids) != len(
        set(ids)
    ):

        seen = set()
        duplicates = []

        for gid in ids:

            if gid in seen:
                duplicates.append(
                    gid
                )
            else:
                seen.add(
                    gid
                )

        raise RuntimeError(
            "Duplicate ENSG IDs in "
            + str(path)
            + " : "
            + ", ".join(
                sorted(
                    set(
                        duplicates
                    )
                )[:20]
            )
        )

    return ids, symbols


def find_components(root):

    files = []

    for dirpath, _, filenames in os.walk(
        root
    ):

        for filename in filenames:

            if filename.startswith(
                ".STEP22A_"
            ):
                continue

            files.append(
                Path(dirpath)
                / filename
            )

    matrices = []

    barcodes = []

    genes = []

    for path in files:

        base = (
            path.name
            .lower()
        )

        if base == "matrix.mtx":

            matrices.append(
                path
            )

        elif base == "barcodes.tsv":

            barcodes.append(
                path
            )

        elif base in {
            "genes.tsv",
            "features.tsv"
        }:

            genes.append(
                path
            )

    return (
        sorted(
            matrices
        ),
        sorted(
            barcodes
        ),
        sorted(
            genes
        ),
        sorted(
            files
        )
    )


def safe_extract(
    tar_path,
    destination
):

    destination.mkdir(
        parents=True,
        exist_ok=True
    )

    root = destination.resolve()

    with tarfile.open(
        tar_path,
        "r:*"
    ) as tar:

        members = tar.getmembers()

        # First pass: safety validation.
        for member in members:

            posix = PurePosixPath(
                member.name
            )

            if (
                posix.is_absolute()
                or
                ".." in posix.parts
            ):

                raise RuntimeError(
                    "Unsafe TAR path: "
                    + member.name
                )

            if (
                member.issym()
                or
                member.islnk()
            ):

                raise RuntimeError(
                    "Symlink/hardlink rejected: "
                    + member.name
                )

            if not (
                member.isfile()
                or
                member.isdir()
            ):

                raise RuntimeError(
                    "Unsupported TAR member: "
                    + member.name
                )

            target = (
                destination
                / member.name
            ).resolve()

            if not (
                target == root
                or
                str(target).startswith(
                    str(root)
                    + os.sep
                )
            ):

                raise RuntimeError(
                    "Path traversal rejected: "
                    + member.name
                )

        # Second pass: extraction.
        for member in members:

            target = (
                destination
                / member.name
            )

            if member.isdir():

                target.mkdir(
                    parents=True,
                    exist_ok=True
                )

                continue

            target.parent.mkdir(
                parents=True,
                exist_ok=True
            )

            src = tar.extractfile(
                member
            )

            if src is None:

                raise RuntimeError(
                    "Cannot extract TAR member: "
                    + member.name
                )

            with src:

                with open(
                    target,
                    "wb"
                ) as dst:

                    shutil.copyfileobj(
                        src,
                        dst
                    )


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

        path = Path(
            raw
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

                expected_i = int(
                    float(
                        expected
                    )
                )

                actual_i = (
                    path.stat().st_size
                )

                if (
                    expected_i
                    != actual_i
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


if (
    gates[
        "technical_preprocessing"
    ]
    != "YES"
):

    fail(
        "Technical preprocessing gate is closed."
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
    "✓ Technical extraction/alignment allowed."
)

print(
    "✓ CORE inference remains blocked."
)

print()


# ============================================================
# 2. Frozen-project integrity
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
# 3. Combined technical inventory
# ============================================================

print("===== 3. LOAD COMBINED 27-TUMOR INVENTORY =====")
print()

combined = read_csv(
    COMBINED_INV
)


if len(
    combined
) != 27:

    fail(
        "Expected 27 rows in combined tumor inventory."
    )


remaining = [
    x
    for x in combined
    if clean(
        x.get(
            "technical_set"
        )
    )
    == "REMAINING_15"
]


paired = [
    x
    for x in combined
    if clean(
        x.get(
            "technical_set"
        )
    )
    == "PAIRED_12"
]


print(
    "Combined tumor observations:",
    len(
        combined
    )
)

print(
    "Paired technical tumors:",
    len(
        paired
    )
)

print(
    "Remaining tumors:",
    len(
        remaining
    )
)

print()


if (
    len(
        paired
    )
    != 12
    or
    len(
        remaining
    )
    != 15
):

    fail(
        "Expected 12 paired + 15 remaining."
    )


# ============================================================
# 4. Extraction-size preflight
# ============================================================

print("===== 4. EXTRACTION SIZE PREFLIGHT =====")
print()

preflight = []

total_expand = 0


for row in remaining:

    filename = clean(
        row.get(
            "filename"
        )
    )

    tar_path = (
        REMAINING_TAR_DIR
        / filename
    )

    if not tar_path.exists():

        fail(
            "Missing downloaded TAR:\n"
            + str(
                tar_path
            )
        )

    with tarfile.open(
        tar_path,
        "r:*"
    ) as tar:

        file_members = [
            m
            for m in tar.getmembers()
            if m.isfile()
        ]

        expanded = sum(
            m.size
            for m in file_members
        )

    total_expand += expanded

    preflight.append({
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

        "sample_name":
            clean(
                row.get(
                    "sample_name"
                )
            ),

        "hrs_accession":
            clean(
                row.get(
                    "hrs_accession"
                )
            ),

        "filename":
            filename,

        "tar_size_MiB":
            round(
                tar_path.stat().st_size
                / 1024**2,
                3
            ),

        "estimated_extracted_MiB":
            round(
                expanded
                / 1024**2,
                3
            ),

        "internal_file_count":
            len(
                file_members
            )
    })


write_csv(
    PREFLIGHT_OUT,
    preflight,
    list(
        preflight[0].keys()
    )
)


free_gib = (
    shutil.disk_usage(
        PROJECT
    ).free
    / 1024**3
)

expand_gib = (
    total_expand
    / 1024**3
)

projected_free = (
    free_gib
    -
    expand_gib
)


print(
    "Estimated extraction size:",
    f"{expand_gib:.3f} GiB"
)

print(
    "Current free disk:",
    f"{free_gib:.2f} GiB"
)

print(
    "Projected free disk:",
    f"{projected_free:.2f} GiB"
)

print()


# Preserve at least 15 GiB free after extraction.
if projected_free < 15:

    fail(
        "Extraction disk guard triggered: "
        "projected free space <15 GiB."
    )


print(
    "✓ Disk guard passed."
)

print()


# ============================================================
# 5. Safe extraction of remaining 15
# ============================================================

print("===== 5. SAFE EXTRACTION OF REMAINING 15 =====")
print()


for i, row in enumerate(
    remaining,
    1
):

    hrs = clean(
        row.get(
            "hrs_accession"
        )
    )

    sample = clean(
        row.get(
            "sample_name"
        )
    )

    filename = clean(
        row.get(
            "filename"
        )
    )

    tar_path = (
        REMAINING_TAR_DIR
        / filename
    )

    sample_dir = (
        EXTRACT_ROOT
        / hrs
    )

    marker = (
        sample_dir
        / ".STEP22A_06B_EXTRACT_COMPLETE"
    )


    print(
        f"[{i:02d}/15] "
        f"{sample:<12s} "
        f"{hrs}"
    )


    tar_sha = sha256_file(
        tar_path
    )


    # Idempotent resume:
    # if marker exists and contains same TAR SHA, retain extraction.
    if marker.exists():

        marker_text = marker.read_text(
            encoding="utf-8",
            errors="replace"
        )

        if tar_sha in marker_text:

            print(
                "  ✓ Existing verified extraction -> reuse."
            )

            continue


    # If directory exists without our verified marker,
    # preserve it rather than silently deleting user data.
    if sample_dir.exists():

        backup = Path(
            str(
                sample_dir
            )
            + ".pre6B_unverified_"
            + datetime.now().strftime(
                "%Y%m%d_%H%M%S"
            )
        )

        print(
            "  Existing unverified directory preserved as:"
        )

        print(
            " ",
            backup.name
        )

        sample_dir.rename(
            backup
        )


    temp_dir = Path(
        str(
            sample_dir
        )
        + ".extracting"
    )


    if temp_dir.exists():

        shutil.rmtree(
            temp_dir
        )


    temp_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    try:

        safe_extract(
            tar_path,
            temp_dir
        )

        matrices, barcodes, genes, _ = (
            find_components(
                temp_dir
            )
        )

        if not (
            len(
                matrices
            )
            == 1
            and
            len(
                barcodes
            )
            == 1
            and
            len(
                genes
            )
            == 1
        ):

            raise RuntimeError(
                "Extracted directory does not contain "
                "exactly one matrix/barcodes/genes set."
            )


        # Atomic-ish promotion after validation.
        temp_dir.rename(
            sample_dir
        )


        marker = (
            sample_dir
            / ".STEP22A_06B_EXTRACT_COMPLETE"
        )


        marker.write_text(
            "\n".join(
                [
                    "STEP22A-6B extraction complete",
                    "source_tar="
                    + filename,
                    "source_sha256="
                    + tar_sha,
                    "recorded_at="
                    + datetime.now().isoformat(
                        timespec="seconds"
                    )
                ]
            )
            + "\n",
            encoding="utf-8"
        )


        print(
            "  ✓ Safe extraction complete."
        )


    except Exception:

        if temp_dir.exists():

            shutil.rmtree(
                temp_dir
            )

        raise


print()

print(
    "✓ Remaining 15 extraction loop complete."
)

print()


# ============================================================
# 6. Matrix structural audit for remaining 15
# ============================================================

print("===== 6. REMAINING-15 MATRIX STRUCTURAL AUDIT =====")
print()

remaining_audit = []

remaining_feature_tables = {}


for row in remaining:

    patient = clean(
        row.get(
            "patient_id"
        )
    )

    timepoint = clean(
        row.get(
            "timepoint"
        )
    )

    sample = clean(
        row.get(
            "sample_name"
        )
    )

    hrs = clean(
        row.get(
            "hrs_accession"
        )
    )

    sample_dir = (
        EXTRACT_ROOT
        / hrs
    )


    matrices, barcodes, genes, files = (
        find_components(
            sample_dir
        )
    )


    if not (
        len(
            matrices
        )
        == 1
        and
        len(
            barcodes
        )
        == 1
        and
        len(
            genes
        )
        == 1
    ):

        fail(
            "Unexpected extracted component count for "
            + sample
        )


    matrix_path = matrices[0]

    barcode_path = barcodes[0]

    gene_path = genes[0]


    nrow, ncol, nnz = (
        read_matrix_dims(
            matrix_path
        )
    )


    barcode_count = count_lines(
        barcode_path
    )


    gene_ids, gene_symbols = (
        read_gene_table(
            gene_path
        )
    )


    gene_count = len(
        gene_ids
    )


    # Barcode uniqueness.
    with open(
        barcode_path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        barcode_values = [
            x.strip()
            for x in f
            if x.strip()
        ]


    unique_barcodes = len(
        set(
            barcode_values
        )
    )


    duplicate_barcodes = (
        len(
            barcode_values
        )
        -
        unique_barcodes
    )


    dimension_status = (
        "PASS"
        if (
            nrow
            == gene_count
            and
            ncol
            == barcode_count
        )
        else "FAIL"
    )


    barcode_status = (
        "PASS"
        if duplicate_barcodes
        == 0
        else "FAIL"
    )


    remaining_feature_tables[
        sample
    ] = {
        "gene_ids":
            gene_ids,

        "gene_symbols":
            gene_symbols,

        "gene_file":
            gene_path
    }


    remaining_audit.append({
        "patient_id":
            patient,

        "timepoint":
            timepoint,

        "sample_name":
            sample,

        "hrs_accession":
            hrs,

        "matrix_file":
            str(
                matrix_path.relative_to(
                    PROJECT
                )
            ),

        "barcodes_file":
            str(
                barcode_path.relative_to(
                    PROJECT
                )
            ),

        "genes_file":
            str(
                gene_path.relative_to(
                    PROJECT
                )
            ),

        "extracted_regular_files":
            len(
                files
            ),

        "matrix_rows":
            nrow,

        "matrix_columns":
            ncol,

        "matrix_nonzero_entries":
            nnz,

        "gene_rows":
            gene_count,

        "barcode_rows":
            barcode_count,

        "unique_barcodes":
            unique_barcodes,

        "duplicate_barcodes":
            duplicate_barcodes,

        "dimension_status":
            dimension_status,

        "barcode_status":
            barcode_status
    })


write_csv(
    REMAINING_MATRIX_OUT,
    remaining_audit,
    list(
        remaining_audit[0].keys()
    )
)


print(
    f"{'SAMPLE':<12}"
    f"{'GENES':>8}"
    f"{'CELLS':>10}"
    f"{'DIM':>8}"
    f"{'BC':>8}"
)

print(
    "-" * 50
)


for x in remaining_audit:

    print(
        f"{x['sample_name']:<12}"
        f"{x['gene_rows']:>8}"
        f"{x['barcode_rows']:>10}"
        f"{x['dimension_status']:>8}"
        f"{x['barcode_status']:>8}"
    )


bad_dimensions = [
    x
    for x in remaining_audit
    if x[
        "dimension_status"
    ]
    != "PASS"
]


bad_barcodes = [
    x
    for x in remaining_audit
    if x[
        "barcode_status"
    ]
    != "PASS"
]


print()

print(
    "Dimension PASS:",
    15
    -
    len(
        bad_dimensions
    ),
    "/ 15"
)

print(
    "Barcode uniqueness PASS:",
    15
    -
    len(
        bad_barcodes
    ),
    "/ 15"
)

print()


if bad_dimensions:

    fail(
        "At least one remaining sample has matrix dimension mismatch."
    )


if bad_barcodes:

    fail(
        "At least one remaining sample has duplicate barcodes."
    )


# ============================================================
# 7. Load existing paired-12 feature tables
# ============================================================

print("===== 7. LOAD ORIGINAL PAIRED-12 FEATURE TABLES =====")
print()


paired_audit = read_csv(
    PAIRED_MATRIX_AUDIT
)


if len(
    paired_audit
) != 12:

    fail(
        "Expected 12 paired matrix-audit rows."
    )


all_feature_tables = {}


feature_count_rows = []


for row in paired_audit:

    sample = clean(
        row.get(
            "sample_name"
        )
    )

    patient = clean(
        row.get(
            "patient_id"
        )
    )

    timepoint = clean(
        row.get(
            "timepoint"
        )
    )

    hrs = clean(
        row.get(
            "hrs_accession"
        )
    )

    gene_path = (
        PROJECT
        / clean(
            row.get(
                "genes_file"
            )
        )
    )


    ids, symbols = read_gene_table(
        gene_path
    )


    all_feature_tables[
        sample
    ] = {
        "gene_ids":
            ids,

        "gene_symbols":
            symbols,

        "gene_file":
            gene_path
    }


    feature_count_rows.append({
        "technical_set":
            "PAIRED_12",

        "patient_id":
            patient,

        "timepoint":
            timepoint,

        "sample_name":
            sample,

        "hrs_accession":
            hrs,

        "feature_rows":
            len(
                ids
            ),

        "unique_ENSG":
            len(
                set(
                    ids
                )
            ),

        "genes_file":
            str(
                gene_path.relative_to(
                    PROJECT
                )
            )
    })


for row in remaining_audit:

    sample = row[
        "sample_name"
    ]

    table = remaining_feature_tables[
        sample
    ]

    all_feature_tables[
        sample
    ] = table


    feature_count_rows.append({
        "technical_set":
            "REMAINING_15",

        "patient_id":
            row[
                "patient_id"
            ],

        "timepoint":
            row[
                "timepoint"
            ],

        "sample_name":
            sample,

        "hrs_accession":
            row[
                "hrs_accession"
            ],

        "feature_rows":
            len(
                table[
                    "gene_ids"
                ]
            ),

        "unique_ENSG":
            len(
                set(
                    table[
                        "gene_ids"
                    ]
                )
            ),

        "genes_file":
            str(
                table[
                    "gene_file"
                ].relative_to(
                    PROJECT
                )
            )
    })


if len(
    all_feature_tables
) != 27:

    fail(
        "Expected feature tables for exactly 27 tumor observations."
    )


write_csv(
    ALL27_FEATURE_COUNT_OUT,
    feature_count_rows,
    [
        "technical_set",
        "patient_id",
        "timepoint",
        "sample_name",
        "hrs_accession",
        "feature_rows",
        "unique_ENSG",
        "genes_file"
    ]
)


print(
    "Paired feature tables:",
    12
)

print(
    "Remaining feature tables:",
    15
)

print(
    "Total feature tables:",
    len(
        all_feature_tables
    )
)

print()


# ============================================================
# 8. Verify frozen old 12-sample common universe
# ============================================================

print("===== 8. VERIFY FROZEN 18,354-FEATURE UNIVERSE =====")
print()


old_common_rows = read_csv(
    OLD_COMMON_FEATURES
)


old_common_ids = [
    clean(
        x.get(
            "gene_id"
        )
    )
    for x in old_common_rows
]


old_common_symbols = [
    clean(
        x.get(
            "gene_symbol"
        )
    )
    for x in old_common_rows
]


if len(
    old_common_ids
) != 18354:

    fail(
        "Expected frozen paired-12 common feature universe = 18,354."
    )


if len(
    set(
        old_common_ids
    )
) != 18354:

    fail(
        "Duplicate ENSG IDs in frozen 18,354 feature universe."
    )


paired_sets = [
    set(
        all_feature_tables[
            clean(
                row.get(
                    "sample_name"
                )
            )
        ][
            "gene_ids"
        ]
    )
    for row in paired_audit
]


recomputed_old = set.intersection(
    *paired_sets
)


print(
    "Frozen paired-12 common ENSG:",
    len(
        old_common_ids
    )
)

print(
    "Recomputed paired-12 intersection:",
    len(
        recomputed_old
    )
)

print()


if set(
    old_common_ids
) != recomputed_old:

    fail(
        "Frozen 18,354 ENSG universe does not match "
        "recomputed paired-12 intersection."
    )


print(
    "✓ Frozen 18,354 feature universe independently verified."
)

print()


# ============================================================
# 9. Compute TRUE all-27 ENSG intersection
# ============================================================

print("===== 9. COMPUTE TRUE ALL-27 ENSG INTERSECTION =====")
print()


all_sets = [
    set(
        x[
            "gene_ids"
        ]
    )
    for x in all_feature_tables.values()
]


all27_common_set = set.intersection(
    *all_sets
)


# Preserve the OLD frozen order.
# Since all27 is an intersection including the original 12,
# it must be a subset of the old 18,354 universe.
all27_common_ids = [
    gid
    for gid in old_common_ids
    if gid in all27_common_set
]


if len(
    all27_common_ids
) != len(
    all27_common_set
):

    unexpected = (
        all27_common_set
        -
        set(
            old_common_ids
        )
    )

    fail(
        "All-27 common universe contains IDs outside "
        "the frozen paired-12 universe, which is impossible.\n"
        + "Unexpected examples: "
        + ", ".join(
            sorted(
                unexpected
            )[:20]
        )
    )


old_symbol_lookup = {
    gid:
        symbol
    for gid, symbol
    in zip(
        old_common_ids,
        old_common_symbols
    )
}


all27_rows = []


for i, gid in enumerate(
    all27_common_ids,
    1
):

    all27_rows.append({
        "all27_common_index":
            i,

        "gene_id":
            gid,

        "gene_symbol":
            old_symbol_lookup.get(
                gid,
                gid
            ),

        "source_order":
            "FROZEN_PAIRED12_COMMON_ORDER"
    })


write_csv(
    ALL27_COMMON_OUT,
    all27_rows,
    [
        "all27_common_index",
        "gene_id",
        "gene_symbol",
        "source_order"
    ]
)


shrink_n = (
    len(
        old_common_ids
    )
    -
    len(
        all27_common_ids
    )
)


print(
    "Paired-12 common ENSG :",
    len(
        old_common_ids
    )
)

print(
    "All-27 common ENSG    :",
    len(
        all27_common_ids
    )
)

print(
    "Features lost after adding 15 tumors:",
    shrink_n
)

print()


if len(
    all27_common_ids
) < 15000:

    fail(
        "All-27 common feature universe <15,000 genes. "
        "Stop for technical review."
    )


print(
    "✓ All-27 common universe computed without union/zero filling."
)

print()


# ============================================================
# 10. Per-sample alignment manifest
#
# No matrices are rewritten.
# We only verify that every sample contains every all-27 ID.
# ============================================================

print("===== 10. ALL-27 ALIGNMENT MANIFEST =====")
print()


alignment_rows = []


combined_lookup = {
    clean(
        x.get(
            "sample_name"
        )
    ):
        x
    for x in combined
}


for sample in sorted(
    all_feature_tables
):

    table = all_feature_tables[
        sample
    ]

    index_lookup = {
        gid:
            i + 1
        for i, gid
        in enumerate(
            table[
                "gene_ids"
            ]
        )
    }


    missing = [
        gid
        for gid in all27_common_ids
        if gid not in index_lookup
    ]


    if missing:

        fail(
            "All-27 common IDs unexpectedly missing from "
            + sample
        )


    mapped_indices = [
        index_lookup[
            gid
        ]
        for gid in all27_common_ids
    ]


    inv = combined_lookup[
        sample
    ]


    alignment_rows.append({
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

        "hrs_accession":
            clean(
                inv.get(
                    "hrs_accession"
                )
            ),

        "original_features":
            len(
                table[
                    "gene_ids"
                ]
            ),

        "all27_common_features":
            len(
                all27_common_ids
            ),

        "all_common_ids_found":
            "YES",

        "row_reordering_required":
            (
                "NO"
                if mapped_indices
                == sorted(
                    mapped_indices
                )
                else "YES"
            ),

        "union_zero_fill_used":
            "NO",

        "aligned_matrix_written":
            "NO",

        "alignment_rule":
            "SUBSET_AND_REORDER_BY_ENSG_AT_ANALYSIS_TIME"
    })


write_csv(
    ALIGNMENT_OUT,
    alignment_rows,
    list(
        alignment_rows[0].keys()
    )
)


print(
    "Alignment manifests:",
    len(
        alignment_rows
    )
)

print(
    "All common IDs present:",
    sum(
        x[
            "all_common_ids_found"
        ]
        == "YES"
        for x in alignment_rows
    ),
    "/ 27"
)

print(
    "Aligned matrix copies written: 0"
)

print(
    "Union + zero filling used: NO"
)

print()


# ============================================================
# 11. Frozen marker coverage
#
# Presence only. No expression values.
# ============================================================

print("===== 11. FROZEN MARKER COVERAGE =====")
print()


marker_groups = {

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


symbol_upper = {
    clean(
        x[
            "gene_symbol"
        ]
    ).upper()
    for x in all27_rows
}


marker_rows = []


for group, genes in marker_groups.items():

    for gene in genes:

        marker_rows.append({
            "marker_class":
                "BROAD_ANNOTATION",

            "marker_group":
                group,

            "gene":
                gene,

            "present_in_all27_common":
                (
                    "YES"
                    if gene.upper()
                    in symbol_upper
                    else "NO"
                )
        })


for gene in core_genes:

    marker_rows.append({
        "marker_class":
            "CORE_PRESENCE_ONLY",

        "marker_group":
            "CORE",

        "gene":
            gene,

        "present_in_all27_common":
            (
                "YES"
                if gene.upper()
                in symbol_upper
                else "NO"
            )
    })


write_csv(
    MARKER_OUT,
    marker_rows,
    [
        "marker_class",
        "marker_group",
        "gene",
        "present_in_all27_common"
    ]
)


broad_rows = [
    x
    for x in marker_rows
    if x[
        "marker_class"
    ]
    == "BROAD_ANNOTATION"
]


core_rows = [
    x
    for x in marker_rows
    if x[
        "marker_class"
    ]
    == "CORE_PRESENCE_ONLY"
]


broad_present = sum(
    x[
        "present_in_all27_common"
    ]
    == "YES"
    for x in broad_rows
)


core_present = sum(
    x[
        "present_in_all27_common"
    ]
    == "YES"
    for x in core_rows
)


print(
    "Frozen broad markers present:",
    broad_present,
    "/",
    len(
        broad_rows
    )
)

print(
    "CORE genes technically present:",
    core_present,
    "/ 4"
)

print()


missing_broad = [
    x[
        "gene"
    ]
    for x in broad_rows
    if x[
        "present_in_all27_common"
    ]
    != "YES"
]


missing_core = [
    x[
        "gene"
    ]
    for x in core_rows
    if x[
        "present_in_all27_common"
    ]
    != "YES"
]


if missing_broad:

    fail(
        "Frozen broad marker(s) missing from all-27 common universe:\n"
        + ", ".join(
            missing_broad
        )
    )


if missing_core:

    fail(
        "CORE gene(s) missing from all-27 common universe:\n"
        + ", ".join(
            missing_core
        )
    )


print(
    "✓ 33/33 frozen broad markers retained."
)

print(
    "✓ CDKN1B / GSN / BGN / TIMP1 retained."
)

print(
    "✓ CORE expression values were NOT read."
)

print()


# ============================================================
# 12. Final disk / frozen checks
# ============================================================

print("===== 12. FINAL SAFETY CHECK =====")
print()


free_after = (
    shutil.disk_usage(
        PROJECT
    ).free
    / 1024**3
)


print(
    "Free disk:",
    f"{free_after:.2f} GiB"
)


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
# 13. Summary
# ============================================================

with open(
    SUMMARY_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-6B SAFE EXTRACTION + ALL27 FEATURE ALIGNMENT\n"
    )

    f.write(
        "=" * 70
        + "\n\n"
    )

    f.write(
        "Remaining tumors safely extracted: 15/15\n"
    )

    f.write(
        "Remaining matrix dimensions PASS: 15/15\n"
    )

    f.write(
        "Remaining barcode uniqueness PASS: 15/15\n"
    )

    f.write(
        f"Paired-12 common ENSG: {len(old_common_ids)}\n"
    )

    f.write(
        f"All-27 common ENSG: {len(all27_common_ids)}\n"
    )

    f.write(
        f"Features lost after expanding 12->27: {shrink_n}\n"
    )

    f.write(
        "Union + zero filling: NO\n"
    )

    f.write(
        "Aligned matrix copies written: NO\n"
    )

    f.write(
        "Frozen broad markers retained: 33/33\n"
    )

    f.write(
        "CORE genes retained: 4/4\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "Matrix expression values inspected: NO\n"
    )

    f.write(
        "CORE expression inspected: NO\n"
    )

    f.write(
        "Cell QC performed: NO\n"
    )

    f.write(
        "Cell annotation performed: NO\n"
    )

    f.write(
        "CORE2 scoring performed: NO\n"
    )

    f.write(
        "Before/After inference performed: NO\n"
    )

    f.write(
        "Histology mapping: INCOMPLETE\n"
    )

    f.write(
        "Primary CORE2 inference: BLOCKED\n"
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

print("=" * 80)
print("STEP22A-6B COMPLETE")
print("=" * 80)
print()

print(
    "REMAINING TUMORS EXTRACTED   : 15 / 15"
)

print(
    "MATRIX DIMENSION PASS        : 15 / 15"
)

print(
    "BARCODE UNIQUENESS PASS      : 15 / 15"
)

print()

print(
    "PAIRED-12 COMMON ENSG        :",
    len(
        old_common_ids
    )
)

print(
    "TRUE ALL-27 COMMON ENSG      :",
    len(
        all27_common_ids
    )
)

print(
    "FEATURES LOST 12 -> 27       :",
    shrink_n
)

print(
    "UNION + ZERO FILLING         : NO"
)

print(
    "ALIGNED MATRIX COPIES WRITTEN: NO"
)

print()

print(
    "FROZEN BROAD MARKERS         :",
    broad_present,
    "/ 33"
)

print(
    "CORE GENES TECHNICALLY PRESENT:",
    core_present,
    "/ 4"
)

print()

print(
    "MATRIX EXPRESSION VALUES READ: NO"
)

print(
    "CORE GENE VALUES READ        : NO"
)

print(
    "CELL QC                      : NOT PERFORMED"
)

print(
    "CELL ANNOTATION              : NOT PERFORMED"
)

print(
    "CORE2 SCORING                : NOT PERFORMED"
)

print(
    "BEFORE/AFTER INFERENCE       : NOT PERFORMED"
)

print()

print(
    "HISTOLOGY MAPPING            : INCOMPLETE"
)

print(
    "EASC IDENTITY                : UNRESOLVED"
)

print(
    "PRIMARY CORE2 INFERENCE      : BLOCKED"
)

print()

print(
    "FREE DISK                    :",
    f"{free_after:.2f} GiB"
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
    "STEP22A-6C = APPLY THE SAME FROZEN"
)

print(
    "SAMPLE-SPECIFIC QC + CONSTRAINED BROAD"
)

print(
    "CELL-TYPE ANNOTATION TO THE REMAINING 15."
)

print()

print(
    "THEN BUILD THE FULL 27-TUMOR"
)

print(
    "FIBROBLAST PSEUDOBULK REFERENCE."
)

print()

print("DO NOT SCORE CORE2.")
print()
