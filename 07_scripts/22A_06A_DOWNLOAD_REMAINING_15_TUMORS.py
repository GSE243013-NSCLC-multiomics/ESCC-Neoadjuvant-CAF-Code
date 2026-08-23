#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-6A
#
# DOWNLOAD REMAINING 15 TUMOR TECHNICAL FILES
#
# This step:
#   - reads the FROZEN target list from Step22A-5G
#   - downloads exactly 15 remaining tumor TARs
#   - supports .part resume
#   - checks size
#   - computes MD5 + SHA256
#   - checks TAR readability and path safety
#   - verifies 10X matrix.mtx / barcodes.tsv / genes|features.tsv
#   - creates combined 27-tumor technical inventory
#
# This step DOES NOT:
#   - extract TAR files
#   - perform QC
#   - perform annotation
#   - read CORE expression
#   - calculate CORE2 / CORE4
#   - perform treatment inference
#
# Histology remains unresolved.
# ============================================================

from pathlib import Path, PurePosixPath
from urllib.parse import urlparse
from datetime import datetime
import csv
import hashlib
import os
import shutil
import subprocess
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

PAIRED_DIR = (
    RAW_ROOT
    / "paired_tumor"
)

REMAINING_DIR = (
    RAW_ROOT
    / "remaining_tumor"
)

LOG_DIR = (
    PROJECT
    / "08_logs"
)

TARGET_FILE = (
    META_DIR
    / "STEP22A_05G_REMAINING_15_TUMOR_TECHNICAL_TARGETS.tsv"
)

PAIRED_TARGET_FILE = (
    META_DIR
    / "STEP22A_03A_PAIRED_DOWNLOAD_TARGETS.tsv"
)

GATE_FILE = (
    META_DIR
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INV = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

STATUS_OUT = (
    META_DIR
    / "STEP22A_06A_DOWNLOAD_STATUS.csv"
)

PROVENANCE_OUT = (
    META_DIR
    / "STEP22A_06A_REMAINING_TUMOR_DOWNLOAD_PROVENANCE.csv"
)

TAR_AUDIT_OUT = (
    META_DIR
    / "STEP22A_06A_REMAINING_TUMOR_TAR_AUDIT.csv"
)

COMBINED_OUT = (
    META_DIR
    / "STEP22A_06A_COMBINED_27_TUMOR_TECHNICAL_INVENTORY.csv"
)

SUMMARY_OUT = (
    META_DIR
    / "STEP22A_06A_DOWNLOAD_SUMMARY.txt"
)


REMAINING_DIR.mkdir(
    parents=True,
    exist_ok=True
)

LOG_DIR.mkdir(
    parents=True,
    exist_ok=True
)


print()
print("=" * 80)
print("STEP22A-6A: DOWNLOAD REMAINING 15 TUMOR TECHNICAL FILES")
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


def read_csv(path, delimiter=","):

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
    fields,
    delimiter=","
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
            delimiter=delimiter,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def clean(x):
    return (
        str(x or "")
        .replace("\r", "")
        .strip()
    )


def to_float(x):

    try:
        return float(
            clean(x)
        )

    except Exception:
        return None


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


def md5_file(path):

    h = hashlib.md5()

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


def quick_tar_readable(path):

    try:

        with tarfile.open(
            path,
            "r:*"
        ) as tar:

            # Force member table read.
            tar.getmembers()

        return True

    except Exception:
        return False


def audit_tar(path):

    result = {
        "tar_readable": "NO",
        "safe_paths": "NO",
        "regular_file_count": 0,
        "matrix_count": 0,
        "barcodes_count": 0,
        "genes_features_count": 0,
        "matrix_files": "",
        "barcodes_files": "",
        "genes_features_files": "",
        "format_status": "FAILED",
        "unsafe_members": ""
    }

    try:

        with tarfile.open(
            path,
            "r:*"
        ) as tar:

            members = tar.getmembers()

            result[
                "tar_readable"
            ] = "YES"

            unsafe = []
            regular = []

            for member in members:

                p = PurePosixPath(
                    member.name
                )

                if (
                    p.is_absolute()
                    or ".." in p.parts
                ):

                    unsafe.append(
                        member.name
                    )

                    continue

                if (
                    member.issym()
                    or member.islnk()
                ):

                    unsafe.append(
                        member.name
                    )

                    continue

                if member.isfile():
                    regular.append(
                        member.name
                    )

            result[
                "regular_file_count"
            ] = len(
                regular
            )

            if unsafe:

                result[
                    "unsafe_members"
                ] = " | ".join(
                    unsafe[:50]
                )

                result[
                    "safe_paths"
                ] = "NO"

                result[
                    "format_status"
                ] = (
                    "UNSAFE_TAR_MEMBERS"
                )

                return result

            result[
                "safe_paths"
            ] = "YES"

            matrices = []

            barcodes = []

            genes = []

            for name in regular:

                base = (
                    PurePosixPath(
                        name
                    )
                    .name
                    .lower()
                )

                if base == "matrix.mtx":
                    matrices.append(
                        name
                    )

                elif base == "barcodes.tsv":
                    barcodes.append(
                        name
                    )

                elif base in {
                    "genes.tsv",
                    "features.tsv"
                }:

                    genes.append(
                        name
                    )

            result[
                "matrix_count"
            ] = len(
                matrices
            )

            result[
                "barcodes_count"
            ] = len(
                barcodes
            )

            result[
                "genes_features_count"
            ] = len(
                genes
            )

            result[
                "matrix_files"
            ] = " | ".join(
                matrices
            )

            result[
                "barcodes_files"
            ] = " | ".join(
                barcodes
            )

            result[
                "genes_features_files"
            ] = " | ".join(
                genes
            )

            if (
                len(matrices) == 1
                and
                len(barcodes) == 1
                and
                len(genes) == 1
            ):

                result[
                    "format_status"
                ] = (
                    "PASS_10X_THREE_FILE"
                )

            else:

                result[
                    "format_status"
                ] = (
                    "REVIEW_NONSTANDARD_10X"
                )

    except Exception as e:

        result[
            "format_status"
        ] = (
            "TAR_READ_FAILED: "
            + str(e)
        )

    return result


def frozen_integrity_check():

    if not FROZEN_INV.exists():
        fail(
            "Missing frozen inventory:\n"
            + str(
                FROZEN_INV
            )
        )

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

        p = Path(
            raw_path
        )

        if not p.is_absolute():
            p = PROJECT / p

        if not p.exists():

            changes.append(
                (
                    str(p),
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
                    p.stat().st_size
                )

                if (
                    actual_i
                    != expected_i
                ):

                    changes.append(
                        (
                            str(p),
                            "SIZE_CHANGED"
                        )
                    )

            except Exception:
                pass

    return changes


# ============================================================
# 1. Analysis-gate check
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
        "Technical preprocessing is blocked."
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
    "✓ Technical download allowed."
)

print(
    "✓ Histology/CORE inference gate remains closed."
)

print()


# ============================================================
# 2. Frozen-project check
# ============================================================

print("===== 2. FROZEN PROJECT CHECK =====")
print()

changes = frozen_integrity_check()

print(
    "Frozen original files changed:",
    len(
        changes
    )
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
# 3. Read remaining-15 target list
# ============================================================

print("===== 3. READ FROZEN 15-TUMOR TARGET LIST =====")
print()

targets_raw = read_csv(
    TARGET_FILE,
    delimiter="\t"
)


required_columns = {
    "patient_id",
    "timepoint",
    "sample_name",
    "hrs_accession",
    "omix_file_id",
    "listed_size_MB",
    "download_url_https"
}


if not targets_raw:

    fail(
        "Remaining target list is empty."
    )


missing_columns = (
    required_columns
    -
    set(
        targets_raw[0].keys()
    )
)


if missing_columns:

    fail(
        "Target TSV missing columns: "
        + ", ".join(
            sorted(
                missing_columns
            )
        )
    )


targets = []


for row in targets_raw:

    x = {
        key:
            clean(
                value
            )
        for key, value
        in row.items()
    }

    x[
        "filename"
    ] = (
        x[
            "omix_file_id"
        ]
        + ".tar"
    )

    targets.append(
        x
    )


print(
    "Remaining target rows:",
    len(
        targets
    )
)


if len(targets) != 15:

    fail(
        "Expected exactly 15 remaining tumor targets."
    )


sample_names = [
    x[
        "sample_name"
    ]
    for x in targets
]


file_ids = [
    x[
        "omix_file_id"
    ]
    for x in targets
]


if (
    len(
        set(
            sample_names
        )
    )
    != 15
):

    fail(
        "Duplicate sample_name in remaining target list."
    )


if (
    len(
        set(
            file_ids
        )
    )
    != 15
):

    fail(
        "Duplicate OMIX file ID in remaining target list."
    )


for x in targets:

    if not x[
        "patient_id"
    ]:

        fail(
            "Missing patient_id."
        )

    if x[
        "timepoint"
    ] not in {
        "Before",
        "After"
    }:

        fail(
            "Unexpected timepoint for "
            + x[
                "sample_name"
            ]
            + ": "
            + x[
                "timepoint"
            ]
        )

    parsed = urlparse(
        x[
            "download_url_https"
        ]
    )

    if (
        parsed.scheme != "https"
        or
        parsed.hostname
        != "download.cncb.ac.cn"
    ):

        fail(
            "Unexpected download URL/host:\n"
            + x[
                "download_url_https"
            ]
        )


print(
    "✓ Exactly 15 unique technical tumor targets."
)

print(
    "✓ All download URLs point to download.cncb.ac.cn."
)

print()


# ============================================================
# 4. Cross-check existing 12 paired targets
# ============================================================

print("===== 4. CROSS-CHECK EXISTING 12 PAIRED TUMORS =====")
print()


paired = read_csv(
    PAIRED_TARGET_FILE,
    delimiter="\t"
)


if len(paired) != 12:

    fail(
        "Expected exactly 12 prior paired target records."
    )


paired_samples = {
    clean(
        x.get(
            "sample_name"
        )
    )
    for x in paired
}


paired_ids = {
    clean(
        x.get(
            "file_id"
        )
    )
    for x in paired
}


remaining_samples = set(
    sample_names
)

remaining_ids = set(
    file_ids
)


if (
    paired_samples
    &
    remaining_samples
):

    fail(
        "Paired and remaining sample sets overlap."
    )


if (
    paired_ids
    &
    remaining_ids
):

    fail(
        "Paired and remaining OMIX file IDs overlap."
    )


if (
    len(
        paired_samples
        |
        remaining_samples
    )
    != 27
):

    fail(
        "Combined tumor sample count is not 27."
    )


paired_tar_files = sorted(
    PAIRED_DIR.glob(
        "*.tar"
    )
)


print(
    "Paired targets:",
    len(
        paired
    )
)

print(
    "Paired TARs currently present:",
    len(
        paired_tar_files
    )
)

print(
    "Combined technical tumor targets:",
    len(
        paired_samples
        |
        remaining_samples
    )
)


if (
    len(
        paired_tar_files
    )
    != 12
):

    fail(
        "Expected 12 paired TAR files before downloading remaining data."
    )


for p in paired_tar_files:

    if not quick_tar_readable(
        p
    ):

        fail(
            "Existing paired TAR is unreadable:\n"
            + str(p)
        )


print(
    "✓ Existing 12 paired TARs remain readable."
)

print(
    "✓ 12 + 15 = 27 unique tumor technical targets."
)

print()


# ============================================================
# 5. Disk-space preflight
# ============================================================

print("===== 5. DISK-SPACE PREFLIGHT =====")
print()


expected_mb = sum(
    x
    for x in [
        to_float(
            row[
                "listed_size_MB"
            ]
        )
        for row in targets
    ]
    if x is not None
)


expected_gib = (
    expected_mb
    / 1024.0
)


disk = shutil.disk_usage(
    PROJECT
)


free_gib = (
    disk.free
    / 1024**3
)


projected_free = (
    free_gib
    -
    expected_gib
)


print(
    "Listed remaining download:",
    f"{expected_mb:.2f} MB"
)

print(
    "Listed remaining download:",
    f"{expected_gib:.3f} GiB"
)

print(
    "Current free disk:",
    f"{free_gib:.2f} GiB"
)

print(
    "Projected free after download:",
    f"{projected_free:.2f} GiB"
)

print()


# Keep 10 GiB reserve after download.
if (
    projected_free
    < 10.0
):

    fail(
        "Disk guard triggered: download would leave <10 GiB free."
    )


print(
    "✓ Disk-space guard passed."
)

print()


# ============================================================
# 6. Download with resumable .part files
# ============================================================

print("===== 6. DOWNLOAD / RESUME 15 TAR FILES =====")
print()


curl_path = shutil.which(
    "curl"
)


if not curl_path:

    fail(
        "curl is not installed or not in PATH."
    )


status_rows = []


def save_status():

    fields = [
        "patient_id",
        "timepoint",
        "sample_name",
        "hrs_accession",
        "omix_file_id",
        "filename",
        "local_path",
        "download_status",
        "partial_bytes",
        "recorded_at"
    ]

    write_csv(
        STATUS_OUT,
        status_rows,
        fields
    )


for i, x in enumerate(
    targets,
    1
):

    final = (
        REMAINING_DIR
        / x[
            "filename"
        ]
    )

    partial = Path(
        str(final)
        + ".part"
    )


    print()
    print("-" * 80)

    print(
        f"[{i:02d}/15] "
        f"{x['patient_id']} "
        f"{x['timepoint']} "
        f"{x['sample_name']}"
    )

    print(
        "HRS   :",
        x[
            "hrs_accession"
        ]
    )

    print(
        "FILE  :",
        x[
            "filename"
        ]
    )

    print(
        "SIZE  :",
        x[
            "listed_size_MB"
        ],
        "MB"
    )


    # --------------------------------------------------------
    # Existing final file:
    # keep it if readable; otherwise preserve and redownload.
    # --------------------------------------------------------

    if final.exists():

        if quick_tar_readable(
            final
        ):

            print(
                "STATUS: existing readable TAR -> SKIP DOWNLOAD"
            )

            status_rows.append({
                "patient_id":
                    x[
                        "patient_id"
                    ],

                "timepoint":
                    x[
                        "timepoint"
                    ],

                "sample_name":
                    x[
                        "sample_name"
                    ],

                "hrs_accession":
                    x[
                        "hrs_accession"
                    ],

                "omix_file_id":
                    x[
                        "omix_file_id"
                    ],

                "filename":
                    x[
                        "filename"
                    ],

                "local_path":
                    str(
                        final.relative_to(
                            PROJECT
                        )
                    ),

                "download_status":
                    "EXISTING_READABLE_SKIPPED",

                "partial_bytes":
                    (
                        partial.stat().st_size
                        if partial.exists()
                        else 0
                    ),

                "recorded_at":
                    datetime.now().isoformat(
                        timespec="seconds"
                    )
            })

            save_status()

            continue

        else:

            backup = Path(
                str(final)
                + ".corrupt_"
                + datetime.now().strftime(
                    "%Y%m%d_%H%M%S"
                )
            )

            print(
                "Existing final TAR is unreadable."
            )

            print(
                "Preserving as:",
                backup.name
            )

            final.rename(
                backup
            )


    if partial.exists():

        print(
            "RESUME FROM:",
            f"{partial.stat().st_size / 1024**2:.2f} MiB"
        )

    else:

        print(
            "START NEW DOWNLOAD"
        )


    cmd = [
        curl_path,
        "-L",
        "--fail",
        "--show-error",
        "--retry",
        "8",
        "--retry-delay",
        "10",
        "--retry-all-errors",
        "--connect-timeout",
        "30",
        "-C",
        "-",
        "-o",
        str(
            partial
        ),
        x[
            "download_url_https"
        ]
    ]


    result = subprocess.run(
        cmd
    )


    if (
        result.returncode
        != 0
    ):

        status_rows.append({
            "patient_id":
                x[
                    "patient_id"
                ],

            "timepoint":
                x[
                    "timepoint"
                ],

            "sample_name":
                x[
                    "sample_name"
                ],

            "hrs_accession":
                x[
                    "hrs_accession"
                ],

            "omix_file_id":
                x[
                    "omix_file_id"
                ],

            "filename":
                x[
                    "filename"
                ],

            "local_path":
                str(
                    final.relative_to(
                        PROJECT
                    )
                ),

            "download_status":
                (
                    "DOWNLOAD_INTERRUPTED_"
                    "PART_PRESERVED"
                ),

            "partial_bytes":
                (
                    partial.stat().st_size
                    if partial.exists()
                    else 0
                ),

            "recorded_at":
                datetime.now().isoformat(
                    timespec="seconds"
                )
        })


        save_status()


        print()
        print(
            "DOWNLOAD INTERRUPTED."
        )

        print(
            "The .part file has been preserved."
        )

        print(
            "Re-run this same command to resume."
        )

        raise SystemExit(
            result.returncode
        )


    if not partial.exists():

        fail(
            "curl returned success but .part file is absent:\n"
            + str(partial)
        )


    os.replace(
        partial,
        final
    )


    print(
        "✓ Download complete:",
        f"{final.stat().st_size / 1024**2:.2f} MiB"
    )


    status_rows.append({
        "patient_id":
            x[
                "patient_id"
            ],

        "timepoint":
            x[
                "timepoint"
            ],

        "sample_name":
            x[
                "sample_name"
            ],

        "hrs_accession":
            x[
                "hrs_accession"
            ],

        "omix_file_id":
            x[
                "omix_file_id"
            ],

        "filename":
            x[
                "filename"
            ],

        "local_path":
            str(
                final.relative_to(
                    PROJECT
                )
            ),

        "download_status":
            "DOWNLOADED_COMPLETE",

        "partial_bytes":
            0,

        "recorded_at":
            datetime.now().isoformat(
                timespec="seconds"
            )
    })


    save_status()


print()
print(
    "✓ Download loop completed."
)

print()


# ============================================================
# 7. Verify exact file count / no .part
# ============================================================

print("===== 7. FILE-COUNT CHECK =====")
print()


remaining_tars = sorted(
    REMAINING_DIR.glob(
        "*.tar"
    )
)

remaining_parts = sorted(
    REMAINING_DIR.glob(
        "*.part"
    )
)


print(
    "Remaining tumor TAR files:",
    len(
        remaining_tars
    )
)

print(
    "Remaining .part files:",
    len(
        remaining_parts
    )
)


if len(
    remaining_tars
) != 15:

    fail(
        "Expected exactly 15 final TARs in remaining_tumor."
    )


if remaining_parts:

    fail(
        "One or more .part files remain."
    )


expected_filenames = {
    x[
        "filename"
    ]
    for x in targets
}


actual_filenames = {
    p.name
    for p in remaining_tars
}


if (
    actual_filenames
    != expected_filenames
):

    print(
        "Missing:",
        sorted(
            expected_filenames
            -
            actual_filenames
        )
    )

    print(
        "Unexpected:",
        sorted(
            actual_filenames
            -
            expected_filenames
        )
    )

    fail(
        "Remaining tumor TAR filenames do not exactly match frozen target list."
    )


print(
    "✓ Exactly the frozen 15 TAR files are present."
)

print()


# ============================================================
# 8. TAR audit + path-safety check
# ============================================================

print("===== 8. TAR READABILITY / SAFETY / 10X AUDIT =====")
print()


tar_rows = []


for i, x in enumerate(
    targets,
    1
):

    path = (
        REMAINING_DIR
        / x[
            "filename"
        ]
    )


    audit = audit_tar(
        path
    )


    row = {
        "patient_id":
            x[
                "patient_id"
            ],

        "timepoint":
            x[
                "timepoint"
            ],

        "sample_name":
            x[
                "sample_name"
            ],

        "hrs_accession":
            x[
                "hrs_accession"
            ],

        "omix_file_id":
            x[
                "omix_file_id"
            ],

        "filename":
            x[
                "filename"
            ],

        **audit
    }


    tar_rows.append(
        row
    )


    print(
        f"[{i:02d}/15] "
        f"{x['sample_name']:<12s} "
        f"{audit['format_status']}"
    )


write_csv(
    TAR_AUDIT_OUT,
    tar_rows,
    [
        "patient_id",
        "timepoint",
        "sample_name",
        "hrs_accession",
        "omix_file_id",
        "filename",
        "tar_readable",
        "safe_paths",
        "regular_file_count",
        "matrix_count",
        "barcodes_count",
        "genes_features_count",
        "matrix_files",
        "barcodes_files",
        "genes_features_files",
        "format_status",
        "unsafe_members"
    ]
)


bad_read = [
    x
    for x in tar_rows
    if x[
        "tar_readable"
    ]
    != "YES"
]


bad_safe = [
    x
    for x in tar_rows
    if x[
        "safe_paths"
    ]
    != "YES"
]


bad_format = [
    x
    for x in tar_rows
    if x[
        "format_status"
    ]
    != "PASS_10X_THREE_FILE"
]


print()

print(
    "Readable TARs:",
    15 - len(
        bad_read
    ),
    "/ 15"
)

print(
    "Safe TAR paths:",
    15 - len(
        bad_safe
    ),
    "/ 15"
)

print(
    "10X three-file archives:",
    15 - len(
        bad_format
    ),
    "/ 15"
)

print()


if bad_read:

    fail(
        "At least one downloaded TAR cannot be read."
    )


if bad_safe:

    fail(
        "At least one TAR contains unsafe paths/links."
    )


if bad_format:

    print(
        "Nonstandard archives:"
    )

    for x in bad_format:

        print(
            " ",
            x[
                "sample_name"
            ],
            "->",
            x[
                "format_status"
            ]
        )

    fail(
        "At least one TAR is not the expected 10X three-file structure. "
        "Stop for review before extraction."
    )


print(
    "✓ All 15 TARs readable."
)

print(
    "✓ All 15 TARs pass path/link safety checks."
)

print(
    "✓ All 15 contain one matrix.mtx + one barcodes.tsv "
    "+ one genes/features.tsv."
)

print()


# ============================================================
# 9. Size + MD5 + SHA256 provenance
# ============================================================

print("===== 9. SIZE + HASH PROVENANCE =====")
print()


provenance = []


for i, x in enumerate(
    targets,
    1
):

    path = (
        REMAINING_DIR
        / x[
            "filename"
        ]
    )


    size_bytes = (
        path.stat().st_size
    )

    actual_mib = (
        size_bytes
        / 1024**2
    )


    expected = to_float(
        x[
            "listed_size_MB"
        ]
    )


    if expected is None:

        size_diff = ""

        size_status = (
            "EXPECTED_SIZE_UNAVAILABLE"
        )

    else:

        diff = (
            actual_mib
            -
            expected
        )

        tolerance = max(
            1.0,
            expected * 0.02
        )

        size_diff = round(
            diff,
            3
        )

        size_status = (
            "PASS"
            if abs(
                diff
            )
            <= tolerance
            else "REVIEW"
        )


    print(
        f"[{i:02d}/15] "
        f"{x['sample_name']}: hashing..."
    )


    md5 = md5_file(
        path
    )

    sha = sha256_file(
        path
    )


    provenance.append({
        "patient_id":
            x[
                "patient_id"
            ],

        "timepoint":
            x[
                "timepoint"
            ],

        "sample_name":
            x[
                "sample_name"
            ],

        "hrs_accession":
            x[
                "hrs_accession"
            ],

        "omix_file_id":
            x[
                "omix_file_id"
            ],

        "filename":
            x[
                "filename"
            ],

        "local_path":
            str(
                path.relative_to(
                    PROJECT
                )
            ),

        "download_url":
            x[
                "download_url_https"
            ],

        "listed_size_MB":
            x[
                "listed_size_MB"
            ],

        "actual_size_bytes":
            size_bytes,

        "actual_size_MiB":
            round(
                actual_mib,
                3
            ),

        "size_difference":
            size_diff,

        "size_status":
            size_status,

        "md5":
            md5,

        "sha256":
            sha,

        "tar_readable":
            "YES",

        "tar_safe":
            "YES",

        "format_status":
            "PASS_10X_THREE_FILE",

        "histology_status":
            (
                "UNRESOLVED_"
                "TECHNICAL_DATA_ONLY"
            ),

        "recorded_at":
            datetime.now().isoformat(
                timespec="seconds"
            )
    })


write_csv(
    PROVENANCE_OUT,
    provenance,
    list(
        provenance[0].keys()
    )
)


size_review = [
    x
    for x in provenance
    if x[
        "size_status"
    ]
    == "REVIEW"
]


total_actual_gib = (
    sum(
        x[
            "actual_size_bytes"
        ]
        for x in provenance
    )
    / 1024**3
)


print()

print(
    "Actual downloaded total:",
    f"{total_actual_gib:.3f} GiB"
)

print(
    "Size checks requiring review:",
    len(
        size_review
    )
)

print()


# ============================================================
# 10. Build combined 27-tumor TECHNICAL inventory
#
# IMPORTANT:
# This is not an ESCC eligibility file.
# Histology remains unresolved.
# ============================================================

print("===== 10. BUILD COMBINED 27-TUMOR TECHNICAL INVENTORY =====")
print()


combined = []


for x in paired:

    sample_name = clean(
        x.get(
            "sample_name"
        )
    )

    filename = clean(
        x.get(
            "filename"
        )
    )

    if not filename:

        filename = (
            clean(
                x.get(
                    "file_id"
                )
            )
            + ".tar"
        )


    local_path = (
        PAIRED_DIR
        / filename
    )


    combined.append({
        "technical_set":
            "PAIRED_12",

        "patient_id":
            clean(
                x.get(
                    "patient_id"
                )
            ),

        "timepoint":
            clean(
                x.get(
                    "timepoint"
                )
            ),

        "sample_name":
            sample_name,

        "hrs_accession":
            clean(
                x.get(
                    "hrs_accession"
                )
            ),

        "omix_file_id":
            clean(
                x.get(
                    "file_id"
                )
            ),

        "filename":
            filename,

        "local_path":
            str(
                local_path.relative_to(
                    PROJECT
                )
            ),

        "local_tar_present":
            (
                "YES"
                if local_path.exists()
                else "NO"
            ),

        "histology_status":
            (
                "UNRESOLVED_"
                "NOT_ESCC_ELIGIBILITY_CLAIM"
            )
    })


for x in targets:

    local_path = (
        REMAINING_DIR
        / x[
            "filename"
        ]
    )


    combined.append({
        "technical_set":
            "REMAINING_15",

        "patient_id":
            x[
                "patient_id"
            ],

        "timepoint":
            x[
                "timepoint"
            ],

        "sample_name":
            x[
                "sample_name"
            ],

        "hrs_accession":
            x[
                "hrs_accession"
            ],

        "omix_file_id":
            x[
                "omix_file_id"
            ],

        "filename":
            x[
                "filename"
            ],

        "local_path":
            str(
                local_path.relative_to(
                    PROJECT
                )
            ),

        "local_tar_present":
            (
                "YES"
                if local_path.exists()
                else "NO"
            ),

        "histology_status":
            (
                "UNRESOLVED_"
                "NOT_ESCC_ELIGIBILITY_CLAIM"
            )
    })


combined.sort(
    key=lambda x: (
        int(
            x[
                "patient_id"
            ].replace(
                "P",
                ""
            )
        ),
        0
        if x[
            "timepoint"
        ]
        == "Before"
        else 1,
        x[
            "sample_name"
        ]
    )
)


if len(
    combined
) != 27:

    fail(
        "Combined technical inventory is not 27 rows."
    )


if any(
    x[
        "local_tar_present"
    ]
    != "YES"
    for x in combined
):

    fail(
        "One or more of the 27 technical tumor TARs are absent."
    )


if (
    len(
        {
            x[
                "sample_name"
            ]
            for x in combined
        }
    )
    != 27
):

    fail(
        "Combined technical inventory has duplicate sample names."
    )


write_csv(
    COMBINED_OUT,
    combined,
    [
        "technical_set",
        "patient_id",
        "timepoint",
        "sample_name",
        "hrs_accession",
        "omix_file_id",
        "filename",
        "local_path",
        "local_tar_present",
        "histology_status"
    ]
)


unique_patients = sorted(
    {
        x[
            "patient_id"
        ]
        for x in combined
    },
    key=lambda x:
        int(
            x.replace(
                "P",
                ""
            )
        )
)


before_n = sum(
    x[
        "timepoint"
    ]
    == "Before"
    for x in combined
)


after_n = sum(
    x[
        "timepoint"
    ]
    == "After"
    for x in combined
)


print(
    "Combined tumor TARs:",
    len(
        combined
    )
)

print(
    "Tumor patients represented:",
    len(
        unique_patients
    )
)

print(
    "Before observations:",
    before_n
)

print(
    "After observations:",
    after_n
)

print()


if (
    before_n != 15
    or after_n != 12
):

    fail(
        "Expected T_B=15 and T_A=12 in the combined 27-tumor inventory."
    )


if len(
    unique_patients
) != 21:

    fail(
        "Expected 21 patients with tumor observations."
    )


print(
    "✓ Combined structure = 27 tumor observations / 21 patients."
)

print(
    "✓ Before=15 / After=12."
)

print()


# ============================================================
# 11. Final disk check
# ============================================================

print("===== 11. FINAL DISK CHECK =====")
print()


disk_after = shutil.disk_usage(
    PROJECT
)


free_after = (
    disk_after.free
    / 1024**3
)


print(
    "Free disk after download:",
    f"{free_after:.2f} GiB"
)

print()


# ============================================================
# 12. Final frozen-project check
# ============================================================

print("===== 12. FINAL FROZEN PROJECT CHECK =====")
print()


changes_after = frozen_integrity_check()


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
# 13. Summary
# ============================================================

with open(
    SUMMARY_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-6A REMAINING TUMOR DOWNLOAD SUMMARY\n"
    )

    f.write(
        "=" * 65
        + "\n\n"
    )

    f.write(
        "Remaining tumor targets: 15\n"
    )

    f.write(
        f"Actual newly available remaining data: "
        f"{total_actual_gib:.3f} GiB\n"
    )

    f.write(
        "Readable remaining TARs: 15/15\n"
    )

    f.write(
        "Safe-path remaining TARs: 15/15\n"
    )

    f.write(
        "10X three-file remaining TARs: 15/15\n"
    )

    f.write(
        f"Size review flags: {len(size_review)}\n"
    )

    f.write(
        "Combined technical tumor TARs: 27\n"
    )

    f.write(
        "Tumor patients represented: 21\n"
    )

    f.write(
        "Before tumor observations: 15\n"
    )

    f.write(
        "After tumor observations: 12\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "TAR extraction performed: NO\n"
    )

    f.write(
        "Cell annotation performed: NO\n"
    )

    f.write(
        "CORE gene values read: NO\n"
    )

    f.write(
        "CORE2 scoring performed: NO\n"
    )

    f.write(
        "CORE4 scoring performed: NO\n"
    )

    f.write(
        "Before/After inference performed: NO\n"
    )

    f.write(
        "Histology mapping: INCOMPLETE\n"
    )

    f.write(
        "EASC identity: UNRESOLVED\n"
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
print("STEP22A-6A COMPLETE")
print("=" * 80)
print()

print(
    "REMAINING TUMOR TARGETS      : 15"
)

print(
    "REMAINING TUMOR TARs         : 15 / 15"
)

print(
    "TAR READABILITY              : 15 / 15 PASS"
)

print(
    "TAR PATH SAFETY              : 15 / 15 PASS"
)

print(
    "10X THREE-FILE STRUCTURE     : 15 / 15 PASS"
)

print(
    "SIZE REVIEW FLAGS            :",
    len(
        size_review
    )
)

print(
    "ACTUAL REMAINING DATA        :",
    f"{total_actual_gib:.3f} GiB"
)

print()

print(
    "COMBINED TECHNICAL TUMORS    : 27"
)

print(
    "TUMOR PATIENTS REPRESENTED   :",
    len(
        unique_patients
    )
)

print(
    "BEFORE TUMOR OBSERVATIONS    :",
    before_n
)

print(
    "AFTER TUMOR OBSERVATIONS     :",
    after_n
)

print()

print(
    "TAR EXTRACTION               : NOT PERFORMED"
)

print(
    "CELL ANNOTATION              : NOT PERFORMED"
)

print(
    "CORE GENE VALUES READ        : NO"
)

print(
    "CORE2 SCORING                : NOT PERFORMED"
)

print(
    "CORE4 SCORING                : NOT PERFORMED"
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
    "STEP22A-6B = SAFE EXTRACT THE REMAINING 15"
)

print(
    "+ ALIGN TO THE FROZEN COMMON ENSG UNIVERSE"
)

print(
    "+ TECHNICAL QC ONLY."
)

print()

print(
    "DO NOT SCORE CORE2."
)

print(
    "DO NOT RUN BEFORE/AFTER STATISTICS."
)

print()
