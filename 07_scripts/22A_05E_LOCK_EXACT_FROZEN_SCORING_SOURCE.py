#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5E
#
# LOCK EXACT FROZEN STEP19L SCORING SOURCE
#
# STATIC SOURCE AUDIT ONLY.
#
# This step DOES NOT:
# - source() or execute Step19L
# - read OMIX CORE-gene expression values
# - normalize OMIX pseudobulk
# - calculate CORE2 / CORE4
# - perform before/after inference
#
# It only:
# - verifies frozen Step19L
# - makes an exact byte-for-byte source copy
# - records SHA256
# - extracts exact cpm()/scale()/rowMeans() calls
# - extracts exact CORE2/CORE4 source context
# - extracts Step19L section 19L3 / 19L8
# - inspects frozen output HEADERS only, never values
# ============================================================

from pathlib import Path
import csv
import hashlib
import os
import re
import shutil
import sys


PROJECT = Path.cwd()

META_DIR = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

STEP19L = (
    PROJECT
    / "07_scripts"
    / "19L_CONSERVED_FIBROBLAST_CORE_VALIDATION.R"
)

GATE_FILE = (
    META_DIR
    / "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

FROZEN_INVENTORY = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

SOURCE_COPY = (
    META_DIR
    / "STEP22A_05E_FROZEN_STEP19L_SOURCE_COPY.R"
)

META_DIR.mkdir(
    parents=True,
    exist_ok=True
)


print()
print("=" * 78)
print("STEP22A-5E: LOCK EXACT FROZEN STEP19L SCORING SOURCE")
print("=" * 78)
print()


# ============================================================
# Helpers
# ============================================================

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


def read_csv(path):

    with open(
        path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        return list(
            csv.DictReader(f)
        )


def strip_r_comments(text):
    """
    Remove R comments while preserving:
    - quoted strings
    - newlines
    - approximate line numbers
    """

    out = []

    in_single = False
    in_double = False
    escaped = False

    i = 0

    while i < len(text):

        c = text[i]

        if escaped:

            out.append(c)
            escaped = False
            i += 1
            continue

        if c == "\\" and (
            in_single
            or in_double
        ):

            out.append(c)
            escaped = True
            i += 1
            continue

        if (
            c == "'"
            and not in_double
        ):

            in_single = not in_single
            out.append(c)
            i += 1
            continue

        if (
            c == '"'
            and not in_single
        ):

            in_double = not in_double
            out.append(c)
            i += 1
            continue

        if (
            c == "#"
            and not in_single
            and not in_double
        ):

            # Preserve columns with spaces until newline.
            while (
                i < len(text)
                and text[i] != "\n"
            ):

                out.append(" ")
                i += 1

            continue

        out.append(c)
        i += 1

    return "".join(out)


def line_number(text, pos):

    return (
        text.count(
            "\n",
            0,
            pos
        )
        + 1
    )


def extract_balanced_calls(
    original_text,
    code_text,
    regex
):
    """
    Find exact function calls in code-only source.
    Parentheses are balanced.
    Positions correspond to original source because
    comments were replaced with spaces.
    """

    calls = []

    for match in re.finditer(
        regex,
        code_text,
        flags=re.I
    ):

        start = match.start()

        open_pos = code_text.find(
            "(",
            match.start(),
            match.end() + 2
        )

        if open_pos < 0:
            continue

        depth = 0
        in_single = False
        in_double = False
        escaped = False

        end = None

        for i in range(
            open_pos,
            len(code_text)
        ):

            c = code_text[i]

            if escaped:

                escaped = False
                continue

            if (
                c == "\\"
                and (
                    in_single
                    or in_double
                )
            ):

                escaped = True
                continue

            if (
                c == "'"
                and not in_double
            ):

                in_single = not in_single
                continue

            if (
                c == '"'
                and not in_single
            ):

                in_double = not in_double
                continue

            if (
                in_single
                or in_double
            ):
                continue

            if c == "(":
                depth += 1

            elif c == ")":

                depth -= 1

                if depth == 0:

                    end = i + 1
                    break

        if end is None:
            continue

        call = original_text[
            start:end
        ]

        calls.append({
            "line_number":
                line_number(
                    original_text,
                    start
                ),

            "call":
                re.sub(
                    r"\s+",
                    " ",
                    call
                ).strip(),

            "raw_call":
                call.strip()
        })

    return calls


def context_around_lines(
    lines,
    hit_lines,
    radius=4
):

    selected = set()

    for line_no in hit_lines:

        for n in range(
            max(
                1,
                line_no - radius
            ),
            min(
                len(lines),
                line_no + radius
            ) + 1
        ):

            selected.add(n)

    return [
        {
            "line_number": n,
            "code": lines[n - 1]
        }
        for n in sorted(selected)
    ]


def find_section(
    lines,
    start_pattern,
    next_pattern
):

    start = None
    end = None

    for i, line in enumerate(
        lines,
        1
    ):

        if (
            start is None
            and re.search(
                start_pattern,
                line,
                flags=re.I
            )
        ):

            start = i
            continue

        if (
            start is not None
            and re.search(
                next_pattern,
                line,
                flags=re.I
            )
        ):

            end = i - 1
            break

    if start is None:
        return None, []

    if end is None:
        end = len(lines)

    return (
        (start, end),
        [
            {
                "line_number": i,
                "code": lines[i - 1]
            }
            for i in range(
                start,
                end + 1
            )
        ]
    )


# ============================================================
# 1. Gate check
# ============================================================

print("===== 1. ANALYSIS GATE CHECK =====")
print()

if not GATE_FILE.exists():

    raise SystemExit(
        "ERROR: missing analysis gate file."
    )

gate_rows = read_csv(
    GATE_FILE
)

gate = {
    x["gate"]:
        x["allowed"]
    for x in gate_rows
}

required_gates = [
    "technical_preprocessing",
    "paired_CORE2_test",
    "manuscript_update"
]

for g in required_gates:

    if g not in gate:

        raise SystemExit(
            "ERROR: missing gate: "
            + g
        )

print(
    "Technical preprocessing :",
    gate[
        "technical_preprocessing"
    ]
)

print(
    "Paired CORE2 inference  :",
    gate[
        "paired_CORE2_test"
    ]
)

print(
    "Manuscript update       :",
    gate[
        "manuscript_update"
    ]
)

print()

if (
    gate[
        "technical_preprocessing"
    ]
    != "YES"
):

    raise SystemExit(
        "TECHNICAL_AUDIT_BLOCKED"
    )

if (
    gate[
        "paired_CORE2_test"
    ]
    == "YES"
):

    raise SystemExit(
        "ERROR: CORE2 inference gate unexpectedly open."
    )

if (
    gate[
        "manuscript_update"
    ]
    == "YES"
):

    raise SystemExit(
        "ERROR: manuscript gate unexpectedly open."
    )

print(
    "✓ Static source audit allowed."
)

print(
    "✓ CORE2 inference remains blocked."
)

print()


# ============================================================
# 2. Verify frozen project before touching anything
# ============================================================

print("===== 2. FROZEN PROJECT INTEGRITY =====")
print()

if not FROZEN_INVENTORY.exists():

    raise SystemExit(
        "ERROR: missing frozen inventory."
    )

frozen = read_csv(
    FROZEN_INVENTORY
)

changes = []

for row in frozen:

    path = Path(
        row.get(
            "path",
            ""
        )
    )

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

    expected_size = row.get(
        "size_bytes",
        ""
    )

    if expected_size:

        try:

            expected_size = int(
                float(
                    expected_size
                )
            )

            actual_size = (
                path.stat().st_size
            )

            if (
                actual_size
                != expected_size
            ):

                changes.append(
                    (
                        str(path),
                        "SIZE_CHANGED"
                    )
                )

        except Exception:
            pass


print(
    "Frozen original files changed:",
    len(changes)
)

if changes:

    for x in changes[:20]:
        print(x)

    raise SystemExit(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )

print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# 3. Freeze exact Step19L source bytes
# ============================================================

print("===== 3. FREEZE EXACT STEP19L SOURCE =====")
print()

if not STEP19L.exists():

    raise SystemExit(
        "ERROR: missing Step19L source:\n"
        + str(
            STEP19L
        )
    )


source_sha = sha256_file(
    STEP19L
)

source_size = (
    STEP19L.stat().st_size
)

# Exact byte-for-byte copy
shutil.copyfile(
    STEP19L,
    SOURCE_COPY
)

copy_sha = sha256_file(
    SOURCE_COPY
)

copy_size = (
    SOURCE_COPY.stat().st_size
)

print(
    "Source:",
    STEP19L
)

print(
    "Source bytes:",
    source_size
)

print(
    "Source SHA256:",
    source_sha
)

print()

print(
    "Frozen copy:",
    SOURCE_COPY
)

print(
    "Copy bytes:",
    copy_size
)

print(
    "Copy SHA256:",
    copy_sha
)

print()

if (
    source_sha != copy_sha
    or
    source_size != copy_size
):

    raise SystemExit(
        "FROZEN_SOURCE_COPY_MISMATCH"
    )

print(
    "✓ Byte-for-byte frozen scoring source locked."
)

print()


# ============================================================
# 4. Read source text
# ============================================================

raw_bytes = STEP19L.read_bytes()

try:

    source_text = raw_bytes.decode(
        "utf-8"
    )

except UnicodeDecodeError:

    source_text = raw_bytes.decode(
        "utf-8",
        errors="replace"
    )

lines = source_text.splitlines()

code_text = strip_r_comments(
    source_text
)

print(
    "Step19L source lines:",
    len(lines)
)

print()


# ============================================================
# 5. Extract frozen specification sections
#
# Step19L3 = sample-level core score
# Step19L8 = strong vs weak two-gene audit
# ============================================================

print("===== 4. EXTRACT FROZEN SCORE SECTIONS =====")
print()

sec3_range, sec3 = find_section(
    lines,
    r"\b19L3\b",
    r"\b19L4\b"
)

sec8_range, sec8 = find_section(
    lines,
    r"\b19L8\b",
    r"\b19L9\b"
)

print(
    "19L3 found:",
    "YES"
    if sec3_range
    else "NO",
    sec3_range
    if sec3_range
    else ""
)

print(
    "19L8 found:",
    "YES"
    if sec8_range
    else "NO",
    sec8_range
    if sec8_range
    else ""
)

print()

if not sec3_range:

    raise SystemExit(
        "STEP19L_SECTION_19L3_NOT_FOUND"
    )

if not sec8_range:

    raise SystemExit(
        "STEP19L_SECTION_19L8_NOT_FOUND"
    )


SECTION_OUT = (
    META_DIR
    / "STEP22A_05E_FROZEN_SCORING_SECTIONS.txt"
)

with open(
    SECTION_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-5E FROZEN STEP19L SCORING SECTIONS\n"
    )

    f.write(
        "SOURCE SHA256: "
        + source_sha
        + "\n\n"
    )

    f.write(
        "========================================\n"
    )

    f.write(
        "SECTION 19L3\n"
    )

    f.write(
        "========================================\n"
    )

    for row in sec3:

        f.write(
            f"{row['line_number']:05d} | "
            f"{row['code']}\n"
        )

    f.write(
        "\n========================================\n"
    )

    f.write(
        "SECTION 19L8\n"
    )

    f.write(
        "========================================\n"
    )

    for row in sec8:

        f.write(
            f"{row['line_number']:05d} | "
            f"{row['code']}\n"
        )


# ============================================================
# 6. Exact cpm() calls
# ============================================================

print("===== 5. EXACT cpm() CALLS =====")
print()

cpm_calls = extract_balanced_calls(
    source_text,
    code_text,
    r"(?<![\w.])(?:edgeR::)?cpm\s*\("
)

print(
    "Executable cpm() calls detected:",
    len(cpm_calls)
)

print()

CPM_OUT = (
    META_DIR
    / "STEP22A_05E_EXACT_CPM_CALLS.csv"
)

cpm_rows = []

for i, call in enumerate(
    cpm_calls,
    1
):

    text = call[
        "call"
    ]

    log_match = re.search(
        r"\blog\s*=\s*(TRUE|FALSE)",
        text,
        flags=re.I
    )

    prior_match = re.search(
        r"\bprior\.count\s*=\s*([^,\)]+)",
        text,
        flags=re.I
    )

    normalized_match = re.search(
        r"\bnormalized\.lib\.sizes\s*=\s*([^,\)]+)",
        text,
        flags=re.I
    )

    row = {
        "call_number":
            i,

        "line_number":
            call[
                "line_number"
            ],

        "log_argument":
            (
                log_match.group(1)
                if log_match
                else "NOT_EXPLICIT"
            ),

        "prior_count_argument":
            (
                prior_match.group(1).strip()
                if prior_match
                else "NOT_EXPLICIT"
            ),

        "normalized_lib_sizes_argument":
            (
                normalized_match.group(1).strip()
                if normalized_match
                else "NOT_EXPLICIT"
            ),

        "exact_call":
            text
    }

    cpm_rows.append(
        row
    )

    print(
        f"CPM CALL {i}"
    )

    print(
        "  line              :",
        call[
            "line_number"
        ]
    )

    print(
        "  log               :",
        row[
            "log_argument"
        ]
    )

    print(
        "  prior.count       :",
        row[
            "prior_count_argument"
        ]
    )

    print(
        "  normalized.lib    :",
        row[
            "normalized_lib_sizes_argument"
        ]
    )

    print(
        "  exact             :",
        text
    )

    print()


write_csv(
    CPM_OUT,
    cpm_rows,
    [
        "call_number",
        "line_number",
        "log_argument",
        "prior_count_argument",
        "normalized_lib_sizes_argument",
        "exact_call"
    ]
)


# ============================================================
# 7. Exact scale() calls
# ============================================================

print("===== 6. EXACT scale() CALLS =====")
print()

scale_calls = extract_balanced_calls(
    source_text,
    code_text,
    r"(?<![\w.])scale\s*\("
)

print(
    "Executable scale() calls detected:",
    len(scale_calls)
)

print()

SCALE_OUT = (
    META_DIR
    / "STEP22A_05E_EXACT_SCALE_CALLS.csv"
)

scale_rows = []

for i, call in enumerate(
    scale_calls,
    1
):

    scale_rows.append({
        "call_number":
            i,

        "line_number":
            call[
                "line_number"
            ],

        "exact_call":
            call[
                "call"
            ]
    })

    print(
        f"SCALE CALL {i}"
    )

    print(
        "  line  :",
        call[
            "line_number"
        ]
    )

    print(
        "  exact :",
        call[
            "call"
        ]
    )

    print()


write_csv(
    SCALE_OUT,
    scale_rows,
    [
        "call_number",
        "line_number",
        "exact_call"
    ]
)


# ============================================================
# 8. Exact rowMeans() calls
# ============================================================

print("===== 7. EXACT rowMeans() CALLS =====")
print()

rowmean_calls = extract_balanced_calls(
    source_text,
    code_text,
    r"(?<![\w.])rowMeans\s*\("
)

print(
    "Executable rowMeans() calls detected:",
    len(rowmean_calls)
)

print()

ROWMEAN_OUT = (
    META_DIR
    / "STEP22A_05E_EXACT_ROWMEANS_CALLS.csv"
)

rowmean_rows = []

for i, call in enumerate(
    rowmean_calls,
    1
):

    rowmean_rows.append({
        "call_number":
            i,

        "line_number":
            call[
                "line_number"
            ],

        "exact_call":
            call[
                "call"
            ]
    })

    print(
        f"ROWMEANS CALL {i}"
    )

    print(
        "  line  :",
        call[
            "line_number"
        ]
    )

    print(
        "  exact :",
        call[
            "call"
        ]
    )

    print()


write_csv(
    ROWMEAN_OUT,
    rowmean_rows,
    [
        "call_number",
        "line_number",
        "exact_call"
    ]
)


# ============================================================
# 9. CORE2 / CORE4 exact source context
# ============================================================

print("===== 8. CORE2 / CORE4 EXACT SOURCE CONTEXT =====")
print()

score_patterns = [
    r"CORE2",
    r"CORE4",
    r"CORE2_STRONG",
    r"CDKN1B",
    r"GSN",
    r"BGN",
    r"TIMP1"
]

hit_lines = set()

for pattern in score_patterns:

    for i, line in enumerate(
        lines,
        1
    ):

        if re.search(
            pattern,
            line,
            flags=re.I
        ):

            hit_lines.add(i)


score_context = context_around_lines(
    lines,
    sorted(
        hit_lines
    ),
    radius=3
)

SCORE_CONTEXT_OUT = (
    META_DIR
    / "STEP22A_05E_CORE_SCORE_SOURCE_CONTEXT.csv"
)

write_csv(
    SCORE_CONTEXT_OUT,
    score_context,
    [
        "line_number",
        "code"
    ]
)


# Print only contexts that look like actual score construction.
construction_lines = []

for row in score_context:

    code = row[
        "code"
    ]

    if re.search(
        r"CORE2|CORE4|scale\s*\(|rowMeans\s*\(",
        code,
        flags=re.I
    ):

        construction_lines.append(
            row
        )


print(
    "Relevant construction/context lines:"
)

print()

for row in construction_lines[:100]:

    print(
        f"{row['line_number']:05d} | "
        f"{row['code']}"
    )

print()


# ============================================================
# 10. Implementation vocabulary
# ============================================================

print("===== 9. IMPLEMENTATION VOCABULARY =====")
print()

implementation_patterns = {
    "cpm":
        r"(?<![\w.])(?:edgeR::)?cpm\s*\(",

    "scale":
        r"(?<![\w.])scale\s*\(",

    "rowMeans":
        r"(?<![\w.])rowMeans\s*\(",

    "log2":
        r"(?<![\w.])log2\s*\(",

    "DGEList":
        r"\bDGEList\b",

    "calcNormFactors":
        r"\bcalcNormFactors\b",

    "voom":
        r"\bvoom\b",

    "TMM":
        r"\bTMM\b",

    "prior.count":
        r"\bprior\.count\b",

    "normalized.lib.sizes":
        r"\bnormalized\.lib\.sizes\b"
}

implementation_rows = []

for name, pattern in (
    implementation_patterns.items()
):

    matches = list(
        re.finditer(
            pattern,
            code_text,
            flags=re.I
        )
    )

    line_numbers = [
        str(
            line_number(
                source_text,
                m.start()
            )
        )
        for m in matches
    ]

    implementation_rows.append({
        "term":
            name,

        "detected":
            "YES"
            if matches
            else "NO",

        "hit_count":
            len(matches),

        "line_numbers":
            ";".join(
                line_numbers
            )
    })


IMPLEMENTATION_OUT = (
    META_DIR
    / "STEP22A_05E_IMPLEMENTATION_VOCABULARY.csv"
)

write_csv(
    IMPLEMENTATION_OUT,
    implementation_rows,
    [
        "term",
        "detected",
        "hit_count",
        "line_numbers"
    ]
)

for row in implementation_rows:

    print(
        f"{row['term']:<24s} "
        f"{row['detected']:<3s} "
        f"hits={row['hit_count']}"
    )

print()


# ============================================================
# 11. Frozen output HEADERS only
#
# NO data rows are read.
# ============================================================

print("===== 10. FROZEN STEP19L OUTPUT HEADERS ONLY =====")
print()

search_roots = [
    PROJECT
    / "04_results"
    / "cross_dataset"
    / "pathway_analysis"
    / "Step19L",

    PROJECT
    / "06_tables"
    / "cross_dataset"
    / "pathway_analysis"
    / "Step19L"
]

wanted_names = {
    "STEP19L_SAMPLE_CORE_MODULE_SCORES.csv",
    "STEP19L_STRONG_VS_WEAK_CORE_MODULE.csv",
    "STEP19L_CORE_MODULE_CANONICAL_EFFECTS.csv"
}

header_rows = []

for root in search_roots:

    if not root.exists():
        continue

    for path in root.rglob(
        "*.csv"
    ):

        if (
            path.name
            not in wanted_names
        ):
            continue

        # Read ONLY first line.
        with open(
            path,
            "r",
            encoding="utf-8-sig",
            errors="replace"
        ) as f:

            header_line = (
                f.readline()
                .rstrip(
                    "\r\n"
                )
            )

        header_rows.append({
            "file":
                str(
                    path.relative_to(
                        PROJECT
                    )
                ),

            "size_bytes":
                path.stat().st_size,

            "sha256":
                sha256_file(
                    path
                ),

            "header_only":
                header_line,

            "data_rows_read":
                "NO"
        })


HEADER_OUT = (
    META_DIR
    / "STEP22A_05E_FROZEN_OUTPUT_HEADERS_ONLY.csv"
)

write_csv(
    HEADER_OUT,
    header_rows,
    [
        "file",
        "size_bytes",
        "sha256",
        "header_only",
        "data_rows_read"
    ]
)

print(
    "Frozen Step19L score-related output files found:",
    len(
        header_rows
    )
)

for row in header_rows:

    print()

    print(
        "FILE:",
        row[
            "file"
        ]
    )

    print(
        "HEADER:",
        row[
            "header_only"
        ]
    )

print()

print(
    "Frozen result data rows read: NO"
)

print()


# ============================================================
# 12. Frozen specification consistency check
#
# Check TEXT only.
#
# The historical frozen specification requires:
#
# CORE4:
# +BGN
# -CDKN1B
# -GSN
# +TIMP1
#
# CORE2_STRONG:
# -CDKN1B
# -GSN
#
# This step still does NOT calculate scores.
# ============================================================

print("===== 11. FROZEN FORMULA SPECIFICATION CHECK =====")
print()

required_tokens = [
    "BGN",
    "CDKN1B",
    "GSN",
    "TIMP1",
    "CORE4",
    "CORE2"
]

token_status = {}

for token in required_tokens:

    token_status[
        token
    ] = bool(
        re.search(
            re.escape(
                token
            ),
            source_text,
            flags=re.I
        )
    )

for token in required_tokens:

    print(
        f"{token:<12s}:",
        "PRESENT"
        if token_status[token]
        else "MISSING"
    )

print()

if not all(
    token_status.values()
):

    raise SystemExit(
        "FROZEN_CORE_FORMULA_TEXT_INCOMPLETE"
    )

print(
    "✓ Frozen CORE gene/module definitions present in source."
)

print(
    "✓ No external score calculated."
)

print()


# ============================================================
# 13. Lock manifest
# ============================================================

LOCK_MANIFEST = (
    META_DIR
    / "STEP22A_05E_FROZEN_SCORING_LOCK_MANIFEST.csv"
)

manifest_rows = [
    {
        "item":
            "frozen_source_path",

        "value":
            str(
                STEP19L.relative_to(
                    PROJECT
                )
            )
    },
    {
        "item":
            "frozen_source_sha256",

        "value":
            source_sha
    },
    {
        "item":
            "source_copy_path",

        "value":
            str(
                SOURCE_COPY.relative_to(
                    PROJECT
                )
            )
    },
    {
        "item":
            "source_copy_sha256",

        "value":
            copy_sha
    },
    {
        "item":
            "source_copy_byte_identical",

        "value":
            "YES"
    },
    {
        "item":
            "section_19L3_found",

        "value":
            "YES"
    },
    {
        "item":
            "section_19L8_found",

        "value":
            "YES"
    },
    {
        "item":
            "executable_cpm_calls",

        "value":
            str(
                len(
                    cpm_calls
                )
            )
    },
    {
        "item":
            "executable_scale_calls",

        "value":
            str(
                len(
                    scale_calls
                )
            )
    },
    {
        "item":
            "executable_rowMeans_calls",

        "value":
            str(
                len(
                    rowmean_calls
                )
            )
    },
    {
        "item":
            "OMIX_normalization_executed",

        "value":
            "NO"
    },
    {
        "item":
            "OMIX_CORE_gene_values_read",

        "value":
            "NO"
    },
    {
        "item":
            "OMIX_CORE2_score_calculated",

        "value":
            "NO"
    },
    {
        "item":
            "OMIX_CORE4_score_calculated",

        "value":
            "NO"
    },
    {
        "item":
            "before_after_inference_run",

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
            "lock_status",

        "value":
            "EXACT_FROZEN_SOURCE_CAPTURED_STATIC_ONLY"
    }
]

write_csv(
    LOCK_MANIFEST,
    manifest_rows,
    [
        "item",
        "value"
    ]
)


# ============================================================
# 14. Human-readable summary
# ============================================================

SUMMARY = (
    META_DIR
    / "STEP22A_05E_EXACT_FROZEN_SCORING_SOURCE_SUMMARY.txt"
)

with open(
    SUMMARY,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-5E EXACT FROZEN SCORING SOURCE AUDIT\n"
    )

    f.write(
        "=" * 65
        + "\n\n"
    )

    f.write(
        f"Frozen Step19L source: {STEP19L}\n"
    )

    f.write(
        f"SHA256: {source_sha}\n"
    )

    f.write(
        "Byte-identical frozen source copy: YES\n"
    )

    f.write(
        f"Section 19L3: {sec3_range}\n"
    )

    f.write(
        f"Section 19L8: {sec8_range}\n"
    )

    f.write(
        f"Executable cpm calls: {len(cpm_calls)}\n"
    )

    f.write(
        f"Executable scale calls: {len(scale_calls)}\n"
    )

    f.write(
        f"Executable rowMeans calls: {len(rowmean_calls)}\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "OMIX normalization executed: NO\n"
    )

    f.write(
        "OMIX CORE gene values read: NO\n"
    )

    f.write(
        "OMIX CORE2 calculated: NO\n"
    )

    f.write(
        "OMIX CORE4 calculated: NO\n"
    )

    f.write(
        "Before/After inference: NO\n"
    )

    f.write(
        "Histology mapping: INCOMPLETE\n"
    )

    f.write(
        "Primary CORE2 inference: BLOCKED\n"
    )

    f.write(
        "\n"
    )

    f.write(
        "LOCK STATUS:\n"
    )

    f.write(
        "EXACT_FROZEN_SOURCE_CAPTURED_STATIC_ONLY\n"
    )


# ============================================================
# 15. Frozen project re-check
# ============================================================

print("===== 12. FINAL FROZEN PROJECT CHECK =====")
print()

changes_after = []

for row in frozen:

    path = Path(
        row.get(
            "path",
            ""
        )
    )

    if not path.is_absolute():

        path = PROJECT / path

    if not path.exists():

        changes_after.append(
            (
                str(path),
                "MISSING"
            )
        )

        continue

    expected_size = row.get(
        "size_bytes",
        ""
    )

    if expected_size:

        try:

            expected_size = int(
                float(
                    expected_size
                )
            )

            actual_size = (
                path.stat().st_size
            )

            if (
                expected_size
                != actual_size
            ):

                changes_after.append(
                    (
                        str(path),
                        "SIZE_CHANGED"
                    )
                )

        except Exception:
            pass


print(
    "Frozen original files changed:",
    len(
        changes_after
    )
)

if changes_after:

    for x in changes_after[:20]:
        print(x)

    raise SystemExit(
        "STEP22A_FROZEN_PROJECT_CHANGED"
    )

print(
    "✓ Frozen Step19-21 inventory unchanged."
)

print()


# ============================================================
# FINAL
# ============================================================

print("=" * 78)
print("STEP22A-5E COMPLETE")
print("=" * 78)
print()

print(
    "FROZEN STEP19L SOURCE        : LOCKED"
)

print(
    "BYTE-IDENTICAL SOURCE COPY   : YES"
)

print(
    "STEP19L SHA256               :",
    source_sha
)

print()

print(
    "SECTION 19L3 FOUND           : YES"
)

print(
    "SECTION 19L8 FOUND           : YES"
)

print()

print(
    "EXECUTABLE cpm() CALLS       :",
    len(
        cpm_calls
    )
)

print(
    "EXECUTABLE scale() CALLS     :",
    len(
        scale_calls
    )
)

print(
    "EXECUTABLE rowMeans() CALLS  :",
    len(
        rowmean_calls
    )
)

print()

print(
    "OMIX NORMALIZATION EXECUTED  : NO"
)

print(
    "OMIX CORE GENE VALUES READ   : NO"
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
    "LOCK STATUS:"
)

print(
    "EXACT_FROZEN_SOURCE_CAPTURED_STATIC_ONLY"
)

print()

print(
    "NEXT REQUIRED STEP:"
)

print(
    "STEP22A-5F = INTERPRET AND FREEZE THE EXACT"
)

print(
    "cpm / z-SCORE / CORE2 / CORE4 IMPLEMENTATION"
)

print(
    "USING ONLY THE CAPTURED STEP19L SOURCE."
)

print()

print("STOP HERE.")
print("DO NOT SCORE OMIX CORE2 YET.")
print()
