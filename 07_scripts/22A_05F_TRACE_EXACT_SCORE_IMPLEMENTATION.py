#!/usr/bin/env python3

# ============================================================
# ESCC Neoadjuvant Project
# STEP22A-5F
#
# TRACE AND LOCK EXACT FROZEN SCORING IMPLEMENTATION
#
# STATIC SOURCE AUDIT ONLY.
#
# Goals:
#   1. trace l2c_tumor source
#   2. identify exact four-gene expression aggregation
#   3. identify exact z-score implementation
#   4. lock CORE4 / CORE2_STRONG / CORE2_WEAK formulas
#
# STRICTLY FORBIDDEN:
#   - execute Step19L
#   - load OMIX CORE-gene values
#   - normalize OMIX pseudobulk
#   - calculate CORE2 / CORE4
#   - inspect Before/After CORE-gene effects
#   - modify frozen Step19-21
# ============================================================

from pathlib import Path
import csv
import hashlib
import re
import sys


PROJECT = Path.cwd()

META_DIR = (
    PROJECT
    / "00_metadata"
    / "OMIX005710"
)

SOURCE_COPY = (
    META_DIR
    / "STEP22A_05E_FROZEN_STEP19L_SOURCE_COPY.R"
)

LOCK_MANIFEST = (
    META_DIR
    / "STEP22A_05E_FROZEN_SCORING_LOCK_MANIFEST.csv"
)

ORIGINAL_STEP19L = (
    PROJECT
    / "07_scripts"
    / "19L_CONSERVED_FIBROBLAST_CORE_VALIDATION.R"
)

FROZEN_INVENTORY = (
    META_DIR
    / "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

META_DIR.mkdir(
    parents=True,
    exist_ok=True
)


print()
print("=" * 78)
print("STEP22A-5F: TRACE EXACT FROZEN SCORING IMPLEMENTATION")
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


def strip_comment(line):
    """
    Remove R comments from one line while preserving quoted '#'.
    """

    out = []
    single = False
    double = False
    escaped = False

    for c in line:

        if escaped:

            out.append(c)
            escaped = False
            continue

        if (
            c == "\\"
            and (
                single
                or double
            )
        ):

            out.append(c)
            escaped = True
            continue

        if (
            c == "'"
            and not double
        ):

            single = not single
            out.append(c)
            continue

        if (
            c == '"'
            and not single
        ):

            double = not double
            out.append(c)
            continue

        if (
            c == "#"
            and not single
            and not double
        ):

            break

        out.append(c)

    return "".join(out)


def code_only_lines(lines):

    return [
        strip_comment(x)
        for x in lines
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
        return None, None

    if end is None:
        end = len(lines)

    return start, end


def get_context(
    lines,
    hit_lines,
    radius=4
):

    use = set()

    for n in hit_lines:

        for k in range(
            max(
                1,
                n - radius
            ),
            min(
                len(lines),
                n + radius
            ) + 1
        ):

            use.add(k)

    return [
        (
            n,
            lines[n - 1]
        )
        for n in sorted(use)
    ]


def grep_lines(
    code_lines,
    pattern
):

    out = []

    for i, line in enumerate(
        code_lines,
        1
    ):

        if re.search(
            pattern,
            line,
            flags=re.I
        ):

            out.append(i)

    return out


def print_context(
    title,
    rows
):

    print(title)
    print("-" * 78)

    if not rows:

        print("<none>")
        print()
        return

    for n, line in rows:

        print(
            f"{n:05d} | {line}"
        )

    print()


# ============================================================
# 1. Verify 5E source lock
# ============================================================

print("===== 1. VERIFY STEP22A-5E SOURCE LOCK =====")
print()

for path in [
    SOURCE_COPY,
    LOCK_MANIFEST,
    ORIGINAL_STEP19L
]:

    if not path.exists():

        raise SystemExit(
            "ERROR: missing required file:\n"
            + str(path)
        )


manifest_rows = read_csv(
    LOCK_MANIFEST
)

manifest = {
    x["item"]:
        x["value"]
    for x in manifest_rows
}


source_sha = sha256_file(
    SOURCE_COPY
)

original_sha = sha256_file(
    ORIGINAL_STEP19L
)

locked_sha = manifest.get(
    "frozen_source_sha256",
    ""
)


print(
    "5E locked SHA256 :",
    locked_sha
)

print(
    "Source-copy SHA  :",
    source_sha
)

print(
    "Original SHA     :",
    original_sha
)

print()


if (
    not locked_sha
    or
    source_sha != locked_sha
    or
    original_sha != locked_sha
):

    raise SystemExit(
        "STEP22A_FROZEN_STEP19L_SHA_MISMATCH"
    )


print(
    "✓ Original Step19L == frozen 5E copy."
)

print(
    "✓ Byte-level scoring source remains frozen."
)

print()


# ============================================================
# 2. Read source
# ============================================================

raw = SOURCE_COPY.read_text(
    encoding="utf-8",
    errors="replace"
)

lines = raw.splitlines()

code_lines = code_only_lines(
    lines
)


print(
    "Frozen Step19L lines:",
    len(lines)
)

print()


# ============================================================
# 3. Locate exact 19L3 section
# ============================================================

print("===== 2. LOCATE SECTION 19L3 =====")
print()

sec_start, sec_end = find_section(
    lines,
    r"\b19L3\b",
    r"\b19L4\b"
)


if (
    sec_start is None
    or sec_end is None
):

    raise SystemExit(
        "STEP19L_SECTION_19L3_NOT_FOUND"
    )


print(
    "19L3 source range:",
    f"{sec_start}-{sec_end}"
)

print()


SECTION_OUT = (
    META_DIR
    / "STEP22A_05F_EXACT_19L3_SOURCE.txt"
)

with open(
    SECTION_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-5F EXACT STEP19L SECTION 19L3\n"
    )

    f.write(
        f"SHA256: {source_sha}\n"
    )

    f.write(
        f"Lines: {sec_start}-{sec_end}\n\n"
    )

    for n in range(
        sec_start,
        sec_end + 1
    ):

        f.write(
            f"{n:05d} | "
            f"{lines[n - 1]}\n"
        )


# ============================================================
# 4. Trace l2c_tumor
# ============================================================

print("===== 3. TRACE l2c_tumor =====")
print()

l2c_tumor_hits = grep_lines(
    code_lines,
    r"\bl2c_tumor\b"
)

l2c_hits = grep_lines(
    code_lines,
    r"\bl2c\b"
)


l2c_tumor_context = get_context(
    lines,
    l2c_tumor_hits,
    radius=6
)

l2c_context = get_context(
    lines,
    l2c_hits,
    radius=5
)


print_context(
    "l2c_tumor source context:",
    l2c_tumor_context
)


L2C_CONTEXT_OUT = (
    META_DIR
    / "STEP22A_05F_L2C_SOURCE_CONTEXT.csv"
)

l2c_rows = []

for n, text in sorted(
    set(
        l2c_tumor_context
        + l2c_context
    )
):

    l2c_rows.append({
        "line_number":
            n,

        "code":
            text
    })


write_csv(
    L2C_CONTEXT_OUT,
    l2c_rows,
    [
        "line_number",
        "code"
    ]
)


# ============================================================
# 5. Look specifically for executable normalization
#
# Comments already removed.
# ============================================================

print("===== 4. EXECUTABLE NORMALIZATION AUDIT =====")
print()

normalization_patterns = {

    "cpm":
        r"(?<![\w.])(?:edgeR::)?cpm\s*\(",

    "DGEList":
        r"\bDGEList\s*\(",

    "calcNormFactors":
        r"\bcalcNormFactors\s*\(",

    "voom":
        r"\bvoom\s*\(",

    "TMM":
        r"\bTMM\b",

    "logCPM":
        r"\blogCPM\b",

    "log2":
        r"(?<![\w.])log2\s*\(",

    "log1p":
        r"(?<![\w.])log1p\s*\("
}


normalization_rows = []

for name, pattern in (
    normalization_patterns.items()
):

    hits = grep_lines(
        code_lines,
        pattern
    )

    normalization_rows.append({
        "operation":
            name,

        "executable_hits":
            len(hits),

        "line_numbers":
            ";".join(
                str(x)
                for x in hits
            )
    })

    print(
        f"{name:<20s}: "
        f"{len(hits)}"
    )


print()


NORMALIZATION_OUT = (
    META_DIR
    / "STEP22A_05F_EXECUTABLE_NORMALIZATION_AUDIT.csv"
)

write_csv(
    NORMALIZATION_OUT,
    normalization_rows,
    [
        "operation",
        "executable_hits",
        "line_numbers"
    ]
)


# ============================================================
# 6. Exact four-gene extraction source
# ============================================================

print("===== 5. FOUR-GENE EXPRESSION CONSTRUCTION =====")
print()

gene_names = [
    "BGN",
    "CDKN1B",
    "GSN",
    "TIMP1"
]


gene_rows = []

for gene in gene_names:

    hits = grep_lines(
        code_lines,
        rf"\b{gene}\b"
    )

    # Restrict primary display to 19L3.
    hits_19l3 = [
        x
        for x in hits
        if (
            sec_start
            <= x
            <= sec_end
        )
    ]

    context = get_context(
        lines,
        hits_19l3,
        radius=3
    )

    print_context(
        f"{gene} context inside 19L3:",
        context
    )

    for n, text in context:

        gene_rows.append({
            "gene":
                gene,

            "line_number":
                n,

            "code":
                text
        })


GENE_SOURCE_OUT = (
    META_DIR
    / "STEP22A_05F_FOUR_GENE_SOURCE_CONTEXT.csv"
)

write_csv(
    GENE_SOURCE_OUT,
    gene_rows,
    [
        "gene",
        "line_number",
        "code"
    ]
)


# ============================================================
# 7. Trace rowMeans
# ============================================================

print("===== 6. rowMeans() IMPLEMENTATION =====")
print()

rowmeans_hits = grep_lines(
    code_lines,
    r"(?<![\w.])rowMeans\s*\("
)

rowmeans_19l3 = [
    x
    for x in rowmeans_hits
    if (
        sec_start
        <= x
        <= sec_end
    )
]


rowmeans_context = get_context(
    lines,
    rowmeans_19l3,
    radius=3
)


print_context(
    "rowMeans() context inside 19L3:",
    rowmeans_context
)


ROWMEANS_OUT = (
    META_DIR
    / "STEP22A_05F_ROWMEANS_SOURCE_CONTEXT.csv"
)

write_csv(
    ROWMEANS_OUT,
    [
        {
            "line_number": n,
            "code": text
        }
        for n, text
        in rowmeans_context
    ],
    [
        "line_number",
        "code"
    ]
)


# ============================================================
# 8. Trace z-score implementation
# ============================================================

print("===== 7. EXACT Z-SCORE IMPLEMENTATION =====")
print()

z_patterns = [
    r"\bBGN_z\b",
    r"\bCDKN1B_z\b",
    r"\bGSN_z\b",
    r"\bTIMP1_z\b",
    r"\bzscore\b",
    r"\bz_score\b",
    r"\bz\.score\b",
    r"(?<![\w.])scale\s*\(",
    r"\bmean\s*\(",
    r"\bsd\s*\("
]


z_hits = set()

for pattern in z_patterns:

    for n in grep_lines(
        code_lines,
        pattern
    ):

        # Include 19L3 plus helper functions before 19L3.
        if n <= sec_end:

            z_hits.add(n)


z_context = get_context(
    lines,
    sorted(
        z_hits
    ),
    radius=4
)


print_context(
    "Z-score source context:",
    z_context
)


ZSCORE_OUT = (
    META_DIR
    / "STEP22A_05F_EXACT_ZSCORE_SOURCE_CONTEXT.csv"
)

write_csv(
    ZSCORE_OUT,
    [
        {
            "line_number": n,
            "code": text
        }
        for n, text
        in z_context
    ],
    [
        "line_number",
        "code"
    ]
)


# ============================================================
# 9. Check whether scale() actually belongs to gene scoring
# ============================================================

print("===== 8. scale() ROLE AUDIT =====")
print()

scale_hits = grep_lines(
    code_lines,
    r"(?<![\w.])scale\s*\("
)


scale_role_rows = []

for n in scale_hits:

    context = " ".join(
        lines[
            max(
                0,
                n - 3
            ):
            min(
                len(lines),
                n + 2
            )
        ]
    )

    gene_related = bool(
        re.search(
            r"BGN|CDKN1B|GSN|TIMP1|CORE2|CORE4",
            context,
            flags=re.I
        )
    )

    covariate_related = bool(
        re.search(
            r"fibro|abundance|fraction|lib|library",
            context,
            flags=re.I
        )
    )

    if gene_related:

        role = (
            "POSSIBLE_GENE_SCORING_SCALE"
        )

    elif covariate_related:

        role = (
            "COVARIATE_SCALE"
        )

    else:

        role = (
            "OTHER_OR_REVIEW"
        )

    scale_role_rows.append({
        "line_number":
            n,

        "role":
            role,

        "code":
            lines[n - 1]
    })


for row in scale_role_rows:

    print(
        f"{row['line_number']:05d} | "
        f"{row['role']:<30s} | "
        f"{row['code']}"
    )


print()


SCALE_ROLE_OUT = (
    META_DIR
    / "STEP22A_05F_SCALE_ROLE_AUDIT.csv"
)

write_csv(
    SCALE_ROLE_OUT,
    scale_role_rows,
    [
        "line_number",
        "role",
        "code"
    ]
)


# ============================================================
# 10. Exact frozen CORE formula lines
# ============================================================

print("===== 9. EXACT CORE FORMULAS =====")
print()

formula_patterns = {

    "CORE4":
        r"\bCORE4_score\b",

    "CORE2_STRONG":
        r"\bCORE2_STRONG_score\b",

    "CORE2_WEAK":
        r"\bCORE2_WEAK_score\b"
}


formula_rows = []

for name, pattern in (
    formula_patterns.items()
):

    hits = grep_lines(
        code_lines,
        pattern
    )

    # Prefer assignment-like lines.
    assignment_hits = [
        n
        for n in hits
        if re.search(
            r"<-|:=|=",
            code_lines[
                n - 1
            ]
        )
    ]

    selected_hits = (
        assignment_hits
        if assignment_hits
        else hits
    )

    print(
        name,
        ":"
    )

    if not selected_hits:

        print(
            "  <not found>"
        )

    for n in selected_hits:

        print(
            f"  {n:05d} | "
            f"{lines[n - 1]}"
        )

        formula_rows.append({
            "score":
                name,

            "line_number":
                n,

            "code":
                lines[n - 1]
        })

    print()


FORMULA_OUT = (
    META_DIR
    / "STEP22A_05F_EXACT_CORE_FORMULAS.csv"
)

write_csv(
    FORMULA_OUT,
    formula_rows,
    [
        "score",
        "line_number",
        "code"
    ]
)


# ============================================================
# 11. Mechanical evidence checks
# ============================================================

print("===== 10. IMPLEMENTATION RESOLUTION CHECK =====")
print()

checks = []


def add_check(name, passed, detail):

    checks.append({
        "check":
            name,

        "status":
            (
                "PASS"
                if passed
                else "REVIEW"
            ),

        "detail":
            detail
    })


add_check(
    "19L3_found",
    sec_start is not None,
    f"{sec_start}-{sec_end}"
)


add_check(
    "l2c_tumor_referenced",
    len(
        l2c_tumor_hits
    ) > 0,
    str(
        l2c_tumor_hits
    )
)


add_check(
    "rowMeans_in_19L3",
    len(
        rowmeans_19l3
    ) > 0,
    str(
        rowmeans_19l3
    )
)


for gene in gene_names:

    hits = grep_lines(
        code_lines,
        rf"\b{gene}_z\b"
    )

    add_check(
        f"{gene}_z_present",
        len(hits) > 0,
        str(hits)
    )


for score in [
    "CORE4",
    "CORE2_STRONG"
]:

    hits = [
        x
        for x in formula_rows
        if x["score"] == score
    ]

    add_check(
        f"{score}_formula_present",
        len(hits) > 0,
        ";".join(
            str(
                x["line_number"]
            )
            for x in hits
        )
    )


exec_cpm = next(
    (
        x["executable_hits"]
        for x in normalization_rows
        if x["operation"]
        == "cpm"
    ),
    0
)


add_check(
    "no_executable_cpm",
    exec_cpm == 0,
    str(
        exec_cpm
    )
)


exec_tmm = sum(
    x["executable_hits"]
    for x in normalization_rows
    if x["operation"]
    in {
        "DGEList",
        "calcNormFactors",
        "voom",
        "TMM"
    }
)


add_check(
    "no_executable_TMM_voom_path",
    exec_tmm == 0,
    str(
        exec_tmm
    )
)


CHECK_OUT = (
    META_DIR
    / "STEP22A_05F_IMPLEMENTATION_RESOLUTION_CHECK.csv"
)

write_csv(
    CHECK_OUT,
    checks,
    [
        "check",
        "status",
        "detail"
    ]
)


for row in checks:

    print(
        f"{row['check']:<34s} "
        f"{row['status']:<7s} "
        f"{row['detail']}"
    )


print()


# ============================================================
# 12. Produce exact human-readable lock file
#
# IMPORTANT:
# No attempt is made to invent missing semantics.
# Exact source context is the source of truth.
# ============================================================

LOCK_OUT = (
    META_DIR
    / "STEP22A_05F_EXACT_SCORING_IMPLEMENTATION_LOCK.txt"
)

with open(
    LOCK_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "STEP22A-5F EXACT FROZEN SCORING IMPLEMENTATION LOCK\n"
    )

    f.write(
        "=" * 70
        + "\n\n"
    )

    f.write(
        f"Frozen Step19L SHA256:\n{source_sha}\n\n"
    )

    f.write(
        f"19L3 source range: {sec_start}-{sec_end}\n\n"
    )

    f.write(
        "============================================================\n"
    )

    f.write(
        "L2C / l2c_tumor SOURCE CONTEXT\n"
    )

    f.write(
        "============================================================\n"
    )

    for row in l2c_rows:

        f.write(
            f"{row['line_number']:05d} | {row['code']}\n"
        )

    f.write(
        "\n============================================================\n"
    )

    f.write(
        "rowMeans SOURCE CONTEXT\n"
    )

    f.write(
        "============================================================\n"
    )

    for n, text in rowmeans_context:

        f.write(
            f"{n:05d} | {text}\n"
        )

    f.write(
        "\n============================================================\n"
    )

    f.write(
        "Z-SCORE SOURCE CONTEXT\n"
    )

    f.write(
        "============================================================\n"
    )

    for n, text in z_context:

        f.write(
            f"{n:05d} | {text}\n"
        )

    f.write(
        "\n============================================================\n"
    )

    f.write(
        "CORE FORMULAS\n"
    )

    f.write(
        "============================================================\n"
    )

    for row in formula_rows:

        f.write(
            f"{row['score']} | "
            f"{row['line_number']:05d} | "
            f"{row['code']}\n"
        )

    f.write(
        "\n============================================================\n"
    )

    f.write(
        "STATIC AUDIT STATUS\n"
    )

    f.write(
        "============================================================\n"
    )

    f.write(
        "OMIX expression values read: NO\n"
    )

    f.write(
        "OMIX normalization executed: NO\n"
    )

    f.write(
        "CORE2 calculated: NO\n"
    )

    f.write(
        "CORE4 calculated: NO\n"
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


# ============================================================
# 13. Final frozen-file integrity check
# ============================================================

print("===== 11. FINAL FROZEN PROJECT CHECK =====")
print()

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

    expected = row.get(
        "size_bytes",
        ""
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

            if actual != expected:

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
# FINAL PRINT
# ============================================================

review_checks = [
    x
    for x in checks
    if x["status"]
    != "PASS"
]


print("=" * 78)
print("STEP22A-5F COMPLETE")
print("=" * 78)
print()

print(
    "FROZEN STEP19L SHA256        :",
    source_sha
)

print(
    "SECTION 19L3                :",
    f"{sec_start}-{sec_end}"
)

print()

print(
    "EXECUTABLE cpm()             :",
    exec_cpm
)

print(
    "EXECUTABLE TMM/voom PATH     :",
    exec_tmm
)

print(
    "rowMeans() IN 19L3           :",
    len(
        rowmeans_19l3
    )
)

print()

print(
    "IMPLEMENTATION CHECKS PASS   :",
    len(checks)
    -
    len(review_checks),
    "/",
    len(checks)
)

print(
    "IMPLEMENTATION CHECKS REVIEW :",
    len(
        review_checks
    )
)

print()

if review_checks:

    print(
        "REVIEW ITEMS:"
    )

    for row in review_checks:

        print(
            " ",
            row["check"],
            "->",
            row["detail"]
        )

    print()


print(
    "OMIX CORE GENE VALUES READ   : NO"
)

print(
    "OMIX NORMALIZATION EXECUTED  : NO"
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
    "EXACT LOCK FILE:"
)

print(
    "  00_metadata/OMIX005710/"
)

print(
    "  STEP22A_05F_EXACT_SCORING_IMPLEMENTATION_LOCK.txt"
)

print()

if review_checks:

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "RESOLVE ONLY THE SOURCE-LINE REVIEW ITEMS ABOVE."
    )

    print(
        "DO NOT SCORE OMIX CORE2."
    )

else:

    print(
        "SCORING IMPLEMENTATION SOURCE LOCK : COMPLETE"
    )

    print()

    print(
        "NEXT REQUIRED STEP:"
    )

    print(
        "STEP22A-5G = FREEZE HISTOLOGY-UNRESOLVED"
    )

    print(
        "SENSITIVITY DESIGN BEFORE ANY CORE2 CALCULATION."
    )

print()

print("STOP HERE.")
print("DO NOT SCORE OMIX CORE2 YET.")
print()
