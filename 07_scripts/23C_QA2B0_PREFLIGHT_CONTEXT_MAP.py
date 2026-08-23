#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import shutil
import re
from difflib import SequenceMatcher
from docx import Document

ROOT = Path.cwd()

META = ROOT / "00_metadata" / "Step23_manuscript_update"
RESULT = ROOT / "04_results" / "manuscript_update" / "Step23"

QA1B = RESULT / "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"

QA2E_LOCK = META / "STEP23C_QA2E_DECISION_LOCK.csv"
QA2G_LOCK = META / "STEP23C_QA2G_REMAINING_DECISION_LOCK.csv"

QA2B_WORKING = (
    RESULT /
    "ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA2B_WORKING.docx"
)

OUT_CSV = META / "STEP23C_QA2B0_CONTEXT_MAP.csv"
OUT_TXT = RESULT / "STEP23C_QA2B0_PREFLIGHT_REPORT.txt"
OUT_SHA = META / "STEP23C_QA2B0_PREFLIGHT_CHECKSUMS.sha256"

EXPECTED_QA1B_SHA = (
    "0e91d04f82c3a3353bc6214b94312f72767f7da50778d40d0f33c016bc65f598"
)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def norm(s):
    s = s or ""
    s = s.replace("\u2013", "-").replace("\u2014", "-")
    s = s.replace("\u2212", "-")
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s


def first_existing(row, names):
    for n in names:
        if n in row and str(row[n]).strip():
            return str(row[n]).strip()
    return ""


print()
print("=" * 92)
print("STEP23C-QA2B0: CONTEXT-MAPPING PREFLIGHT")
print("=" * 92)
print()


# ============================================================
# 1. Verify frozen source
# ============================================================

for p in [QA1B, QA2E_LOCK, QA2G_LOCK]:
    if not p.exists():
        raise SystemExit(f"STOP: missing required input: {p}")

qa1b_sha = sha256(QA1B)

print("QA1B SHA256:")
print(" ", qa1b_sha)
print()

if qa1b_sha != EXPECTED_QA1B_SHA:
    raise SystemExit(
        "STOP: QA1B SHA256 no longer matches frozen manuscript."
    )

print("✓ Frozen QA1B manuscript verified.")
print()


# ============================================================
# 2. Create clean QA2B working copy
# ============================================================

if QA2B_WORKING.exists():
    QA2B_WORKING.unlink()

shutil.copy2(QA1B, QA2B_WORKING)

if sha256(QA2B_WORKING) != qa1b_sha:
    raise SystemExit(
        "STOP: QA2B working copy differs immediately after copy."
    )

print("✓ New QA2B working copy created:")
print(" ", QA2B_WORKING.relative_to(ROOT))
print()


# ============================================================
# 3. Read both human decision locks
# ============================================================

e_rows = read_csv(QA2E_LOCK)
g_rows = read_csv(QA2G_LOCK)

if len(e_rows) != 15:
    raise SystemExit(
        f"STOP: expected 15 QA2E rows, found {len(e_rows)}"
    )

if len(g_rows) != 39:
    raise SystemExit(
        f"STOP: expected 39 QA2G rows, found {len(g_rows)}"
    )

print("QA2E locked rows :", len(e_rows))
print("QA2G locked rows :", len(g_rows))
print()


def canonicalize(row, source):
    return {
        "source_lock": source,

        "item_id": first_existing(
            row,
            ["item_id", "qa2c_item", "id"]
        ),

        "section": first_existing(
            row,
            ["section", "manuscript_section"]
        ),

        "paragraph": first_existing(
            row,
            ["paragraph", "paragraph_number", "paragraph_id"]
        ),

        "original_text": first_existing(
            row,
            [
                "original_text",
                "text",
                "source_text",
                "original"
            ]
        ),

        "decision": first_existing(
            row,
            [
                "human_decision",
                "decision",
                "locked_decision",
                "disposition"
            ]
        ),

        "approved": first_existing(
            row,
            [
                "approved_text_or_instruction",
                "approved_text",
                "human_revision",
                "locked_revision",
                "revision"
            ]
        ),

        "edit_group": first_existing(
            row,
            ["edit_group", "group"]
        ),
    }


all_rows = (
    [canonicalize(r, "QA2E") for r in e_rows] +
    [canonicalize(r, "QA2G") for r in g_rows]
)

ids = [r["item_id"] for r in all_rows]

if len(set(ids)) != 54:
    dup = sorted(
        x for x in set(ids)
        if ids.count(x) > 1
    )
    raise SystemExit(
        f"STOP: 54 lock rows are not unique. Duplicates: {dup}"
    )

print("✓ Combined human decisions: 54 unique items.")
print()


# ============================================================
# 4. Read actual DOCX paragraphs
# ============================================================

doc = Document(QA2B_WORKING)

paras = []
for i, p in enumerate(doc.paragraphs):
    text = p.text or ""
    paras.append({
        "index0": i,
        "index1": i + 1,
        "text": text,
        "norm": norm(text),
    })

print("DOCX paragraphs:", len(paras))
print()


# ============================================================
# 5. Map decision rows to actual paragraphs
# ============================================================

mapped = []

for r in all_rows:
    item_id = r["item_id"]
    original = r["original_text"]
    original_n = norm(original)

    candidates = []

    # --------------------------------------------------------
    # A. Exact/substring matching across all paragraphs
    # --------------------------------------------------------

    if original_n:
        for p in paras:
            pn = p["norm"]

            if not pn:
                continue

            if original_n == pn:
                candidates.append(
                    (1.0000, p, "EXACT")
                )

            elif original_n in pn:
                candidates.append(
                    (0.9950, p, "ORIGINAL_IN_PARAGRAPH")
                )

            elif pn in original_n and len(pn) >= 20:
                candidates.append(
                    (0.9900, p, "PARAGRAPH_IN_ORIGINAL")
                )

    # --------------------------------------------------------
    # B. Use frozen paragraph number as a local hint
    #    check ±3 because previous audit numbering may differ
    # --------------------------------------------------------

    para_hint = None

    try:
        para_hint = int(r["paragraph"])
    except Exception:
        pass

    if original_n and para_hint is not None:
        for idx1 in range(
            max(1, para_hint - 3),
            min(len(paras), para_hint + 3) + 1
        ):
            p = paras[idx1 - 1]

            score = SequenceMatcher(
                None,
                original_n,
                p["norm"]
            ).ratio()

            if score >= 0.35:
                candidates.append(
                    (score, p, "PARAGRAPH_HINT_FUZZY")
                )

    # --------------------------------------------------------
    # C. Global fuzzy fallback
    # --------------------------------------------------------

    if original_n and not candidates:
        scored = []

        for p in paras:
            if not p["norm"]:
                continue

            score = SequenceMatcher(
                None,
                original_n,
                p["norm"]
            ).ratio()

            scored.append((score, p))

        scored.sort(
            key=lambda x: x[0],
            reverse=True
        )

        for score, p in scored[:3]:
            if score >= 0.30:
                candidates.append(
                    (score, p, "GLOBAL_FUZZY")
                )

    # Deduplicate candidates by paragraph
    best_by_para = {}

    for score, p, method in candidates:
        key = p["index0"]

        if (
            key not in best_by_para or
            score > best_by_para[key][0]
        ):
            best_by_para[key] = (
                score,
                p,
                method
            )

    candidates = sorted(
        best_by_para.values(),
        key=lambda x: x[0],
        reverse=True
    )

    # Packaging items legitimately have no original paragraph.
    package_like = (
        not original_n and
        (
            "PACKAG" in r["decision"].upper()
            or
            not r["paragraph"]
        )
    )

    if package_like:
        status = "PACKAGE_LEVEL"
        best_score = ""
        best_index = ""
        best_text = ""
        method = ""

    elif not candidates:
        status = "NO_MATCH"
        best_score = ""
        best_index = ""
        best_text = ""
        method = ""

    else:
        best_score, best_p, method = candidates[0]

        best_index = best_p["index1"]
        best_text = best_p["text"]

        second_score = (
            candidates[1][0]
            if len(candidates) > 1
            else None
        )

        if best_score >= 0.98:
            status = "PASS_EXACT_OR_SUBSTRING"

        elif (
            best_score >= 0.70 and
            (
                second_score is None or
                best_score - second_score >= 0.08
            )
        ):
            status = "PASS_UNIQUE_FUZZY"

        else:
            status = "REVIEW_REQUIRED"

    mapped.append({
        **r,
        "map_status": status,
        "match_method": method,
        "match_score": (
            f"{best_score:.4f}"
            if isinstance(best_score, float)
            else ""
        ),
        "docx_paragraph_index1": best_index,
        "docx_paragraph_text": best_text,
    })


# ============================================================
# 6. Write mapping CSV
# ============================================================

fields = [
    "source_lock",
    "item_id",
    "section",
    "paragraph",
    "decision",
    "edit_group",
    "original_text",
    "approved",
    "map_status",
    "match_method",
    "match_score",
    "docx_paragraph_index1",
    "docx_paragraph_text",
]

with open(
    OUT_CSV,
    "w",
    encoding="utf-8-sig",
    newline=""
) as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        extrasaction="ignore"
    )
    w.writeheader()
    w.writerows(mapped)


# ============================================================
# 7. Summarize
# ============================================================

from collections import Counter

counts = Counter(
    r["map_status"]
    for r in mapped
)

review_rows = [
    r for r in mapped
    if r["map_status"] in {
        "NO_MATCH",
        "REVIEW_REQUIRED",
    }
]

package_rows = [
    r for r in mapped
    if r["map_status"] == "PACKAGE_LEVEL"
]

report = []

report.append(
    "STEP23C-QA2B0 CONTEXT-MAPPING PREFLIGHT REPORT"
)
report.append("=" * 92)
report.append("")
report.append(f"QA1B SHA256: {qa1b_sha}")
report.append("")
report.append("QA1B modified: NO")
report.append("QA2B working copy edited: NO")
report.append("Step22A recomputed: NO")
report.append("")
report.append("TOTAL HUMAN DECISIONS: 54")
report.append("")

for k in sorted(counts):
    report.append(
        f"{k:<30s}: {counts[k]}"
    )

report.append("")
report.append(
    f"PACKAGE-LEVEL ITEMS: {len(package_rows)}"
)
report.append(
    f"UNRESOLVED/REVIEW ITEMS: {len(review_rows)}"
)
report.append("")

if review_rows:
    report.append(
        "ITEMS REQUIRING MANUAL CONTEXT CHECK:"
    )

    for r in review_rows:
        report.append(
            f"  {r['item_id']} | "
            f"{r['section']} | "
            f"frozen paragraph={r['paragraph']} | "
            f"status={r['map_status']} | "
            f"score={r['match_score']}"
        )

OUT_TXT.write_text(
    "\n".join(report),
    encoding="utf-8"
)


# ============================================================
# 8. Ensure QA2B working copy still byte-identical
# ============================================================

qa2b_sha = sha256(QA2B_WORKING)

if qa2b_sha != qa1b_sha:
    raise SystemExit(
        "CRITICAL STOP: QA2B working copy changed during preflight."
    )


# ============================================================
# 9. Checksums
# ============================================================

with open(
    OUT_SHA,
    "w",
    encoding="utf-8"
) as f:
    for p in [OUT_CSV, OUT_TXT]:
        f.write(
            f"{sha256(p)}  {p.relative_to(ROOT)}\n"
        )


# ============================================================
# 10. Console summary
# ============================================================

print("=" * 92)
print("STEP23C-QA2B0 PREFLIGHT COMPLETE")
print("=" * 92)
print()

print("QA1B DOCX MODIFIED        : NO")
print("QA2B WORKING COPY EDITED  : NO")
print("STEP22A RECOMPUTED        : NO")
print()

print("TOTAL HUMAN DECISIONS     : 54")

for k in sorted(counts):
    print(f"{k:<28s}: {counts[k]}")

print()
print(
    "UNRESOLVED/REVIEW ITEMS   :",
    len(review_rows)
)

print()
print("OUTPUTS:")
print(" ", OUT_CSV.relative_to(ROOT))
print(" ", OUT_TXT.relative_to(ROOT))
print(" ", OUT_SHA.relative_to(ROOT))
print()
print("WORKING COPY:")
print(" ", QA2B_WORKING.relative_to(ROOT))
print()

if review_rows:
    print("STATUS:")
    print(
        "STEP23C_QA2B0_CONTEXT_MAP_NEEDS_REVIEW"
    )
    print()
    print("REVIEW ITEMS:")
    for r in review_rows:
        print(
            f"  {r['item_id']}: "
            f"{r['map_status']} "
            f"(score={r['match_score']})"
        )
else:
    print("STATUS:")
    print(
        "STEP23C_QA2B0_CONTEXT_MAP_PASS"
    )

print()
print("NEXT:")
print(
    "APPLY LOCKED QA2E + QA2G DECISIONS "
    "TO QA2B WORKING COPY ONLY."
)
print()

