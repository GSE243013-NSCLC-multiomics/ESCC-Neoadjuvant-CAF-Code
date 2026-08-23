#!/bin/bash

set -euo pipefail

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"
META_DIR="$PROJECT/00_metadata/OMIX005710"
LOG_DIR="$PROJECT/08_logs"

cd "$PROJECT"

echo ""
echo "============================================================"
echo "STEP22A-4C: FREEZE HISTOLOGY GATE"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# 1. Required previous outputs
# ------------------------------------------------------------

echo "===== 1. VERIFY STEP22A-4A / 4B ====="
echo ""

REQ1="$META_DIR/STEP22A_04_AUTHOR_SOURCE_HISTOLOGY_SUMMARY.txt"
REQ2="$META_DIR/STEP22A_04B_HISTOLOGY_MAPPING_DECISION.txt"

for f in "$REQ1" "$REQ2"
do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing:"
        echo "$f"
        exit 1
    fi

    echo "✓ $f"
done

echo ""

# ------------------------------------------------------------
# 2. Record scientific decision
# ------------------------------------------------------------

cat > "$META_DIR/STEP22A_HISTOLOGY_GATE.txt" <<'TXT'
STEP22A HISTOLOGY GATE
============================================================

SOURCE PAPER FACT:
The source publication reports one EASC patient among the
22-patient cohort.

PUBLIC PATIENT-LEVEL MAPPING:
UNRESOLVED

PUBLIC SUPPLEMENT AUDIT:
COMPLETED

PUBLIC HRI AUDIT:
COMPLETED

PUBLIC HISTOLOGY SOURCES EXHAUSTED:
YES

EASC PATIENT:
UNKNOWN

IMPORTANT ANALYSIS RULES
------------------------------------------------------------

ALLOWED NOW:

1. Technical matrix preprocessing.
2. Per-sample QC audit.
3. Broad major-cell-type annotation.
4. Identification of a broad Fibroblast_CAF compartment.
5. Fibroblast cell-count QC.
6. Creation of sample/patient-level fibroblast pseudobulk
   objects for technical preparation only.

NOT ALLOWED YET:

1. Primary CORE2 before/after inferential validation.
2. Wilcoxon paired CORE2 test.
3. Exact sign-flip CORE2 test.
4. CORE2 effect-size interpretation.
5. CORE4 external comparison.
6. Three-cohort replication claim.
7. Manuscript modification.
8. Claim that all analyzed patients are ESCC.

CELL ANNOTATION RULE
------------------------------------------------------------

Because no usable public author cell-level annotation was found,
use only the frozen constrained major-cell-type annotation.

Fibroblast-positive markers may include:

COL1A1
COL1A2
DCN
LUM
COL3A1
COL6A1
COL6A2
COL6A3
DPT
CFD
PDGFRA
FAP

Exclude epithelial-dominant cells using:

EPCAM
KRT8
KRT18
KRT19
KRT5
KRT14

Exclude immune cells using:

PTPRC
CD3D
CD3E
NKG7
LST1
MS4A1
CD79A

Exclude endothelial cells using:

PECAM1
VWF
EMCN
KDR

Distinguish mural/pericyte-like cells using:

RGS5
CSPG4
MCAM
NOTCH3

FORBIDDEN FOR CELL SELECTION:

CDKN1B
GSN
BGN
TIMP1
CORE2
CORE4

Reason:
Using CORE genes to define fibroblasts would create circular
external validation.

CURRENT HISTOLOGY STATUS:

HISTOLOGY_MAPPING_INCOMPLETE

CURRENT ANALYSIS STATUS:

TECHNICAL_PREPROCESSING_ALLOWED
PRIMARY_CORE2_INFERENCE_BLOCKED
TXT

echo "✓ Histology gate written."
echo ""

# ------------------------------------------------------------
# 3. Machine-readable status
# ------------------------------------------------------------

python3 - <<'PY'
import csv
from datetime import datetime

out = (
    "00_metadata/OMIX005710/"
    "STEP22A_ANALYSIS_GATE_STATUS.csv"
)

rows = [
    {
        "gate": "histology_mapping",
        "status": "INCOMPLETE",
        "allowed": "NO",
        "reason": "EASC patient cannot be mapped from public data"
    },
    {
        "gate": "technical_preprocessing",
        "status": "OPEN",
        "allowed": "YES",
        "reason": "Does not require histology-specific inference"
    },
    {
        "gate": "broad_cell_annotation",
        "status": "OPEN",
        "allowed": "YES",
        "reason": "Frozen constrained annotation is permitted"
    },
    {
        "gate": "fibroblast_QC",
        "status": "OPEN",
        "allowed": "YES",
        "reason": "Technical sample-level preparation only"
    },
    {
        "gate": "CORE2_scoring",
        "status": "BLOCKED_FOR_INFERENCE",
        "allowed": "NO",
        "reason": "Histology mapping incomplete"
    },
    {
        "gate": "paired_CORE2_test",
        "status": "BLOCKED",
        "allowed": "NO",
        "reason": "Histology mapping incomplete"
    },
    {
        "gate": "manuscript_update",
        "status": "BLOCKED",
        "allowed": "NO",
        "reason": "External validation incomplete"
    },
]

for row in rows:
    row["recorded_at"] = datetime.now().isoformat(
        timespec="seconds"
    )

with open(
    out,
    "w",
    newline="",
    encoding="utf-8-sig"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "gate",
            "status",
            "allowed",
            "reason",
            "recorded_at"
        ]
    )

    writer.writeheader()
    writer.writerows(rows)

print("Saved:")
print(out)
PY

echo ""

# ------------------------------------------------------------
# 4. Verify exactly 12 paired TARs still present
# ------------------------------------------------------------

echo "===== 2. DATA SAFETY CHECK ====="
echo ""

TAR_COUNT=$(find \
  01_raw/OMIX005710/paired_tumor \
  -maxdepth 1 \
  -type f \
  -name '*.tar' \
  | wc -l \
  | tr -d ' ')

echo "Paired TAR files: $TAR_COUNT"

if [ "$TAR_COUNT" -ne 12 ]; then
    echo "ERROR: expected 12."
    exit 1
fi

echo "✓ No additional tumor files required."
echo ""

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

echo "============================================================"
echo "STEP22A-4C COMPLETE"
echo "============================================================"
echo ""
echo "HISTOLOGY MAPPING            : INCOMPLETE"
echo "EASC PATIENT                 : UNRESOLVED"
echo ""
echo "TECHNICAL PREPROCESSING      : ALLOWED"
echo "BROAD CELL ANNOTATION        : ALLOWED"
echo "FIBROBLAST QC                : ALLOWED"
echo ""
echo "CORE2 PRIMARY INFERENCE      : BLOCKED"
echo "CORE4 EXTERNAL COMPARISON    : BLOCKED"
echo "MANUSCRIPT UPDATE            : BLOCKED"
echo ""
echo "CORE GENES USED FOR ANNOTATION: NO"
echo "SAMPLES EXCLUDED             : 0"
echo "NEW OMIX TAR DOWNLOADED      : 0"
echo ""
echo "READY FOR STEP22A-5:"
echo "CONSTRAINED MAJOR-CELL-TYPE ANNOTATION"
echo ""
