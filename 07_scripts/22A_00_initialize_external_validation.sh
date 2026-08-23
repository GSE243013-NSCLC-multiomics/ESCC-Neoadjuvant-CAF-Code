#!/bin/bash

set -euo pipefail

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"

cd "$PROJECT"

echo ""
echo "============================================================"
echo "STEP22A-0: INITIALIZE THIRD-COHORT EXTERNAL VALIDATION"
echo "============================================================"
echo ""

# ============================================================
# 1. Confirm project root
# ============================================================

echo "===== 1. CONFIRM PROJECT ROOT ====="
echo ""

echo "Current project:"
pwd
echo ""

REQUIRED_DIRS=(
  "04_results"
  "05_figures"
  "06_tables"
  "07_scripts"
  "08_logs"
)

for d in "${REQUIRED_DIRS[@]}"
do
    if [ ! -d "$d" ]; then
        echo "ERROR: missing project directory:"
        echo "$d"
        echo ""
        echo "PROJECT_ROOT_NOT_CONFIRMED"
        exit 1
    fi
done

echo "✓ Core project directories detected."
echo ""

# ============================================================
# 2. Confirm frozen Step19L
# ============================================================

echo "===== 2. CHECK FROZEN STEP19L ====="
echo ""

STEP19L=$(find 07_scripts \
  -maxdepth 1 \
  -type f \
  -iname '*19L*' \
  -print \
  | head -1 || true)

if [ -z "$STEP19L" ]; then
    echo "ERROR: Step19L script not found."
    echo "PROJECT_ROOT_NOT_CONFIRMED"
    exit 1
fi

echo "✓ Step19L found:"
echo "$STEP19L"
echo ""

# ============================================================
# 3. Confirm frozen Step20
# ============================================================

echo "===== 3. CHECK FROZEN STEP20 ====="
echo ""

STEP20=$(find 07_scripts \
  -maxdepth 1 \
  -type f \
  -iname '*20*' \
  -print \
  | head -1 || true)

if [ -z "$STEP20" ]; then
    echo "ERROR: Step20 script not found."
    echo "PROJECT_ROOT_NOT_CONFIRMED"
    exit 1
fi

echo "✓ Step20 found:"
echo "$STEP20"
echo ""

# ============================================================
# 4. Confirm Step21
# ============================================================

echo "===== 4. CHECK STEP21 ====="
echo ""

STEP21=$(find 07_scripts \
  -maxdepth 1 \
  -type f \
  -iname '*21*' \
  -print \
  | head -1 || true)

if [ -z "$STEP21" ]; then
    echo "ERROR: Step21 script not found."
    echo "PROJECT_ROOT_NOT_CONFIRMED"
    exit 1
fi

echo "✓ Step21 found:"
echo "$STEP21"
echo ""

# ============================================================
# 5. Automatically create Step22A directory structure
# ============================================================

echo "===== 5. CREATE STEP22A DIRECTORIES ====="
echo ""

mkdir -p \
  "00_metadata/OMIX005710" \
  "01_raw/OMIX005710" \
  "01_raw/OMIX005710/extracted" \
  "02_processed/OMIX005710" \
  "03_objects/OMIX005710" \
  "04_results/external_validation/Step22A" \
  "05_figures/external_validation/Step22A" \
  "06_tables/external_validation/Step22A"

echo "✓ Step22A directories created."
echo ""

# ============================================================
# 6. Record frozen pre-Step22A project state
# ============================================================

echo "===== 6. RECORD FROZEN ORIGINAL PROJECT STATE ====="
echo ""

python3 - <<'PY'
import os
import csv
from datetime import datetime

roots = [
    "07_scripts",
    "04_results/cross_dataset",
    "05_figures/cross_dataset",
    "06_tables/cross_dataset",
]

keywords = (
    "step19",
    "step20",
    "step21",
    "19l_",
    "20_",
    "21a_",
    "21b_",
    "21c_",
    "21d_",
    "21e_",
)

rows = []

for root in roots:

    if not os.path.exists(root):
        continue

    for dirpath, dirnames, filenames in os.walk(root):

        for filename in filenames:

            path = os.path.join(dirpath, filename)
            low = path.lower()

            if not any(k in low for k in keywords):
                continue

            stat = os.stat(path)

            rows.append({
                "path": path,
                "size_bytes": stat.st_size,
                "mtime": datetime.fromtimestamp(
                    stat.st_mtime
                ).isoformat(timespec="seconds")
            })

rows = sorted(
    rows,
    key=lambda x: x["path"].lower()
)

outfile = (
    "00_metadata/OMIX005710/"
    "STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
)

with open(
    outfile,
    "w",
    newline="",
    encoding="utf-8-sig"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "path",
            "size_bytes",
            "mtime"
        ]
    )

    writer.writeheader()
    writer.writerows(rows)

print("Frozen Step19-21 files recorded:", len(rows))
print("Saved:")
print(outfile)
PY

echo ""

# ============================================================
# 7. Create protection record
# ============================================================

cat > \
"00_metadata/OMIX005710/STEP22A_ANALYSIS_GUARDRAILS.txt" \
<<'TXT'
STEP22A THIRD-COHORT EXTERNAL VALIDATION

External dataset:
OMIX005710 / PRJCA016745

Rules:

1. Run from the original ESCC project root.

2. Step19 outputs are FROZEN.

3. Step20 outputs are FROZEN.

4. Step21 outputs are FROZEN.

5. Existing manuscript is FROZEN.

6. Existing submission package is FROZEN.

7. GSE197677 will NOT be rerun.

8. GSE221561 will NOT be rerun.

9. OMIX005710 is used only as an independent external
   validation cohort.

10. CORE2_STRONG definition must NOT be changed.

11. CORE4 definition must NOT be changed.

12. No new signature may be optimized using OMIX005710.

13. All OMIX005710 outputs must remain inside the
    Step22A external-validation branch.

14. Manuscript must NOT be modified until Step22A
    external validation has been completed.
TXT

echo "✓ Guardrail file created."
echo ""

# ============================================================
# 8. Show new directory structure
# ============================================================

echo "===== 7. STEP22A DIRECTORY STRUCTURE ====="
echo ""

for d in \
  "00_metadata/OMIX005710" \
  "01_raw/OMIX005710" \
  "02_processed/OMIX005710" \
  "03_objects/OMIX005710" \
  "04_results/external_validation/Step22A" \
  "05_figures/external_validation/Step22A" \
  "06_tables/external_validation/Step22A"
do
    if [ -d "$d" ]; then
        echo "✓ $d"
    else
        echo "✗ $d"
        exit 1
    fi
done

# ============================================================
# 9. Disk space
# ============================================================

echo ""
echo "===== 8. AVAILABLE DISK SPACE ====="
echo ""

df -h "$PROJECT"

# ============================================================
# Final
# ============================================================

echo ""
echo "============================================================"
echo "STEP22A-0 COMPLETE"
echo "============================================================"
echo ""
echo "PROJECT ROOT CONFIRMED: YES"
echo "STEP19L VISIBLE: YES"
echo "STEP20 VISIBLE: YES"
echo "STEP21 VISIBLE: YES"
echo "STEP22A DIRECTORIES CREATED: YES"
echo ""
echo "FROZEN ORIGINAL OUTPUTS MODIFIED: NO"
echo "MANUSCRIPT MODIFIED: NO"
echo "EXTERNAL DATA DOWNLOADED: NO"
echo ""
echo "READY FOR STEP22A-1"
echo ""
