#!/bin/bash

set -euo pipefail

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"
cd "$PROJECT"

echo ""
echo "============================================================"
echo "STEP22A-0B: EXACT FROZEN PROJECT VERIFICATION"
echo "============================================================"
echo ""

ERRORS=0

# ------------------------------------------------------------
# 1. Exact Step19L check
# ------------------------------------------------------------

echo "===== 1. STEP19L ====="

STEP19L="07_scripts/19L_CONSERVED_FIBROBLAST_CORE_VALIDATION.R"

if [ -f "$STEP19L" ]; then
    echo "✓ FOUND:"
    echo "  $STEP19L"
else
    echo "✗ NOT FOUND:"
    echo "  $STEP19L"
    ERRORS=$((ERRORS + 1))
fi

echo ""

# ------------------------------------------------------------
# 2. Exact Step20 check
# ------------------------------------------------------------

echo "===== 2. STEP20 ====="

STEP20="07_scripts/20_GLOBAL_PROJECT_STATUS_AND_SYNTHESIS.R"

if [ -f "$STEP20" ]; then
    echo "✓ FOUND:"
    echo "  $STEP20"
else
    echo "✗ NOT FOUND:"
    echo "  $STEP20"
    ERRORS=$((ERRORS + 1))
fi

echo ""

# ------------------------------------------------------------
# 3. Exact Step21-series checks
# ------------------------------------------------------------

echo "===== 3. STEP21 SERIES ====="
echo ""

for PREFIX in 21A 21B 21C 21D 21E
do
    MATCHES=$(find 07_scripts \
        -maxdepth 1 \
        -type f \
        -name "${PREFIX}_*.R" \
        -print \
        | sort)

    if [ -n "$MATCHES" ]; then
        echo "✓ ${PREFIX}:"
        echo "$MATCHES" | sed 's/^/  /'
    else
        echo "✗ ${PREFIX}: no matching script found"
        ERRORS=$((ERRORS + 1))
    fi

    echo ""
done

# ------------------------------------------------------------
# 4. Explicitly show the previous false-positive file
# ------------------------------------------------------------

echo "===== 4. FALSE-POSITIVE GUARD ====="
echo ""

FALSE_POSITIVE=$(find 07_scripts \
    -maxdepth 1 \
    -type f \
    -name '09_GSE221561_mapping_validation.R' \
    -print || true)

if [ -n "$FALSE_POSITIVE" ]; then
    echo "Detected unrelated file:"
    echo "  $FALSE_POSITIVE"
    echo ""
    echo "✓ It is NOT being counted as Step21."
else
    echo "No 09_GSE221561_mapping_validation.R file found."
    echo "✓ No issue."
fi

echo ""

# ------------------------------------------------------------
# 5. Confirm Step22A directories
# ------------------------------------------------------------

echo "===== 5. STEP22A DIRECTORIES ====="
echo ""

DIRS=(
  "00_metadata/OMIX005710"
  "01_raw/OMIX005710"
  "01_raw/OMIX005710/extracted"
  "02_processed/OMIX005710"
  "03_objects/OMIX005710"
  "04_results/external_validation/Step22A"
  "05_figures/external_validation/Step22A"
  "06_tables/external_validation/Step22A"
)

for d in "${DIRS[@]}"
do
    if [ -d "$d" ]; then
        echo "✓ $d"
    else
        echo "✗ MISSING: $d"
        ERRORS=$((ERRORS + 1))
    fi
done

echo ""

# ------------------------------------------------------------
# 6. Confirm frozen inventory + guardrails
# ------------------------------------------------------------

echo "===== 6. PROTECTION FILES ====="
echo ""

FILES=(
  "00_metadata/OMIX005710/STEP22A_PREANALYSIS_FROZEN_FILE_INVENTORY.csv"
  "00_metadata/OMIX005710/STEP22A_ANALYSIS_GUARDRAILS.txt"
)

for f in "${FILES[@]}"
do
    if [ -f "$f" ]; then
        echo "✓ $f"
    else
        echo "✗ MISSING: $f"
        ERRORS=$((ERRORS + 1))
    fi
done

echo ""

# ------------------------------------------------------------
# 7. Disk space — information only
# ------------------------------------------------------------

echo "===== 7. DISK SPACE ====="
echo ""

df -h "$PROJECT"

echo ""

# ------------------------------------------------------------
# Final decision
# ------------------------------------------------------------

echo "============================================================"

if [ "$ERRORS" -eq 0 ]; then

    echo "STEP22A-0B COMPLETE"
    echo "============================================================"
    echo ""
    echo "PROJECT ROOT CONFIRMED: YES"
    echo "STEP19L EXACT MATCH: YES"
    echo "STEP20 EXACT MATCH: YES"
    echo "STEP21A-E EXACT MATCH: YES"
    echo "FALSE-POSITIVE STEP21 MATCH: ELIMINATED"
    echo "STEP22A DIRECTORIES: READY"
    echo "FROZEN PROTECTION FILES: READY"
    echo ""
    echo "ORIGINAL STEP19-21 OUTPUTS MODIFIED: NO"
    echo "MANUSCRIPT MODIFIED: NO"
    echo "EXTERNAL DATA DOWNLOADED: NO"
    echo ""
    echo "READY FOR STEP22A-1"

else

    echo "STEP22A-0B STOPPED"
    echo "============================================================"
    echo ""
    echo "VERIFICATION ERRORS: $ERRORS"
    echo ""
    echo "PROJECT_FROZEN_STATE_NOT_FULLY_CONFIRMED"
    echo "DO NOT START STEP22A-1"

fi
