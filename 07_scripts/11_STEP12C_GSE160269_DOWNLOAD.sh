#!/bin/bash
# Step 12C: Download GSE160269 supplementary files (sequential with delays)
# NCBI rate limiting causes 997-byte error pages for concurrent downloads
# Solution: sequential downloads with 10-second delays between files

set -euo pipefail

OUTDIR="${ESCC_CAF_PROJECT_ROOT:-.}/01_raw/GEO/GSE160269"
LOGDIR="${ESCC_CAF_PROJECT_ROOT:-.}/08_logs"
LOGFILE="${LOGDIR}/gse160269_download_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUTDIR" "$LOGDIR"

echo "=== GSE160269 Download ===" | tee "$LOGFILE"
echo "Start: $(date)" | tee -a "$LOGFILE"
echo "Output: $OUTDIR" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160269/suppl/"

# Sequential list of files
FILES=(
    "GSE160269_CD45neg_UMIs.txt.gz"
    "GSE160269_CD45neg_cells.txt.gz"
    "GSE160269_CD45pos_UMIs.txt.gz"
    "GSE160269_CD45pos_cells.txt.gz"
    "GSE160269_TCR_contig_annotations.csv.gz"
    "GSE160269_UMI_matrix_Bcell.txt.gz"
    "GSE160269_UMI_matrix_Endothelial.txt.gz"
    "GSE160269_UMI_matrix_Epithelia.txt.gz"
    "GSE160269_UMI_matrix_FRC.txt.gz"
    "GSE160269_UMI_matrix_Fibroblast.txt.gz"
    "GSE160269_UMI_matrix_Myeloid.txt.gz"
    "GSE160269_UMI_matrix_Pericytes.txt.gz"
    "GSE160269_UMI_matrix_Tcell.txt.gz"
)

SUCCESS=0
FAIL=0
SKIP=0

for fname in "${FILES[@]}"; do
    outfile="${OUTDIR}/${fname}"
    
    # Skip if already downloaded and valid size
    if [[ -f "$outfile" ]]; then
        actual_size=$(stat -f%z "$outfile" 2>/dev/null || echo 0)
        if [[ "$actual_size" -gt 1000 ]]; then
            echo "[SKIP] $fname (already ${actual_size} bytes)" | tee -a "$LOGFILE"
            SKIP=$((SKIP + 1))
            continue
        fi
        echo "[RETRY] $fname (was ${actual_size} bytes, re-downloading)" | tee -a "$LOGFILE"
    fi
    
    echo "[DOWNLOAD] $fname ..." | tee -a "$LOGFILE"
    curl -L --retry 3 --retry-delay 5 -o "$outfile" "${BASE_URL}${fname}" 2>> "$LOGFILE"
    
    actual_size=$(stat -f%z "$outfile" 2>/dev/null || echo 0)
    
    if [[ "$actual_size" -gt 1000 ]]; then
        echo "  OK: $actual_size bytes" | tee -a "$LOGFILE"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "  FAIL: $actual_size bytes (likely error page)" | tee -a "$LOGFILE"
        rm -f "$outfile"
        FAIL=$((FAIL + 1))
    fi
    
    # Delay between downloads to avoid rate limiting
    echo "  Waiting 10s..." | tee -a "$LOGFILE"
    sleep 10
done

echo "" | tee -a "$LOGFILE"
echo "=== Summary ===" | tee -a "$LOGFILE"
echo "Success: $SUCCESS" | tee -a "$LOGFILE"
echo "Skipped (already OK): $SKIP" | tee -a "$LOGFILE"
echo "Failed: $FAIL" | tee -a "$LOGFILE"
echo "End: $(date)" | tee -a "$LOGFILE"
