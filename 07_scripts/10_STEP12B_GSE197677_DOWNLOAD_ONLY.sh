#!/bin/bash

set -u

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"
DATA_DIR="$PROJECT/01_raw/GEO/GSE197677"
LOG_DIR="$PROJECT/08_logs"
TARGETS="$PROJECT/00_metadata/05_GSE197677_download_targets.tsv"

mkdir -p "$DATA_DIR" "$LOG_DIR"
cd "$PROJECT" || exit 1

echo "============================================================"
echo "STEP 12B ONLY: DOWNLOAD GSE197677"
echo "============================================================"
echo ""

# 1. Extract targets from plan
python3 - <<'PY'
import csv, os, sys

plan = "00_metadata/01_GEO_download_plan.csv"
out  = "00_metadata/05_GSE197677_download_targets.tsv"

if not os.path.exists(plan):
    sys.exit(f"ERROR: cannot find {plan}")

rows = []
with open(plan, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r.get("GSE", "").strip() == "GSE197677":
            filename = (r.get("filename") or "").strip()
            url = (r.get("source_url") or "").strip()
            if filename and url:
                rows.append((filename, url))

if not rows:
    sys.exit("ERROR: no GSE197677 targets found")

with open(out, "w") as f:
    for fn, u in rows:
        f.write(f"{fn}\t{u}\n")

print("GSE197677 download targets:")
for fn, u in rows:
    print(f" - {fn}")
PY

# 2. Download function
download_one () {
    filename="$1"
    url="$2"
    final="$DATA_DIR/$filename"
    part="$final.part"

    echo ""
    echo "------------------------------------------------------------"
    echo "FILE: $filename"
    echo "------------------------------------------------------------"

    if [ -f "$final" ]; then
        echo "Already complete:"
        ls -lh "$final"
        return 0
    fi

    echo "Starting/resuming download..."
    curl -L --fail --retry 20 --retry-delay 10 --connect-timeout 60 -C - -o "$part" "$url" \
        > "$LOG_DIR/GSE197677_${filename}.curl.log" 2>&1 &
    CURL_PID=$!
    echo "curl PID: $CURL_PID"

    while kill -0 "$CURL_PID" 2>/dev/null; do
        if [ -f "$part" ]; then
            SIZE=$(stat -f%z "$part" 2>/dev/null || echo 0)
            python3 -c "
import sys, datetime
size = int(sys.argv[1])
print(f'[{datetime.datetime.now().strftime(\"%H:%M:%S\")}] {sys.argv[2]}: {size/1024**2:,.1f} MB ({size/1024**3:.2f} GB)')
" "$SIZE" "$filename"
        fi
        sleep 60
    done

    wait "$CURL_PID"
    status=$?

    if [ "$status" -ne 0 ]; then
        echo "Download failed. Partial file preserved. Re-run to resume."
        return 1
    fi

    if [ ! -f "$part" ]; then
        echo "ERROR: .part file missing."
        return 1
    fi

    mv "$part" "$final"
    echo "DOWNLOAD COMPLETE:"
    ls -lh "$final"
    return 0
}

# 3. Download each target
while IFS=$'\t' read -r filename url; do
    download_one "$filename" "$url"
    if [ $? -ne 0 ]; then
        echo "STEP 12B PAUSED - re-run to resume"
        exit 1
    fi
done < "$TARGETS"

# 4. Status check
echo ""
echo "============================================================"
echo "GSE197677 DOWNLOAD STATUS"
echo "============================================================"
ls -lh "$DATA_DIR"
find "$DATA_DIR" -name "*.part" -print
echo ""
echo "STEP 12B DOWNLOAD FINISHED"
echo "STOP HERE - no further analysis"
echo "============================================================"

exit 0
