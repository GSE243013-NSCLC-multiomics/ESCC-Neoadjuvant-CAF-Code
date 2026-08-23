#!/bin/bash
set -euo pipefail

echo "============================================================"
echo "STEP23C-QA1B-R: RERENDER PLACEMENT-FIXED MANUSCRIPT"
echo "============================================================"
echo ""

PROJECT="$(pwd)"

DOCX="04_results/manuscript_update/Step23/ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.docx"

OUTDIR="04_results/manuscript_update/Step23/STEP23C_QA1B_RENDER"

PDF="$OUTDIR/ESCC_Fibroblast_CORE2_Manuscript_STEP23C_QA1B.pdf"

SOFFICE="/Applications/LibreOffice.app/Contents/MacOS/soffice"

STATUS="$OUTDIR/STEP23C_QA1B_RENDER_STATUS.txt"


# ============================================================
# 1. Verify input / tools
# ============================================================

echo "===== 1. VERIFY INPUT / TOOLS ====="
echo ""

if [ ! -s "$DOCX" ]; then
    echo "STOP: QA1B DOCX missing:"
    echo "$DOCX"
    exit 1
fi

if [ ! -x "$SOFFICE" ]; then
    echo "STOP: LibreOffice executable missing:"
    echo "$SOFFICE"
    exit 1
fi

if ! command -v pdfinfo >/dev/null 2>&1; then
    echo "STOP: pdfinfo missing."
    exit 1
fi

if ! command -v pdftoppm >/dev/null 2>&1; then
    echo "STOP: pdftoppm missing."
    exit 1
fi

if ! python3 -c 'from PIL import Image' >/dev/null 2>&1; then
    echo "STOP: Pillow missing."
    exit 1
fi

echo "✓ QA1B DOCX found"
echo "✓ LibreOffice found"
echo "✓ poppler found"
echo "✓ Pillow found"
echo ""


# ============================================================
# 2. Clean separate QA1B render directory
# ============================================================

echo "===== 2. PREPARE RENDER DIRECTORY ====="
echo ""

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "✓ Output directory:"
echo "  $OUTDIR"
echo ""


# ============================================================
# 3. Render DOCX -> PDF
# ============================================================

echo "===== 3. DOCX -> PDF ====="
echo ""

LO_PROFILE="$(mktemp -d)"
PROFILE_URI="file://$LO_PROFILE"

"$SOFFICE" \
  "-env:UserInstallation=$PROFILE_URI" \
  --headless \
  --convert-to pdf \
  --outdir "$OUTDIR" \
  "$DOCX"

rm -rf "$LO_PROFILE"

if [ ! -s "$PDF" ]; then
    echo ""
    echo "STOP: PDF was not generated."
    exit 1
fi

echo ""
echo "✓ PDF generated"
echo ""


# ============================================================
# 4. PDF basic audit
# ============================================================

echo "===== 4. PDF BASIC AUDIT ====="
echo ""

PAGES="$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')"

if [ -z "$PAGES" ]; then
    echo "STOP: PDF page count unavailable."
    exit 1
fi

echo "PDF pages : $PAGES"
echo "PDF size  : $(du -h "$PDF" | awk '{print $1}')"
echo ""


# ============================================================
# 5. PDF -> page PNG
# ============================================================

echo "===== 5. PDF -> PNG ====="
echo ""

pdftoppm \
  -png \
  -r 140 \
  "$PDF" \
  "$OUTDIR/page"

PNG_COUNT="$(
  find "$OUTDIR" \
    -maxdepth 1 \
    -name 'page-*.png' \
    | wc -l \
    | tr -d ' '
)"

echo "PDF pages : $PAGES"
echo "PNG pages : $PNG_COUNT"

if [ "$PNG_COUNT" != "$PAGES" ]; then
    echo ""
    echo "STOP: page count mismatch."
    exit 1
fi

echo "✓ Every page rasterized"
echo ""


# ============================================================
# 6. Contact sheets
#    4 manuscript pages / contact sheet
# ============================================================

echo "===== 6. CONTACT SHEETS ====="
echo ""

python3 - "$OUTDIR" <<'PY'
from pathlib import Path
from PIL import Image, ImageDraw
import re
import sys

outdir = Path(sys.argv[1])

def page_no(path):
    m = re.search(r"page-(\d+)\.png$", path.name)
    return int(m.group(1)) if m else 999999

pages = sorted(
    outdir.glob("page-*.png"),
    key=page_no
)

if not pages:
    raise SystemExit("No PNG pages found.")

PER_SHEET = 4
THUMB_W = 1000
MARGIN = 35
LABEL_H = 55

outputs = []

for start in range(0, len(pages), PER_SHEET):

    batch = pages[start:start + PER_SHEET]

    prepared = []

    for path in batch:

        image = Image.open(path).convert("RGB")

        ratio = THUMB_W / image.width

        thumb_h = round(
            image.height * ratio
        )

        image = image.resize(
            (THUMB_W, thumb_h),
            Image.Resampling.LANCZOS
        )

        panel = Image.new(
            "RGB",
            (THUMB_W, thumb_h + LABEL_H),
            "white"
        )

        panel.paste(
            image,
            (0, LABEL_H)
        )

        draw = ImageDraw.Draw(
            panel
        )

        draw.text(
            (18, 16),
            f"Page {page_no(path)}",
            fill="black"
        )

        prepared.append(
            panel
        )

    max_h = max(
        x.height
        for x in prepared
    )

    sheet_w = (
        MARGIN
        + 2 * (THUMB_W + MARGIN)
    )

    sheet_h = (
        MARGIN
        + 2 * (max_h + MARGIN)
    )

    sheet = Image.new(
        "RGB",
        (sheet_w, sheet_h),
        "white"
    )

    for i, panel in enumerate(prepared):

        row = i // 2
        col = i % 2

        x = (
            MARGIN
            + col * (THUMB_W + MARGIN)
        )

        y = (
            MARGIN
            + row * (max_h + MARGIN)
        )

        sheet.paste(
            panel,
            (x, y)
        )

    first = page_no(
        batch[0]
    )

    last = page_no(
        batch[-1]
    )

    outfile = (
        outdir
        / f"contact_pages_{first:02d}-{last:02d}.png"
    )

    sheet.save(
        outfile,
        optimize=True
    )

    outputs.append(
        outfile
    )

print(
    "Contact sheets generated:",
    len(outputs)
)

for path in outputs:
    print(path)
PY

CONTACT_COUNT="$(
  find "$OUTDIR" \
    -maxdepth 1 \
    -name 'contact_pages_*.png' \
    | wc -l \
    | tr -d ' '
)"

echo ""
echo "Contact sheets: $CONTACT_COUNT"
echo ""


# ============================================================
# 7. Create one ZIP for optional easy upload
# ============================================================

echo "===== 7. PACKAGE CONTACT SHEETS ====="
echo ""

ZIP="04_results/manuscript_update/Step23/STEP23C_QA1B_CONTACT_SHEETS.zip"

rm -f "$ZIP"

zip -j \
  "$ZIP" \
  "$OUTDIR"/contact_pages_*.png \
  >/dev/null

echo "✓ Contact-sheet ZIP:"
echo "  $ZIP"
echo ""


# ============================================================
# 8. Final render lock
# ============================================================

DOCX_SHA="$(
  shasum -a 256 "$DOCX" \
  | awk '{print $1}'
)"

PDF_SHA="$(
  shasum -a 256 "$PDF" \
  | awk '{print $1}'
)"

cat > "$STATUS" <<STATUS_EOF
STEP23C-QA1B RENDER STATUS

Input DOCX:
$DOCX

Input DOCX SHA256:
$DOCX_SHA

Rendered PDF:
$PDF

Rendered PDF SHA256:
$PDF_SHA

PDF pages:
$PAGES

PNG pages:
$PNG_COUNT

Contact sheets:
$CONTACT_COUNT

STATUS:
QA1B_RERENDER_COMPLETE_READY_FOR_VISUAL_REVIEW
STATUS_EOF


echo "============================================================"
echo "STEP23C-QA1B-R COMPLETE"
echo "============================================================"
echo ""

echo "INPUT DOCX                 : QA1B placement-fixed copy"
echo "DOCX MODIFIED DURING RENDER: NO"
echo ""

echo "PDF RENDER                 : PASS"
echo "PDF PAGES                  : $PAGES"
echo "PAGE PNGs                  : $PNG_COUNT / $PAGES"
echo "CONTACT SHEETS             : $CONTACT_COUNT"
echo ""

echo "QA DIRECTORY:"
echo "  $OUTDIR"
echo ""

echo "UPLOAD ZIP:"
echo "  $ZIP"
echo ""

echo "STATUS:"
echo "QA1B_RERENDER_COMPLETE_READY_FOR_VISUAL_REVIEW"
echo ""

echo "DO NOT EDIT OR FINALIZE THE MANUSCRIPT YET."
echo ""
