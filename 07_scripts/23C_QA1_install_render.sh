#!/bin/bash

set -e

echo "============================================================"
echo "STEP23C-QA1: INSTALL RENDERER + VISUAL QA RENDER"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# 1. Ensure Homebrew
# ------------------------------------------------------------

echo "===== 1. Check Homebrew ====="

if command -v brew >/dev/null 2>&1; then
    echo "✓ Homebrew found: $(command -v brew)"
else
    echo "Homebrew not found. Installing with official installer..."

    NONINTERACTIVE=1 /bin/bash -c \
      "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    # Apple Silicon
    if [ -x /opt/homebrew/bin/brew ]; then
        eval "$(/opt/homebrew/bin/brew shellenv)"
    fi

    # Intel Mac
    if [ -x /usr/local/bin/brew ]; then
        eval "$(/usr/local/bin/brew shellenv)"
    fi
fi

if ! command -v brew >/dev/null 2>&1; then
    echo ""
    echo "STOP: Homebrew installation did not complete."
    echo "Send me this terminal screen; do not continue manually."
    exit 1
fi

echo "✓ Homebrew ready"
echo ""

# ------------------------------------------------------------
# 2. Install LibreOffice
# ------------------------------------------------------------

echo "===== 2. Check LibreOffice ====="

SOFFICE="/Applications/LibreOffice.app/Contents/MacOS/soffice"

if [ -x "$SOFFICE" ]; then
    echo "✓ LibreOffice already installed"
else
    echo "Installing LibreOffice..."
    brew install --cask libreoffice
fi

if [ ! -x "$SOFFICE" ]; then
    echo "STOP: LibreOffice executable not found after installation."
    exit 1
fi

"$SOFFICE" --version
echo ""

# ------------------------------------------------------------
# 3. Install Poppler / pdftoppm
# ------------------------------------------------------------

echo "===== 3. Check PDF renderer ====="

if command -v pdftoppm >/dev/null 2>&1 && \
   command -v pdfinfo >/dev/null 2>&1; then

    echo "✓ Poppler tools already available"

else
    echo "Installing poppler..."
    brew install poppler
fi

command -v pdftoppm
command -v pdfinfo
echo ""

# ------------------------------------------------------------
# 4. Paths
# ------------------------------------------------------------

WORKING="04_results/manuscript_update/Step23/ESCC_Fibroblast_CORE2_Manuscript_STEP23C_WORKING.docx"

OUTDIR="04_results/manuscript_update/Step23/STEP23C_QA_RENDER"

PDF="$OUTDIR/ESCC_Fibroblast_CORE2_Manuscript_STEP23C_WORKING.pdf"

if [ ! -f "$WORKING" ]; then
    echo "STOP: working DOCX not found:"
    echo "$WORKING"
    exit 1
fi

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "===== 4. Render DOCX -> PDF ====="
echo "DOCX: $WORKING"
echo ""

# Dedicated LibreOffice profile avoids profile-lock problems.
LO_PROFILE="$(mktemp -d)"
PROFILE_URI="file://$LO_PROFILE"

"$SOFFICE" \
  "-env:UserInstallation=$PROFILE_URI" \
  --headless \
  --convert-to pdf \
  --outdir "$OUTDIR" \
  "$WORKING"

rm -rf "$LO_PROFILE"

if [ ! -s "$PDF" ]; then
    echo ""
    echo "STOP: PDF render was not produced."
    exit 1
fi

echo ""
echo "✓ PDF generated:"
echo "$PDF"
echo ""

# ------------------------------------------------------------
# 5. PDF basic QA
# ------------------------------------------------------------

echo "===== 5. PDF basic QA ====="

PAGES="$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')"

if [ -z "$PAGES" ]; then
    echo "STOP: could not determine PDF page count."
    exit 1
fi

echo "PDF pages: $PAGES"
echo "PDF size : $(du -h "$PDF" | awk '{print $1}')"
echo ""

# ------------------------------------------------------------
# 6. PDF -> page PNG
# ------------------------------------------------------------

echo "===== 6. Render every page -> PNG ====="

pdftoppm \
  -png \
  -r 120 \
  "$PDF" \
  "$OUTDIR/page"

PNG_COUNT="$(find "$OUTDIR" -maxdepth 1 -name 'page-*.png' | wc -l | tr -d ' ')"

echo "PNG pages generated: $PNG_COUNT"

if [ "$PNG_COUNT" != "$PAGES" ]; then
    echo ""
    echo "STOP: PDF page count and PNG count disagree."
    echo "PDF pages : $PAGES"
    echo "PNG pages : $PNG_COUNT"
    exit 1
fi

echo "✓ Every PDF page rasterized"
echo ""

# ------------------------------------------------------------
# 7. Ensure Pillow
# ------------------------------------------------------------

echo "===== 7. Prepare contact sheets ====="

if ! python3 -c 'from PIL import Image' >/dev/null 2>&1; then
    echo "Installing Pillow..."
    python3 -m pip install --user pillow
fi

# ------------------------------------------------------------
# 8. Generate contact sheets: 4 manuscript pages per image
# ------------------------------------------------------------

python3 - "$OUTDIR" <<'PY'
from pathlib import Path
from PIL import Image, ImageOps, ImageDraw
import re
import sys
import math

outdir = Path(sys.argv[1])

def page_number(path):
    m = re.search(r'page-(\d+)\.png$', path.name)
    return int(m.group(1)) if m else 10**9

pages = sorted(
    outdir.glob("page-*.png"),
    key=page_number
)

if not pages:
    raise SystemExit("No page PNG files found")

# Four pages per contact sheet: 2 x 2.
per_sheet = 4
thumb_w = 850
margin = 30
label_h = 45

contacts = []

for start in range(0, len(pages), per_sheet):
    batch = pages[start:start + per_sheet]

    prepared = []

    for path in batch:
        img = Image.open(path).convert("RGB")

        ratio = thumb_w / img.width
        thumb_h = round(img.height * ratio)

        img = img.resize(
            (thumb_w, thumb_h),
            Image.Resampling.LANCZOS
        )

        canvas = Image.new(
            "RGB",
            (thumb_w, thumb_h + label_h),
            "white"
        )

        canvas.paste(img, (0, label_h))

        draw = ImageDraw.Draw(canvas)
        n = page_number(path)

        draw.text(
            (15, 12),
            f"Page {n}",
            fill="black"
        )

        prepared.append(canvas)

    max_h = max(x.height for x in prepared)

    sheet_w = margin + 2 * (thumb_w + margin)
    sheet_h = margin + 2 * (max_h + margin)

    sheet = Image.new(
        "RGB",
        (sheet_w, sheet_h),
        "white"
    )

    for i, img in enumerate(prepared):
        row = i // 2
        col = i % 2

        x = margin + col * (thumb_w + margin)
        y = margin + row * (max_h + margin)

        sheet.paste(img, (x, y))

    first_page = page_number(batch[0])
    last_page = page_number(batch[-1])

    outfile = outdir / (
        f"contact_pages_{first_page:02d}-{last_page:02d}.png"
    )

    sheet.save(outfile, optimize=True)
    contacts.append(outfile)

print(f"Contact sheets generated: {len(contacts)}")
for x in contacts:
    print(x)
PY

echo ""

# ------------------------------------------------------------
# 9. Final report
# ------------------------------------------------------------

CONTACT_COUNT="$(find "$OUTDIR" -maxdepth 1 -name 'contact_pages_*.png' | wc -l | tr -d ' ')"

cat > "$OUTDIR/STEP23C_QA1_RENDER_STATUS.txt" <<EOF
STEP23C-QA1 VISUAL RENDER STATUS

Working DOCX:
$WORKING

Rendered PDF:
$PDF

PDF pages:
$PAGES

Page PNGs:
$PNG_COUNT

Contact sheets:
$CONTACT_COUNT

LibreOffice:
$("$SOFFICE" --version | head -1)

STATUS:
RENDER_COMPLETE_READY_FOR_VISUAL_QA
EOF

echo "============================================================"
echo "STEP23C-QA1 COMPLETE"
echo "============================================================"
echo ""
echo "WORKING DOCX              : $WORKING"
echo "PDF RENDER                : PASS"
echo "PDF PAGES                 : $PAGES"
echo "PAGE PNGs                 : $PNG_COUNT / $PAGES"
echo "CONTACT SHEETS            : $CONTACT_COUNT"
echo ""
echo "QA DIRECTORY:"
echo "  $OUTDIR"
echo ""
echo "STATUS:"
echo "RENDER_COMPLETE_READY_FOR_VISUAL_QA"
echo ""
echo "DO NOT FINALIZE THE MANUSCRIPT YET."
echo ""
