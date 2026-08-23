#!/usr/bin/env python3

# ============================================================
# STEP22A-7E
# FINAL PROJECT-LEVEL ARCHIVAL / INDEX
#
# NON-ANALYTIC CLOSEOUT ONLY.
#
# DOES NOT:
#   - open RDS
#   - read expression values
#   - recompute any score
#   - run statistics
#   - alter prior output
#   - modify manuscript
#
# Produces:
#   STEP22A_FINAL_PROJECT_README.md
#   STEP22A_FINAL_FILE_MANIFEST.csv
#   STEP22A_FINAL_CHECKSUMS.sha256
#   STEP22A_FINAL_CLOSEOUT.txt
# ============================================================

from pathlib import Path
import csv
import hashlib
import datetime
import sys


PROJECT = Path.cwd()

META = PROJECT / "00_metadata" / "OMIX005710"
RESULT = PROJECT / "04_results" / "external_validation" / "Step22A"
TABLE = PROJECT / "06_tables" / "external_validation" / "Step22A"
SCRIPTS = PROJECT / "07_scripts"
LOGS = PROJECT / "08_logs"

FINAL_7D_LOCK = (
    META
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_LOCK.csv"
)

FINAL_7D_DOSSIER = (
    META
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_DOSSIER.txt"
)

FINAL_7D_SUMMARY = (
    RESULT
    / "STEP22A_07D_FINAL_EXTERNAL_VALIDATION_SUMMARY.txt"
)

REPORTING_LANGUAGE = (
    RESULT
    / "STEP22A_07D_REPORTING_LANGUAGE.txt"
)

README_OUT = (
    META
    / "STEP22A_FINAL_PROJECT_README.md"
)

MANIFEST_OUT = (
    META
    / "STEP22A_FINAL_FILE_MANIFEST.csv"
)

CHECKSUM_OUT = (
    META
    / "STEP22A_FINAL_CHECKSUMS.sha256"
)

CLOSEOUT_OUT = (
    META
    / "STEP22A_FINAL_CLOSEOUT.txt"
)


print()
print("=" * 78)
print("STEP22A-7E: FINAL PROJECT ARCHIVAL")
print("=" * 78)
print()


# ============================================================
# Helpers
# ============================================================

def fail(msg):
    print()
    print("ERROR:")
    print(msg)
    print()
    raise SystemExit(1)


def clean(x):
    return str(x or "").replace("\r", "").strip()


def sha256_file(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        while True:
            block = f.read(1024 * 1024)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def read_csv(path):
    with open(
        path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:
        return list(csv.DictReader(f))


def relative(path):
    try:
        return str(path.relative_to(PROJECT))
    except Exception:
        return str(path)


# ============================================================
# 1. Verify final Step22A-7D lock
# ============================================================

print("===== 1. VERIFY FINAL 7D LOCK =====")
print()

for path in [
    FINAL_7D_LOCK,
    FINAL_7D_DOSSIER,
    FINAL_7D_SUMMARY,
    REPORTING_LANGUAGE
]:
    if not path.exists():
        fail(
            "Missing final Step22A-7D artifact:\n"
            + str(path)
        )

lock_rows = read_csv(
    FINAL_7D_LOCK
)

lock = {
    clean(x.get("item")):
        clean(x.get("value"))
    for x in lock_rows
}

expected = {
    "primary_external_validation":
        "PRIMARY_ESCC_ONLY_VALIDATION_UNAVAILABLE_HISTOLOGY_UNRESOLVED",

    "secondary_histology_sensitivity":
        "SECONDARY_HISTOLOGY_SENSITIVITY_DIRECTION_ROBUSTLY_OPPOSITE",

    "hypothetical_EASC_scenarios_completed":
        "22_OF_22",

    "scenario_selection":
        "NONE",

    "CORE2_mean_negative_scenarios":
        "22_OF_22",

    "CORE4_mean_negative_scenarios":
        "22_OF_22",

    "paired_LOPO_runs":
        "105",

    "paired_LOPO_mean_negative":
        "105_OF_105",

    "posthoc_CAF_subtype_explanation_performed":
        "NO",

    "CORE4_used_as_rescue_endpoint":
        "NO",

    "final_lock_status":
        "STEP22A_OMIX005710_EXTERNAL_VALIDATION_CONCLUSION_FROZEN"
}

for key, value in expected.items():

    observed = lock.get(
        key,
        "<missing>"
    )

    print(
        f"{key:<48s}: {observed}"
    )

    if observed != value:
        fail(
            f"Final lock mismatch: {key}\n"
            f"Expected: {value}\n"
            f"Observed: {observed}"
        )

print()
print("✓ Step22A-7D final conclusion verified.")
print()


# ============================================================
# 2. Collect Step22A project artifacts
#
# We deliberately DO NOT include original raw external data
# or unrelated project files.
# ============================================================

print("===== 2. COLLECT STEP22A ARTIFACTS =====")
print()

candidate_files = []

# Metadata / locks / audit files.
if META.exists():
    candidate_files.extend(
        x
        for x in META.iterdir()
        if (
            x.is_file()
            and
            x.name.startswith("STEP22A")
        )
    )

# Results.
if RESULT.exists():
    candidate_files.extend(
        x
        for x in RESULT.iterdir()
        if (
            x.is_file()
            and
            x.name.startswith("STEP22A")
        )
    )

# Tables.
if TABLE.exists():
    candidate_files.extend(
        x
        for x in TABLE.iterdir()
        if (
            x.is_file()
            and
            x.name.startswith("STEP22A")
        )
    )

# Scripts.
if SCRIPTS.exists():
    candidate_files.extend(
        x
        for x in SCRIPTS.iterdir()
        if (
            x.is_file()
            and
            (
                x.name.startswith("22A_")
                or
                x.name.startswith("STEP22A")
            )
        )
    )

# Logs/session info.
if LOGS.exists():
    candidate_files.extend(
        x
        for x in LOGS.iterdir()
        if (
            x.is_file()
            and
            "STEP22A" in x.name
        )
    )

# Remove files being generated by this archival run.
exclude = {
    README_OUT.resolve(),
    MANIFEST_OUT.resolve(),
    CHECKSUM_OUT.resolve(),
    CLOSEOUT_OUT.resolve()
}

files = sorted(
    {
        x.resolve()
        for x in candidate_files
        if x.resolve() not in exclude
    },
    key=lambda p: str(p)
)

print(
    "Step22A files found:",
    len(files)
)

print()


# ============================================================
# 3. Build manifest + SHA256
# ============================================================

print("===== 3. GENERATE MANIFEST / CHECKSUMS =====")
print()

manifest_rows = []
checksum_lines = []

for i, path in enumerate(files, start=1):

    stat = path.stat()

    digest = sha256_file(
        path
    )

    rel = relative(
        path
    )

    manifest_rows.append({
        "index":
            i,

        "path":
            rel,

        "filename":
            path.name,

        "size_bytes":
            stat.st_size,

        "modified_time":
            datetime.datetime.fromtimestamp(
                stat.st_mtime
            ).isoformat(
                timespec="seconds"
            ),

        "sha256":
            digest
    })

    checksum_lines.append(
        f"{digest}  {rel}"
    )


with open(
    MANIFEST_OUT,
    "w",
    newline="",
    encoding="utf-8-sig"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "index",
            "path",
            "filename",
            "size_bytes",
            "modified_time",
            "sha256"
        ]
    )

    writer.writeheader()
    writer.writerows(
        manifest_rows
    )


with open(
    CHECKSUM_OUT,
    "w",
    encoding="utf-8"
) as f:

    for line in checksum_lines:
        f.write(
            line
            + "\n"
        )


print(
    "Manifest entries:",
    len(
        manifest_rows
    )
)

print(
    "SHA256 entries  :",
    len(
        checksum_lines
    )
)

print()


# ============================================================
# 4. Final project README
# ============================================================

print("===== 4. WRITE FINAL PROJECT README =====")
print()

readme = """# STEP22A — OMIX005710 External Validation

## Status

**CLOSED / FROZEN**

Final scientific status:

`PRIMARY_ESCC_ONLY_UNAVAILABLE__SECONDARY_DIRECTION_ROBUSTLY_OPPOSITE`

No additional confirmatory analysis should be appended to this
external-validation workflow without explicitly opening a new,
separately labeled exploratory analysis.

---

## 1. Dataset and technical processing

- 27 tumor patient-timepoints
- 21 patients with tumor observations
- True all-27 common feature universe: 16,790 ENSG genes
- QC-retained cells: 155,501
- Frozen Fibroblast_CAF cells: 33,385
- CORE genes were excluded from QC / HVG / PCA / clustering /
  broad cell annotation to prevent circular cell selection.
- CORE genes were restored only at raw-count pseudobulk stage.

---

## 2. Frozen fibroblast technical rule

Minimum fibroblast count per patient-timepoint:

`>=20 cells`

Final technical eligibility:

- PASS: 26/27 tumor patient-timepoints
- FAIL: P15 Before = 18 fibroblast cells
- Threshold was not relaxed.
- No 18-cell rescue branch was created.
- P15 Before remains retained in raw provenance but is excluded
  from the technical reference.

Complete paired patients:

- P6
- P9
- P15
- P16
- P18
- P19

Technically evaluable paired patients:

- P6
- P9
- P16
- P18
- P19

P15 was excluded from paired analysis because its Before sample
failed the frozen >=20-cell criterion.

---

## 3. Histology limitation

The source cohort contains one EASC patient among 22 patients.

The EASC patient's individual identity could not be resolved from:

- public supplementary materials;
- public HRI patient pages;
- available OMIX/GSA-Human metadata.

Therefore a definitive patient-level **ESCC-only primary reference**
cannot be constructed without making an unsupported assumption.

Primary ESCC-only external validation is consequently:

**UNAVAILABLE**

This is a histology-identification limitation, not a computational
failure.

---

## 4. Pre-specified secondary histology sensitivity

Before reading CORE scores, 23 histology states were frozen:

- 1 unresolved descriptive baseline;
- 22 hypothetical EASC identities.

All 22 hypothetical EASC scenarios were subsequently evaluated.

No scenario was selected.

No scenario was promoted to primary.

The scenarios form a sensitivity envelope and are not independent
validation cohorts.

---

## 5. Frozen scoring implementation

Expression transformation:

`raw counts -> library-size CPM -> log2(CPM + 1)`

Z-score reference:

- OMIX005710 internal eligible tumor patient-timepoints;
- gene-wise mean and SD;
- recomputed within each frozen histology scenario.

Frozen scores:

`CORE2_STRONG = (-CDKN1B_z - GSN_z) / 2`

`CORE2_WEAK = (BGN_z + TIMP1_z) / 2`

`CORE4 = (BGN_z - CDKN1B_z - GSN_z + TIMP1_z) / 4`

Treatment effect:

`After - Before`

Positive change was the frozen replication direction.

---

## 6. Secondary external-validation result

Across all 22 hypothetical-EASC scenarios:

- CORE2_STRONG mean delta negative: **22/22**
- CORE2_STRONG median delta negative: **22/22**
- Majority paired patients negative: **22/22**
- CORE4 mean delta negative: **22/22**
- CDKN1B contribution negative: **22/22**
- GSN contribution negative: **22/22**
- BGN/TIMP1 component:
  - negative: **18/22**
  - positive: **4/22**

Thus the frozen treatment-associated direction observed in the
discovery cohorts was **not replicated** in OMIX005710's
pre-specified secondary sensitivity framework.

The observed direction was consistently opposite.

---

## 7. Technical replication audit

Independent recomputation reproduced the scoring chain within
floating-point precision.

Audited items included:

- sample/timepoint labels;
- Before vs After orientation;
- P15 exclusion;
- reference membership;
- library-size CPM;
- log2(CPM+1);
- gene-wise z-scores;
- CORE2 formula;
- paired delta calculation.

No sign, pairing, reference-membership, or score-implementation
artifact was detected.

---

## 8. Leave-one-paired-patient-out audit

Full LOPO recomputation was performed.

Total LOPO runs:

**105**

Results:

- mean CORE2 delta negative: **105/105**
- median CORE2 delta negative: **105/105**
- majority remaining patients negative: **105/105**
- all remaining paired patients negative: **105/105**

Therefore the negative secondary result is not attributable to one
paired patient.

---

## 9. Interpretation boundary

Supported conclusion:

> In the pre-specified secondary histology-sensitivity analysis,
> OMIX005710 did not replicate the frozen treatment-associated CORE
> direction observed in the discovery cohorts and instead showed a
> robustly opposite direction.

Not supported:

- definitive ESCC-only external validation;
- a claim that OMIX005710 proves the opposite biological mechanism;
- selection of a favorable histology scenario;
- post-hoc modification of the fibroblast threshold;
- CORE4 as a rescue endpoint;
- post-hoc CAF subtype explanations presented as confirmatory.

Any future investigation of the biological cause of this discordance
must be opened as a separate **exploratory / hypothesis-generating**
analysis.

---

## 10. Manuscript status

The analysis pipeline did not automatically modify the manuscript.

Suggested Results / Discussion wording is stored in:

`04_results/external_validation/Step22A/STEP22A_07D_REPORTING_LANGUAGE.txt`

---

## 11. Final archival files

- `STEP22A_FINAL_PROJECT_README.md`
- `STEP22A_FINAL_FILE_MANIFEST.csv`
- `STEP22A_FINAL_CHECKSUMS.sha256`
- `STEP22A_FINAL_CLOSEOUT.txt`

The manifest and checksum file provide a frozen inventory of the
Step22A analytical artifacts at closeout.

---

## Final status

**STEP22A EXTERNAL-VALIDATION PIPELINE CLOSED.**
"""

with open(
    README_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        readme
    )

print(
    "✓ Final README written."
)

print()


# ============================================================
# 5. Closeout lock
# ============================================================

timestamp = datetime.datetime.now().astimezone().isoformat(
    timespec="seconds"
)

manifest_sha = sha256_file(
    MANIFEST_OUT
)

checksum_list_sha = sha256_file(
    CHECKSUM_OUT
)

readme_sha = sha256_file(
    README_OUT
)

closeout = f"""STEP22A FINAL CLOSEOUT
================================================================

Closeout timestamp:
{timestamp}

ANALYTIC STATUS:
CLOSED

FINAL SCIENTIFIC STATUS:
PRIMARY_ESCC_ONLY_UNAVAILABLE__SECONDARY_DIRECTION_ROBUSTLY_OPPOSITE

PRIMARY ESCC-ONLY:
UNAVAILABLE

PRIMARY LIMITATION:
EASC patient identity unresolved in public patient-level metadata

SECONDARY HISTOLOGY SENSITIVITY:
22 / 22 hypothetical EASC scenarios completed

CORE2 DIRECTION:
NEGATIVE 22 / 22 scenarios

CORE4 DIRECTION:
NEGATIVE 22 / 22 scenarios

PAIRED LOPO:
105 / 105 mean negative
105 / 105 median negative
105 / 105 majority negative
105 / 105 all remaining pairs negative

TECHNICAL ARTIFACT DETECTED:
NO

SINGLE-PATIENT DRIVER DETECTED:
NO

SCENARIO SELECTION:
NO

THRESHOLD RELAXATION:
NO

CORE4 USED AS RESCUE:
NO

POST-HOC BIOLOGICAL EXPLANATION:
NO

MANUSCRIPT AUTOMATICALLY MODIFIED:
NO

STEP22A ARTIFACT COUNT:
{len(manifest_rows)}

FINAL MANIFEST SHA256:
{manifest_sha}

FINAL CHECKSUM-LIST SHA256:
{checksum_list_sha}

FINAL README SHA256:
{readme_sha}

FINAL RULE:
No additional confirmatory analysis may be appended to Step22A.
Any biological mechanism investigation must be opened separately
and explicitly labeled exploratory / hypothesis-generating.

FINAL CLOSEOUT STATUS:
STEP22A_EXTERNAL_VALIDATION_ARCHIVED_AND_CLOSED
"""

with open(
    CLOSEOUT_OUT,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        closeout
    )


# ============================================================
# FINAL
# ============================================================

print("=" * 78)
print("STEP22A-7E COMPLETE")
print("=" * 78)
print()

print(
    "ANALYTIC PIPELINE              : CLOSED"
)

print(
    "FINAL SCIENTIFIC STATUS       :"
)

print(
    "PRIMARY_ESCC_ONLY_UNAVAILABLE__SECONDARY_DIRECTION_ROBUSTLY_OPPOSITE"
)

print()

print(
    "STEP22A ARTIFACTS INDEXED     :",
    len(
        manifest_rows
    )
)

print()

print(
    "FINAL FILES:"
)

print(
    "  00_metadata/OMIX005710/STEP22A_FINAL_PROJECT_README.md"
)

print(
    "  00_metadata/OMIX005710/STEP22A_FINAL_FILE_MANIFEST.csv"
)

print(
    "  00_metadata/OMIX005710/STEP22A_FINAL_CHECKSUMS.sha256"
)

print(
    "  00_metadata/OMIX005710/STEP22A_FINAL_CLOSEOUT.txt"
)

print()

print(
    "NEW ANALYSIS PERFORMED         : NO"
)

print(
    "EXPRESSION DATA RECOMPUTED     : NO"
)

print(
    "CORE SCORES RECOMPUTED         : NO"
)

print(
    "MANUSCRIPT MODIFIED            : NO"
)

print()

print(
    "FINAL CLOSEOUT STATUS:"
)

print(
    "STEP22A_EXTERNAL_VALIDATION_ARCHIVED_AND_CLOSED"
)

print()

print(
    "NO FURTHER STEP22A ANALYSIS IS REQUIRED."
)

print()
