# Script map

This document inventories the 191 workflow scripts in `07_scripts/` so they can
be interpreted without renaming or moving any file. It is a navigation aid, not
an execution guarantee — see the reproducibility notes below and in
[README.md](README.md).

## Scope

- **191 workflow scripts total**
- **143 R**
- **36 Python**
- **12 shell**

These counts refer to the scripts tracked under `07_scripts/` at the initial
commit. On-disk local variants coexisting via `git update-index --skip-worktree`
are outside the scope of this map and are never modified here.

## How to read filenames

Scripts use a stage-prefix convention:

- **`NN` / `NNA` / `NNB` / …** — two-digit stage number, optionally suffixed to
  disambiguate multiple scripts within one stage (e.g., `17H1`, `17H2`,
  `22A_05A` vs `22A_05A_v2`).
- **`MAIN`, `FIX`, `REV`, `FREEZE`, `FINAL`, `LOCK`** — these encode
  *workflow/revision semantics*, not independent biological conclusions. They
  indicate that a script is a corrected variant, a revision, a frozen-state
  producer/consumer, or a final reconciled output of an earlier analysis. They
  should be read as provenance markers, not as separate scientific claims.

The stage prefixes (`00`–`23`) track the project's analytical evolution and the
approximate order in which results were produced.

## Stage map

The stage boundaries below are the same as those used in
[README.md](README.md). Stages 17–23 are described in additional detail because
they contain the largest number of supporting, pilot, and QA scripts.

| Stage | Scope |
|-------|-------|
| 00–03 | Project setup, GEO metadata audit, download planning, preflight. |
| 04–09 | GSE221561 acquisition, extraction/audit, Seurat construction, QC, metadata standardization, mapping validation. |
| 10–11 | GSE197677 and GSE160269 acquisition/download helpers. |
| 12–13 | GSE160269 object construction, corrected metadata, fibroblast/reference annotation. |
| 14–15 | GSE197677 lineage QC, cell-type annotation, tumor eligibility correction. |
| 16 | Unified cross-dataset lineage/eligibility tables. |
| 17 | Malignant-state/CNV analyses and supporting CopyKAT/inferCNV audits (many `17H*` per-sample CNV pilots and batch runs). |
| 18 | Cross-dataset consistency and frozen CNV/reference layers (incl. CopyKAT reference-sensitivity pilots). |
| 19 | Frozen cross-dataset analytical core: pseudobulk inventory (19A), DGE/QC (19B), GSE197677 DE and sensitivity (19C), replication and orientation correction (19D), hallmark/pathway robustness (19E/19F), epithelial localization/mixed-program audits (19G–19I), CORE definitions and direction correction (19J), fibroblast core/abundance/heterogeneity/regulator audits (19K–19N), targeted ligand-receptor signaling audits and reconciliation (19O–19P), and final signaling synthesis (19Q). |
| 20 | Global project status, step status, evidence hierarchy, and optional-branch tracking. |
| 21 | Frozen result-report/manuscript/figure/table/export assembly (incl. DOCX/PDF consistency tooling). |
| 22A | External validation: frozen score source, 27-sample tumor QC, histology gate, technical CORE2/CORE4 replication audits, and final governance/conclusion. |
| 23 | Manuscript integration, evidence-map freeze, QA edits, verified source-reference insertion, and human-review gates. |

### Stage 17–23 detail (why so many scripts)

- **Stage 17 / 18** are dominated by per-sample CNV pilots (`17H*`), method
  comparisons (inferCNV vs CopyKAT), and reference-sensitivity audits. Many
  files are pilots or per-sample batch runs rather than final entrypoints.
- **Stage 19** is the analytical core; it is subdivided by sub-letter
  (`19A`, `19B`, …, `19Q`). Sub-steps such as `19O_FIX221` and `19J_CORR*`
  are correction/reconciliation variants of earlier sub-steps.
- **Stage 21** mixes final-assembly scripts with revision variants
  (`21D_REV1*`, `21D_REV2*`) and the DOCX/PDF export/consistency utilities.
- **Stage 22A** contains both the external-validation analyses and the
  freeze/verification/histology-gate scripts that lock the scoring source.
- **Stage 23** is almost entirely manuscript-integration and QA tooling
  (context maps, safe-locked edits, review packets, human-decision freezes);
  these are publication utilities, not new analyses.

## Canonical vs supporting scripts

This map deliberately does **not** designate a fixed set of "canonical
entrypoints" for every stage. When evidence (frozen/final markers, downstream
consumption, or explicit freeze records) is insufficient, no claim is made. The
191 scripts fall into workflow-role categories identifiable from filenames:

- **primary / frozen / final analytical scripts** — e.g., those carrying
  `STEP*_FROZEN*`, `*_FREEZE*`, `*_FINAL*`, or the core sub-step letters
  (19A–19Q). These produce or consume the frozen evidence tables.
- **pilot / sensitivity scripts** — e.g., filenames containing `PILOT`,
  `SENSITIVITY`, `LOSO`, or per-sample CNV pilots (many `17H*`, `18E*`,
  `19C12`, `19G*`, `19K*`, `19L*`, `19N*`, `19P*`).
- **audit / QC scripts** — e.g., filenames containing `AUDIT`, `QC`,
  `PREFILGHT`, `REVIEW`, `VERIFY`, `INVENTORY`.
- **fix / reconciliation scripts** — e.g., filenames containing `FIX`,
  `CORR`, `RECONCILE`, `CORRECTED`, `DIRECTION`.
- **render / export / manuscript utilities** — e.g., the `*.py` export/assembly
  tools and the `23C_QA*` / `23C_*_render*` manuscript-integration helpers.

Filenames such as `PILOT`, `AUDIT`, `REVIEW`, `FIX`, `REV`, `FREEZE`, `FINAL`,
and `LOCK` carry **workflow provenance**. They should *not* be deleted merely
because they are not the final analysis entrypoints — they document how frozen
results were reached and are part of the auditable record.

## Execution

The numeric stage prefixes reflect the project's evolution and an *approximate*
dependency order. They are **not** automatically equivalent to a single
sequential executable pipeline: scripts read inputs produced by earlier stages
(many of which are intermediate objects excluded from this repository), and
several stages contain parallel or exploratory variants rather than one linear
chain.

No input/output dependency is asserted here beyond what is explicitly stated in
the scripts themselves (file paths, design matrices, and freeze records under
`06_tables/` and `00_metadata/`). Do not infer additional wiring from filename
proximity alone.

## Reproducibility

- [README.md](README.md) — overview, environment, and reproduction scope
- [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) — public data sources and exclusions
- [CODE_AVAILABILITY.md](CODE_AVAILABILITY.md) — code-release status and placeholders
- [RESOURCE_PROVENANCE.md](RESOURCE_PROVENANCE.md) — external resource provenance and to-be-verified items
