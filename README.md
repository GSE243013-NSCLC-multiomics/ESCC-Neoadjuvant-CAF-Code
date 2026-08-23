# ESCC Neoadjuvant Epithelial–Myeloid–CAF Analysis

## Overview

This repository contains the reproducible analysis code, selected derived
results, and provenance records for a cross-dataset study of the conserved
fibroblast **CORE2** response in esophageal squamous cell carcinoma (ESCC)
under neoadjuvant treatment. The analysis integrates single-cell and
bulk-derived pseudobulk transcriptomes from multiple independent public
datasets (see [Public datasets](#public-datasets)), performs cross-dataset
differential expression and replication under a frozen analysis plan, and
builds a layered evidence hierarchy for the reported biological claims.

The repository is organized as a *code-plus-selected-output* release: it ships
the scripts, the final and intermediate tables that support the manuscript,
the publication figures, metadata/provenance records, and the R environment
lockfile. Large raw datasets and heavy intermediate objects are intentionally
**not** redistributed here (see [Data availability](#data-availability)).

## Study design

- **Disease/system:** ESCC tumor microenvironment, with focus on fibroblast
  (CAF) state changes induced by neoadjuvant chemotherapy.
- **Primary outcome:** the conserved fibroblast CORE2 signature response —
  a state defined and frozen in earlier project steps —
  with CORE4 and epithelial/immune compartments investigated as
  interactive or comparative axes.
- **Analysis layers:**
  - per-dataset preprocessing, QC, metadata/histology audits and Seurat
    object construction;
  - a frozen Step 19 cross-dataset pipeline: epithelial pseudobulk
    construction, DGEList/TMM normalization, QC decision tables,
    within-dataset PCA/MDS, sample outlier and covariate audits,
    GSE197677 DE, leave-one-sample-out sensitivity, GSE221561 replication
    with direction concordance, and final freeze tables;
  - hallmark pathway analysis (preranked fgsea/msigdbr) with robustness
    audit;
  - a ligand-receptor (LR) signaling audit using a small, script-preaspecified
    curated LR resource (see `RESOURCE_PROVENANCE.md`), with explicit
    documentation of resource limitations (100 curated pairs; no external
    CellChatDB/NicheNet/OmniPath/LIANA resources available at analysis time;
    receptor evaluability constrained by the frozen 386-gene fibroblast
    pseudobulk; no ligand-target matrix);
  - Step 22A external validation of the CORE2 technical score across an
    independent 27-sample paired-tumor collection, driven by a frozen
    original-source scoring reference and histology gate;
  - manuscript assembly (Steps 21/23) with DOCX/PDF export consistency checks.

Key design principle: the analysis plan and evidence-frozen tables are
versioned and auditable. Each claim at the end of the pipeline is traceable to
specific script steps and frozen table files. Branches that were explored but
not adopted (e.g., Step 20 optional branches) are recorded rather than
silently dropped.

## Public datasets

The following public data sources were used or considered in this study.
Raw data are not redistributed in this repository; they remain available from
their original repositories/providers.

| ID / name             | Role in this study                                                    | Status in this release |
|-----------------------|------------------------------------------------------------------------|------------------------|
| GSE221561             | Independent ESCC scRNA-seq cohort; download, QC, and Seurat build (Steps 01-09); pseudobulk DE replication (Steps 19D+) | Used, accession recorded in metadata |
| GSE197677             | ESCC scRNA-seq cohort; lineage/broad-celltype annotation, corrected tumor eligibility; discovery DE in Step 19 | Used, accession recorded in metadata |
| GSE160269             | ESCC scRNA-seq cohort; reference for fibroblast annotation and sensitivity audits (Steps 12-13, 19G) | Used, accession recorded in metadata |
| OMIX005710            | Independent ESCC scRNA-seq cohort (CNCB-NGDC); secondary external-validation / histology-sensitivity analysis (Step 22A) | Used, accession recorded in metadata |
| TCGA-ESCA             | Referenced context for ESCC gene-expression-level survival analyses (keyword-/source audit in Step 22A metadata) | Considered/referenced |
| GSE53625              | Optional external bulk-validation branch considered in Step 20 planning | Considered (optional branch; not adopted) |
| DepMap                | Dependency/projection context considered for optional branches | Considered |

Datasets are referenced in `00_metadata/` by their series accessions as
downloaded/processed at analysis time. Exact retrieval versions, URLs, and
DOIs are listed per-release in `DATA_AVAILABILITY.md` and
`RESOURCE_PROVENANCE.md` using the requested public-release placeholders until
final publication details are confirmed.

**Project root / paths.** Analysis scripts resolve the project root via the
`ESCC_CAF_PROJECT_ROOT` environment variable, defaulting to the current working
directory. Reproduce by either running the scripts from the repository root, or
by exporting `ESCC_CAF_PROJECT_ROOT` to the directory that contains the expected
input/working subdirectories (e.g. `01_raw/`, `00_metadata/`, `06_tables/`,
`08_logs/`). Public datasets must be downloaded from their original sources into
the appropriate input locations; they are not bundled in this repository.

## Repository structure

```
00_metadata/   Sample metadata, GEO series audit, download plans/manifests,
               OMIX source tables and frozen provenance records (161 files)
01_raw/        NOT TRACKED — raw source data (kept locally)
02_processed/  NOT TRACKED — processed intermediates (kept locally)
03_objects/    NOT TRACKED — large R objects (kept locally)
04_results/    Final external-validation result summaries (5 files)
05_figures/    Publication figures incl. Step21D manuscript finals (70 files)
06_tables/     All intermediate and frozen result tables (351 files)
07_scripts/    191 workflow scripts (143 R, 36 Python, 12 shell)
08_logs/       Logs, sessionInfo, and audit artifacts (113 files)
09_tmp/        NOT TRACKED — scratch space (kept locally)
renv/          R environment scaffolding; `renv.lock` records R 4.6.1, Bioconductor 3.23, and the CRAN repository only — it does not pin all analysis packages (a complete package lock should be generated/verified before public release).
```

Scripts, figures and tables are cross-referenced by Step number; the frozen
"CORE2 effect evidence" of the manuscript is assembled from `06_tables`,
`05_figures`, and `08_logs` at Steps 21D/21E.

## Software environment

- R `4.6.1`, Bioconductor `3.23` (recorded in `renv.lock`).
- Package management with `renv` (`renv.lock`): CRAN repos pinned to
  `https://cloud.r-project.org`.
- Python 3 used for Step 22A / Step 23 manuscript-integration tooling
  (standard-library plus the specific imports declared in each script).
- Key R analysis packages include `edgeR`, `Seurat`, `fgsea`, `msigdbr`,
  and `data.table`; exact versions are recorded in experiment-level logs
  (`08_logs/*_sessionInfo.txt`) where the analysis ran. The lockfile pins
  the R runtime and repository; package-level pins beyond the bootstrap
  environment should be regenerated before public release (see
  `RESOURCE_PROVENANCE.md`, item "full environment snapshot").

## Reproducibility

Scripts locate the repository root at runtime instead of using hard-coded
absolute paths. A small self-contained detector is embedded in each affected
script:

```r
find_escc_project_root <- function() {
  .escc_args_ <- commandArgs(trailingOnly = FALSE)
  .escc_file_ <- sub("^--file=", "", .escc_args_[grep("^--file=", .escc_args_)])
  .escc_start_ <- if (length(.escc_file_) && nzchar(.escc_file_)) normalizePath(.escc_file_)
                  else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) normalizePath(rstudioapi::getSourceEditorContext()$path)
                  else normalizePath(getwd())
  .p_ <- .escc_start_
  repeat {
    .d_ <- dirname(.p_)
    if (dir.exists(file.path(.d_, "00_metadata")) && dir.exists(file.path(.d_, "07_scripts")) && dir.exists(file.path(.d_, "06_tables"))) return(.d_)
    if (.d_ == .p_) break
    .p_ <- .d_
  }
  normalizePath(getwd())
}
base_dir <- find_escc_project_root()
```

The root is the directory containing `00_metadata/`, `07_scripts/`, and
`06_tables/`; all data paths are then built relative to it via `file.path()`.
This keeps scripts runnable after clone/fork without editing paths.

**Reproducibility caveat.** Large raw datasets and selected intermediate
objects are intentionally excluded. Some scripts therefore require retrieval or
regeneration of upstream public data or intermediate objects before execution.
The repository should not claim a single-command end-to-end rebuild unless that
workflow has actually been validated.

**Reproduction scope note.** This release does *not* claim one-command
end-to-end reproduction. Raw data downloads require access to the original
repositories (network + dataset availability), large preprocessing objects are
not shipped, and several steps have frozen inputs that must be produced in
order. What is reproducible from this repository with the shipped inputs:

- re-running any Step's analysis on the shipped `06_tables`/`00_metadata`
  inputs (dependency order documented per script header);
- regenerating derived tables and figures recorded as script outputs;
- reproducing the frozen QC/evidence tables that were used verbatim by later
  steps.

See `CODE_AVAILABILITY.md` and `DATA_AVAILABILITY.md` for boundaries.

## Analysis workflow / script organization

The 191 workflow scripts (143 R, 36 Python, 12 shell) in `07_scripts/` are a
mix of **canonical analysis scripts**, **supporting audits**, **pilot/sensitivity
analyses**, **freeze/reconciliation steps**, and **publication/manuscript QA
utilities**. Because many scripts depend on upstream public-data retrieval or
intermediate objects that are intentionally excluded from this repository, the
collection should *not* be assumed to be executable independently, one script
per fresh clone; dependency order is documented in each script header.

See [SCRIPT_MAP.md](SCRIPT_MAP.md) for the detailed script inventory and interpretation.

| Stage | Scope |
|-------|-------|
| 00–03 | Project setup, GEO metadata audit, download planning, preflight. |
| 04–09 | GSE221561 acquisition, extraction/audit, Seurat construction, QC, metadata standardization, mapping validation. |
| 10–11 | GSE197677 and GSE160269 acquisition/download helpers. |
| 12–13 | GSE160269 object construction, corrected metadata, fibroblast/reference annotation. |
| 14–15 | GSE197677 lineage QC, cell-type annotation, tumor eligibility correction. |
| 16 | Unified cross-dataset lineage/eligibility tables. |
| 17 | Malignant-state/CNV analyses and supporting CopyKAT/inferCNV audits. |
| 18 | Cross-dataset consistency and frozen CNV/reference layers. |
| 19 | Frozen cross-dataset analytical core: pseudobulk inventory, DGE/QC, GSE197677 DE and sensitivity, replication/orientation correction, hallmark/pathway robustness, epithelial localization/mixed-program audits, fibroblast core/abundance/heterogeneity/regulator audits, targeted ligand-receptor signaling audits, and final signaling synthesis. |
| 20 | Global project status, evidence hierarchy, and optional branch tracking. |
| 21 | Frozen result-report/manuscript/figure/table/export assembly. |
| 22A | External validation, including frozen score source, 27-sample tumor QC, histology gate, technical CORE2/CORE4 replication, and final governance/conclusion. |
| 23 | Manuscript integration, evidence-map freeze, QA edits, verified source-reference insertion, and human-review gates. |

Naming convention: `NN_STEP…_<DATASET>_<PURPOSE>`; `*_MAIN`/`*_FIX`/`*_REV*`
mark corrected or revision-level variants. Frozen inputs/outputs are marked
`STEP*_FROZEN*`, `*_FREEZE*`, `*_FINAL*`, or `*_LOCK*` and are consumed verbatim
by downstream stages.

## Figures and tables

- **Figures (`05_figures/`)** — 70 files: publication figure builds (Step21D
  manuscript finals: Figures 1-6, Supplementary Figures 1-5), per-dataset
  lineage/QC heatmaps and UMAPs, pathway and regulator plots, and Step22A
  external-validation figures.
- **Tables (`06_tables/`)** — 351 files: frozen evidence tables
  (e.g., `STEP19D_*_FREEZE*/FINAL*`, `STEP20_FINAL_EVIDENCE_HIERARCHY.csv`),
  eligible/final LR-pair lists, hallmark pathway results, and per-dataset
  QC/audit tables consumed by later steps.
- **Metadata (`00_metadata/`)** — 161 files: GEO series audit, download
  manifests, sample metadata, OMIX source tables, and provenance / keyword
  audit records used by the scoring source (Step 22A).

The 11 manuscript-referenced Step21D figures are tracked at their canonical
paths; superseded FINAL/REV variants and pre-revision `main/` copies remain on
disk locally but are excluded from the index (see `.gitignore`, D4b).

## Data availability

See `DATA_AVAILABILITY.md`.

## Code availability

See `CODE_AVAILABILITY.md`.

## Citation

This repository is citable via `CITATION.cff` (machine-readable) and the
metadata it provides. `version`, DOIs, and the final release repository URL
are placeholders pending public release — update them before citing.

## License

License selection is pending before public release. The scope of the
repository contents (analysis code + selected derived outputs) and the
boundary around third-party data rights are described in `LICENSE_SCOPE.md`.

## Contact

[Corresponding author / maintainer contact to be added after publication —
replace this line before release.]