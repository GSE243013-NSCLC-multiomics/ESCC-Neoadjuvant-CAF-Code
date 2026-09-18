# ESCC Neoadjuvant CAF Score Transportability Analysis

## Associated manuscript

**Treatment-associated direction of a fixed fibroblast CDKN1B-GSN score varies across sampling contexts in neoadjuvant-treated esophageal squamous cell carcinoma**

**Authors:** Xiaowei Xie, Ying Guo, Lanmei Yu  
**Affiliation:** Department of Cardiothoracic Surgery, The First Hospital of Putian, Putian, Fujian, China  
**Corresponding author:** Xiaowei Xie — xiexiaoweivip@outlook.com

## Overview

This repository contains analysis code, workflow documentation, environment metadata, and provenance materials supporting a public-data reanalysis of fibroblast programs in esophageal squamous cell carcinoma (ESCC) treated in neoadjuvant settings.

The final Cancers manuscript focuses on a fixed fibroblast score:

- **CORE2 = -z(CDKN1B) - z(GSN)**
- **CORE4 = CORE2 + z(BGN) + z(TIMP1)**, used as a sensitivity comparator.

The genes, weights, and orientation of the score were locked before inspection of the paired OMIX005710 results. “Locked” does not mean that CORE2 was prospectively prespecified before the original discovery workflow.

The principal manuscript finding is intentionally descriptive: CORE2 had a positive treated-minus-control point-estimate direction in two cross-sectional discovery cohorts, whereas post-treatment minus pretreatment CORE2 change was negative across all prespecified OMIX005710 histology-sensitivity settings. Because cohort, treatment regimen, specimen route, platform, patient matching, and biological time window differ together, the study does **not** attribute this discordance to any single factor and does **not** present CORE2 as a validated predictive, pharmacodynamic, prognostic, or treatment-selection biomarker.

## Principal manuscript datasets

| Dataset | Final manuscript role | Final analytic contrast |
|---|---|---|
| **GSE197677** | Cross-sectional discovery cohort | 6 neoadjuvant chemotherapy-treated vs. 4 untreated ESCC tumors |
| **GSE221561** | Cross-sectional discovery cohort | 7 post-neoadjuvant residual vs. 2 surgery-alone ESCC tumors |
| **OMIX005710** | Independent paired sensitivity analysis | 26 technically eligible patient-timepoints from 21 patients; 4-5 complete pre/post pairs per histology-sensitivity setting; post-treatment minus pretreatment |

These cohorts estimate clinically non-equivalent contrasts. They are therefore compared as a **transportability stress test**, not pooled as if they estimated the same treatment effect.

### Historical/developmental resources

Historical scripts and provenance records may refer to additional resources considered or used during earlier project-development branches, including **GSE160269, TCGA-ESCA, GSE53625, and DepMap**. Their presence in historical code does **not** mean that they contribute to the final Cancers manuscript's principal three-cohort inference. Such branches are retained for provenance rather than presented as additional independent validation.

## Analysis principles

- The sample or patient-timepoint, not the individual cell, is the principal inferential unit.
- Fibroblast expression is summarized by sample/patient-level pseudobulk.
- CORE2 and CORE4 are discovery-derived, post-selection constructs.
- GSE221561 is directionally supportive but imprecise because the surgery-alone control group contains only two samples.
- OMIX005710 histology-sensitivity settings reuse the same four to five paired patients and are **not** independent replications.
- Leave-one-patient-out calculations are influence analyses and do not increase effective sample size.
- Regulator and source-ligand analyses are hypothesis-generating and do not establish a confirmed mechanism.
- No meta-analysis is used to combine clinically non-equivalent cross-sectional and paired contrasts.

## Public data availability

Raw public data are **not redistributed** in this repository.

- **GSE197677 — NCBI Gene Expression Omnibus**  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197677
- **GSE221561 — NCBI Gene Expression Omnibus**  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221561
- **OMIX005710 — National Genomics Data Center OMIX**  
  https://ngdc.cncb.ac.cn/omix/release/OMIX005710

Analysis-specific eligibility, provenance, sensitivity outputs, and numerical audits used in the final manuscript are supplied with the journal submission as **Supplementary Files S1 and S2**.

See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for details.

## Public repository scope

The public GitHub repository is primarily a **code-and-documentation release**. The current public root contains:

```text
07_scripts/              workflow scripts and analysis/manuscript-QA utilities
renv/                    renv project scaffolding
.gitignore
CITATION.cff
CODE_AVAILABILITY.md
DATA_AVAILABILITY.md
LICENSE
README.md
RESOURCE_PROVENANCE.md
SCRIPT_MAP.md
renv.lock
```

Large raw datasets, large processed objects, and heavy computational intermediates are intentionally not redistributed. Historical scripts may reference working directories such as `00_metadata/`, `04_results/`, `05_figures/`, `06_tables/`, and `08_logs/`; those paths document the original analysis workflow and should not be interpreted as files currently shipped in the public code-only repository.

The public repository does **not** claim a single-command end-to-end rebuild from a fresh clone. Upstream public inputs and some intermediate objects must first be retrieved or regenerated according to the relevant scripts.

## Software environment

The uploaded/public `renv.lock` records:

- **R 4.6.1**
- **Bioconductor 3.23**
- **192 package records**
- `renv` 1.2.4
- `BiocManager` 1.30.27
- `Seurat` 5.5.1
- `edgeR` 4.10.1
- `fgsea` 1.38.0
- `msigdbr` 26.1.0
- `data.table` 1.18.4
- CRAN repository: `https://cloud.r-project.org`

The lockfile therefore contains package-level pins for the released R environment. Reproduction remains subject to platform compatibility and continued availability of the recorded package sources.

SHA-256 of the reviewed `renv.lock`:

```text
b52009f6d6efadadfa75659189974153f0dbdf90cb3995c13e1765abd9ad8275
```

## Script organization

The scripts under `07_scripts/` retain their historical stage numbering and revision labels for provenance. Labels such as `MAIN`, `FIX`, `REV`, `FREEZE`, `FINAL`, and `LOCK` describe workflow state; they are **not** separate biological claims.

The high-level stages are:

| Stage | Scope |
|---|---|
| 00-03 | Project setup, GEO metadata audit, download planning, preflight |
| 04-09 | GSE221561 acquisition, QC, Seurat construction, metadata standardization |
| 10-15 | GSE197677/GSE160269 acquisition and lineage/reference-development branches |
| 16-18 | Cross-dataset lineage/CNV/reference audits |
| 19 | Frozen cross-dataset analytical core and fibroblast/regulator/signaling audits |
| 20 | Project status, evidence hierarchy, optional-branch tracking |
| 21 | Result/manuscript/figure/table assembly |
| 22A | OMIX005710 paired sensitivity analysis, technical eligibility, histology-sensitivity and influence audits |
| 23 | Manuscript integration, evidence-map freeze, QA and verified-reference tooling |

See [SCRIPT_MAP.md](SCRIPT_MAP.md) for the detailed inventory.

## Reproducibility boundary

The repository supports auditability of the published computational workflow, but it does not imply that every historical script is an independent final analysis or that every intermediate can be regenerated without first restoring upstream inputs.

In particular:

- historical/pilot/fix scripts are retained intentionally;
- large raw and intermediate data are not redistributed;
- final numerical manuscript claims should be checked against the journal supplementary files and the frozen evidence tables used by the final workflow;
- historical script names such as “validation” or “replication” are workflow labels and should not be read as stronger manuscript claims.

## Code availability

GitHub repository:  
https://github.com/GSE243013-NSCLC-multiomics/ESCC-Neoadjuvant-CAF-Code

Zenodo archived release **v1.0.0**:  
https://doi.org/10.5281/zenodo.22083690

The GitHub hosting-organization name is a legacy account label and is unrelated to the disease analyzed in this ESCC project.

See [CODE_AVAILABILITY.md](CODE_AVAILABILITY.md) for details.

## Evidence boundary

The current evidence does **not** establish:

- a clinically validated CORE2 biomarker;
- treatment-selection utility;
- a validated predictive or pharmacodynamic test;
- a confirmed master regulator or signaling mechanism;
- causality for the observed cross-context direction difference.

The translational message is that stromal scores should be validated within a prespecified clinical context of use with aligned specimen timing, treatment contrast, and patient-level inference.

## Citation

If using this software release, please cite:

**Xie X, Guo Y, Yu L. ESCC Neoadjuvant CAF Score Transportability Analysis. Version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.22083690**

Please also cite the associated Cancers article after publication.

Machine-readable citation metadata are provided in [CITATION.cff](CITATION.cff).

## License

Repository code is released under the **MIT License** as specified in `LICENSE`. Original public datasets remain subject to the terms of their source repositories and are not relicensed here.

## Contact

**Xiaowei Xie**  
Department of Cardiothoracic Surgery, The First Hospital of Putian, Putian, Fujian, China  
Email: xiexiaoweivip@outlook.com  
ORCID: https://orcid.org/0009-0009-4221-2777
