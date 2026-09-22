# Data Availability

## Principal public datasets used in the final Cancers manuscript

The final manuscript reanalyzes three public single-cell RNA-sequencing datasets. Raw source data are not redistributed in this repository.

| Dataset | Repository | Public access | Final manuscript role |
|---|---|---|---|
| **GSE197677** | NCBI Gene Expression Omnibus (GEO) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197677 | Cross-sectional discovery cohort |
| **GSE221561** | NCBI Gene Expression Omnibus (GEO) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221561 | Cross-sectional discovery cohort |
| **OMIX005710** | National Genomics Data Center OMIX | https://ngdc.cncb.ac.cn/omix/release/OMIX005710 | Independent paired sensitivity analysis |

Primary source publications verified in the project provenance record are:

- GSE197677 — PMID 37091252; DOI `10.1016/j.isci.2023.106480`
- GSE221561 — PMID 37563120; DOI `10.1038/s41392-023-01518-0`
- OMIX005710 — PMID 38566201; DOI `10.1186/s13073-024-01320-9`

## Derived analytical outputs

This study generates new **derived analytical results** from the public source data, including sample/patient-level pseudobulk summaries, fixed CORE2/CORE4 score analyses, histology-sensitivity summaries, influence audits, figures, and numerical result tables.

The final manuscript's analysis-specific eligibility, provenance, sensitivity outputs, and numerical audits are supplied with the journal submission as **Supplementary Files S1 and S2**.

The public GitHub repository is primarily a code-and-documentation release and should not be interpreted as a complete mirror of every working-directory result or intermediate generated during project development.

## Historical/developmental resources

Historical scripts may refer to additional resources such as GSE160269, TCGA-ESCA, GSE53625, and DepMap. These records are retained for workflow provenance. They do **not** contribute to the final Cancers manuscript's principal three-cohort inference unless explicitly stated in the manuscript.

## Redistribution and access

- Large raw source datasets are **not redistributed** here.
- Large Seurat objects, matrix objects, and other heavy intermediates are not redistributed.
- No new patient recruitment, specimen collection, or intervention was performed for this study.
- Access conditions for original public datasets remain those of their source repositories.

## Code archive

Analysis code and workflow documentation:

- GitHub: https://github.com/GSE243013-NSCLC-multiomics/ESCC-Neoadjuvant-CAF-Code
- Zenodo software archive (all versions): https://doi.org/10.5281/zenodo.22083689


For the associated manuscript, the version-specific DOI of the final archived release should be cited to identify the exact software snapshot used.
