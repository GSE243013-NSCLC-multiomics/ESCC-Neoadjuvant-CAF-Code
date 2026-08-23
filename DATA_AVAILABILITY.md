# Data Availability

## Public data sources

This study used and referenced the following public datasets. All raw source
data remain available from their original repositories/providers and are
subject to their original terms of use:

| Dataset        | Original provider/repository | Availability of raw data            |
|----------------|------------------------------|-------------------------------------|
| GSE221561      | Public genomics repository (series accession as recorded in `00_metadata/`) | From original repository (unchanged terms) |
| GSE197677      | Public genomics repository (series accession as recorded in `00_metadata/`) | From original repository (unchanged terms) |
| GSE160269      | Public genomics repository (series accession as recorded in `00_metadata/`) | From original repository (unchanged terms) |
| TCGA-ESCA      | TCGA program repository       | From original TCGA/GDC-controlled access (subject to TCGA data-use terms) |
| GSE53625       | Public genomics repository (series accession referenced in Step 20 planning) | From original repository (unchanged terms) |
| DepMap         | Project DepMap (Broad Institute) | From Project DepMap portal (subject to its terms) |

Exact download URLs, retrieval dates, checksums, and any accession-level
details are documented in `RESOURCE_PROVENANCE.md`; items not yet verified are
explicitly marked `TO BE VERIFIED BEFORE PUBLIC RELEASE`.

## What this repository does and does not include

- **Not redistributed:** large raw datasets (e.g., FASTQ, raw UMI count
  matrices, and other source-level files) are **not** redistributed here.
  They remain available from their original repositories/providers.
- **Included for transparency:** selected **derived tables and metadata**
  necessary to reproduce the reported analyses and to audit provenance
  (e.g., sample metadata, pseudobulk counts used at Steps 19A+,
  eligible/final ligand-receptor pair lists, freeze/final evidence tables,
  step-status and evidence-hierarchy tables).
- **Excluded regenerable intermediates:** large intermediate analysis objects
  (Seurat objects, big matrix objects, intermediate R binaries) are excluded.
  Where applicable, they can be regenerated from the source data using the
  provided scripts in dependency order.

## Controlled accession

No controlled-access statement is made in this release. Access restrictions
(if any) that apply to specific datasets derive from the datasets' own
original providers and are not additional to, or different from, those terms.
This statement will be updated to reference any finalized accession-specific
data-use agreements before publication if needed.

---

*Status: pre-publication draft. URLs, DOIs, retrieval dates, and accession
details to be finalized before public release.*