# Resource Provenance

> Template for the ligand–receptor (LR) resource audit and, more broadly, for
> every external resource referenced by the analysis pipeline.
>
> Rule: **no value is filled from memory.** Only values that are verifiable in
> this repository (scripts, tables, metadata, lockfile) are entered; everything
> else is `TO BE VERIFIED BEFORE PUBLIC RELEASE`.

Legend: **Verified in repo** = value is directly traceable to a committed file
in this repository. **TO BE VERIFIED** = value must be confirmed against the
original resource before public release.

---

## Provenance record template

| Field | Value / instruction |
|-------|---------------------|
| **Resource** | Name of the resource (database/package/dataset). |
| **Purpose** | What the resource does in this analysis (input? reference? validation? derived-source?). |
| **Exact object/file name** | Exact filename/path as used in scripts/tables. |
| **Package/database** | The package or database that provides it (if part of one). |
| **Version** | Version identifier. Only record if present in a committed file. |
| **Freeze/retrieval date** | Date the resource state was frozen/retrieved. Only record if present in a committed file. |
| **Source URL** | Original distribution URL. Only record if present in a committed file. |
| **DOI** | Digital Object Identifier of the resource/release. Only record if present in a committed file. |
| **Checksum** | Content checksum (e.g., MD5/SHA). Only record if present in a committed file; otherwise TO BE VERIFIED. |
| **Used by script(s)** | Exact script path(s) in `07_scripts/` that reference it. |
| **Verification status** | `Verified in repo` / `TO BE VERIFIED` (+ what to do at release). |
| **Notes** | Additional constraints, caveats, or links to status tables. |

---

## Resource records (current, evidence-based)

### 1. curated_minimal LR resource (primary LR pair set)

| Field | Value |
|-------|-------|
| **Resource** | `curated_minimal` — project-local, script-embedded minimal ligand–receptor pair set |
| **Purpose** | Source of LR pairs for Steps 19O–19P targeted fibroblast signaling audit |
| **Exact object/file name** | Defined in-script (`LR resource used: curated_minimal`); outputs: `06_tables/cross_dataset/pathway_analysis/Step19O/STEP19O_ELIGIBLE_LIGAND_RECEPTOR_PAIRS.csv`, `.../STEP19O_FINAL_LIGAND_RECEPTOR_PAIRS.csv`, `.../Step19P/STEP19P_LIGAND_TARGET_RESOURCE_STATUS.csv` |
| **Package/database** | None (project-local definition; no external LR database source identified in repository history) |
| **Version** | Project-local script definition (repository/status-table label: `prespecified_in_script`) |
| **Freeze/retrieval date** | No external retrieval date; earliest repository evidence is root commit `45fab97` (2026-08-20) |
| **Source URL** | None (project-local definition) |
| **DOI** | None identified in repository history |
| **Checksum** | n/a (script-derived); pair tables are committed as CSVs |
| **Used by script(s)** | `07_scripts/19O_CORR_CLASSIFICATION_RECONCILIATION.R`, `07_scripts/19O_FIX221_GSE221561_SIGNALING_COMPLETION.R`, `07_scripts/19O_TARGETED_FIBROBLAST_SIGNALING_AUDIT.R`, `07_scripts/19P_FOCUSED_REPLICATED_SIGNAL_VALIDATION.R` |
| **Verification status** | Exact pair definitions and downstream use verified in repository |
| **Content-level selection provenance** | Not recoverable from repository history; the available history does not document how the original 100 pairs were selected |
| **Notes** | `STEP19O_RESOURCE_STATUS.csv`: `available=TRUE`, `n_pairs=100`. `STEP19O_CORR_RESOURCE_LIMITATIONS.csv` records: 100 LR pairs, `MINIMAL_CURATED_RESOURCE`, no ligand-target matrix. Eligible/final pair tables contain 64 rows with families CXCL_CCL, EGF, NOTCH, TGF_beta, TNF (verified in the CSVs). |

### 2. CellChatDB

| Field | Value |
|-------|-------|
| **Resource** | CellChatDB |
| **Purpose** | External LR database considered for cross-checking; **not available** at analysis time |
| **Exact object/file name** | n/a (marked unavailable) |
| **Package/database** | CellChatDB (external; not retrieved) |
| **Version** | N/A — resource not retrieved or used; repository status tables record `NA` |
| **Freeze/retrieval date** | N/A — resource not retrieved |
| **Source URL** | N/A — no external resource was retrieved or used |
| **DOI** | N/A — no external resource was retrieved or used |
| **Checksum** | N/A — no local resource file |
| **Used by script(s)** | `07_scripts/19O_TARGETED_FIBROBLAST_SIGNALING_AUDIT.R` (availability preflight) |
| **Verification status** | Verified from repository evidence as unavailable at analysis time; no external resource version/retrieval metadata applies |
| **Notes** | `STEP19O_RESOURCE_STATUS.csv`: `available=FALSE, n_pairs=0`. `STEP19O_CORR_RESOURCE_LIMITATIONS.csv`: `external_LR_databases NOT_AVAILABLE (CellChatDB/NicheNet/OmniPath/LIANA)`. |

### 3. NicheNet resources

| Field | Value |
|-------|-------|
| **Resource** | `NicheNet_LR` / `NicheNet_LR_network` / `NicheNet_ligand_target` |
| **Purpose** | LR network and ligand–target validation; **not available** at analysis time |
| **Exact object/file name** | n/a (marked unavailable) |
| **Package/database** | NicheNet resources (external; not retrieved) |
| **Version** | N/A — resource not retrieved or used; repository status tables record `NA` |
| **Freeze/retrieval date** | N/A — resource not retrieved |
| **Source URL** | N/A — no external resource was retrieved or used |
| **DOI** | N/A — no external resource was retrieved or used |
| **Checksum** | N/A — no local resource file |
| **Used by script(s)** | `07_scripts/19P_FOCUSED_REPLICATED_SIGNAL_VALIDATION.R` (ligand-target preflight), `07_scripts/19O_TARGETED_FIBROBLAST_SIGNALING_AUDIT.R` |
| **Verification status** | Verified from repository evidence as unavailable at analysis time; no external resource version/retrieval metadata applies |
| **Notes** | `STEP19P_LIGAND_TARGET_RESOURCE_STATUS.csv`: NicheNet_ligand_target / NicheNet_LR_network / project_local_LT / other_local_LT all `available=FALSE`, `version=NA`. No ligand-target matrix used anywhere in the pipeline. |

### 4. OmniPath & LIANA

| Field | Value |
|-------|-------|
| **Resource** | OmniPath; LIANA |
| **Purpose** | External LR/interaction resources considered; **not available** at analysis time |
| **Exact object/file name** | n/a (marked unavailable) |
| **Package/database** | OmniPath / LIANA (external; not retrieved) |
| **Version** | N/A — resource not retrieved or used; repository status tables record `NA` |
| **Freeze/retrieval date** | N/A — resource not retrieved |
| **Source URL** | N/A — no external resource was retrieved or used |
| **DOI** | N/A — no external resource was retrieved or used |
| **Checksum** | N/A — no local resource file |
| **Used by script(s)** | `07_scripts/19O_TARGETED_FIBROBLAST_SIGNALING_AUDIT.R` (availability preflight) |
| **Verification status** | Verified from repository evidence as unavailable at analysis time; no external resource version/retrieval metadata applies |
| **Notes** | `STEP19O_RESOURCE_STATUS.csv`: `available=FALSE`. |

### 5. R environment / lockfile (reproducibility resource)

| Field | Value |
|-------|-------|
| **Resource** | renv lockfile + R runtime |
| **Purpose** | Anonymous/safe reproduction of the analysis environment |
| **Exact object/file name** | `renv.lock`, `renv/activate.R`, `renv/settings.json`, `.Rprofile` |
| **Package/database** | renv; R; Bioconductor |
| **Version** | R `4.6.1`; Bioconductor `3.23`; renv `1.2.4`; BiocManager `1.30.27` (all verified in `renv.lock`) |
| **Freeze/retrieval date** | TO BE VERIFIED BEFORE PUBLIC RELEASE (recorded per-analysis in `08_logs/*_sessionInfo.txt`) |
| **Source URL** | CRAN repo `https://cloud.r-project.org` (in `renv.lock`) |
| **DOI** | None |
| **Checksum** | TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Used by script(s)** | All R scripts (environment bootstrap) |
| **Verification status** | Verified in repo (lockfile content); full package-level pinning beyond bootstrap = TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Notes** | Starting point: regeneration should run `renv::restore()` plus mandated package installs; analysis packages (e.g., edgeR, Seurat, fgsea, msigdbr) must be pinned by a full snapshot before release. |

### 6. Public datasets (provenance pointer)

| Field | Value |
|-------|-------|
| **Resource** | GSE221561, GSE197677, GSE160269, TCGA-ESCA, GSE53625, DepMap |
| **Purpose** | Input data / validation context / optional branches (see README) |
| **Exact object/file name** | Download targets/manifests in `00_metadata/` (`*_download_targets.tsv`, `_supplementary_files.csv`, etc.) |
| **Package/database** | Original dataset repositories (unchanged terms) |
| **Version** | Accession-level only; exact archived versions TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Freeze/retrieval date** | TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Source URL** | TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **DOI** | TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Checksum** | TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Used by script(s)** | Per-dataset pipeline scripts (Steps 01–22A; see `README.md` workflow table) |
| **Verification status** | Accessions verified in repo metadata; URLs/DOIs/dates/checksums TO BE VERIFIED BEFORE PUBLIC RELEASE |
| **Notes** | Raw data are not redistributed; see `DATA_AVAILABILITY.md`. TCGA-ESCA / GSE53625 / DepMap are referenced in Step 20/22A planning/audit records (`STEP20_OPTIONAL_NEXT_BRANCHES.csv`, `STEP22A_04_KEYWORD_HITS.csv`). |

---

## STEP24N LR resource provenance freeze (2026-08-23)

Provenance reconciliation for the targeted signaling ligand–receptor (LR)
resource. Git archaeology of the repository is complete; there is no
pre-repository origin from which an earlier LR-resource provenance can be
reconstructed (root/initial commit `45fab97` is the first repository commit).

### Freeze record

| Field | Value |
|-------|-------|
| STEP24N_LR_RESOURCE_PROVENANCE | PASS |
| RESOURCE_NAME | curated_minimal |
| RESOURCE_ORIGIN | project-local |
| RESOURCE_IMPLEMENTATION | script-embedded |
| RESOURCE_STATUS | prespecified minimal panel |
| LR_PAIRS | 100 |
| SIGNALING_FAMILIES | 13 |
| FIRST_REPOSITORY_COMMIT | 45fab97 |
| FIRST_REPOSITORY_COMMIT_IS_ROOT | YES |
| EXTERNAL_LR_DATABASE_USED | NO |
| PAIR_LIST_REPRODUCIBLE_FROM_SCRIPT | YES |
| PAIR_LEVEL_SOURCE_CITATIONS_DOCUMENTED | NO |
| PAIR_SELECTION_ALGORITHM_DOCUMENTED | NO |
| EXTERNAL_DATABASE_ATTRIBUTION_ALLOWED | NO |

### Resolved vs unresolved provenance

- **(a) Computational / resource provenance — RESOLVED.** The 100-pair panel
  is defined in-script as `curated_minimal`, is project-local and
  script-embedded, is reproducible from the committed scripts, and was the
  only LR interaction resource used. At the root commit `45fab97` the pipeline
  explicitly checked CellChatDB, NicheNet_LR, NicheNet_ligand_target, OmniPath,
  LIANA, and project_local_LR; the external resources were unavailable and were
  not used. Commit `aeeb6ea` later clarified the project-local provenance but
  did not establish any external database origin.
- **(b) Pre-repository pair-selection provenance — NOT DOCUMENTED.** The
  repository does not contain a pair-level literature citation map, an external
  database accession/version for the 100 pairs, or a reproducible algorithm
  explaining exactly why these 100 pairs were selected. This is a gap in
  content-level selection provenance, **not** evidence that the interactions
  themselves lack biological support.

### Attribution guardrails (do not publish)

- Do **not** attribute the 100-pair panel to CellChatDB, NicheNet, OmniPath,
  LIANA, CellPhoneDB, or any other external comprehensive LR database.
- Do **not** invent pair-level literature provenance.
- Do **not** describe the absence of documented provenance as evidence that
  the interactions lack biological support.

### Manuscript-facing consistency status

- Manuscript drafts (README; Step19O/19P interpretation + guardrails;
  Step20 safe-language; Step21A/B limitations + abstract) consistently
  describe the resource as a project-local, script-embedded, minimal curated
  LR panel (~100 pairs / 13 families), explicitly note that external LR
  databases were unavailable, and state the panel is not comprehensive.
- No unsupported external-database attribution was found in the reviewed
  manuscript-facing files. No manuscript edit is required at this time.

---

## STEP24P Step23C-QA2B1C bibliography hold closure (2026-08-23)

Closure of the Step23C-QA2B1C bibliographic-verification HOLD
(QA2C-054). The three dataset-source publication links were
independently verified against authoritative primary accession /
publication records (STEP24O dataset bibliography audit = PASS).

### Closure record

| Field | Value |
|-------|-------|
| STEP23C_QA2B1C_STATUS | CLOSED |
| GSE197677 | VERIFIED_PRIMARY_PUBLICATION — PMID 37091252, DOI 10.1016/j.isci.2023.106480 |
| GSE221561 | VERIFIED_PRIMARY_PUBLICATION — PMID 37563120, DOI 10.1038/s41392-023-01518-0 |
| OMIX005710 | VERIFIED_PRIMARY_PUBLICATION — PMID 38566201, DOI 10.1186/s13073-024-01320-9 |
| ALL_THREE_ACCESSION_IDENTITIES_CONFIRMED | YES |
| ALL_THREE_BIBLIOGRAPHIC_RECORDS_DEFENSIBLE | YES |
| BIBLIOGRAPHIC_CORRECTION_REQUIRED | NO |
| MANUSCRIPT_CORRECTION_REQUIRED | NO |
| SCIENTIFIC_RESULTS_CHANGED | NO |
| SCIENTIFIC_RERUN_REQUIRED | NO |

Note: the canonical working hold record
`04_results/manuscript_update/Step23/Step23C_QA2B1C_REFERENCE_VERIFICATION_HOLD.txt`
resides under a gitignored working directory
(`/04_results/manuscript_update/`); its local copy carries the same
closure text. This tracked provenance record is the committed source of
truth for the closure.

---

## Open verification tasks (do before public release)

1. Resolve exact version/URL/DOI/retrieval date/checksum for any external
   resource the final manuscript relies on (CellChatDB, NicheNet, OmniPath,
   LIANA) — or confirm the final manuscript makes **no** such claim.
2. If an external record of the original curation source for the `curated_minimal` 100 LR pairs becomes available outside this repository, add it before release; otherwise state explicitly that the exact project-local pair set is reproducible from code but the original content-level selection provenance is not recoverable from repository history.
3. Generate a full `renv.lock` snapshot pinning all analysis packages, or
   record exact package versions from `08_logs/*_sessionInfo.txt` per release.
4. Freeze dataset retrieval dates + checksums for the datasets actually used.
5. After publication (DOI/release), back-fill the archive DOI, GitHub URL, and
   release date into this file, `CITATION.cff`, `CODE_AVAILABILITY.md`, and
   `DATA_AVAILABILITY.md`.