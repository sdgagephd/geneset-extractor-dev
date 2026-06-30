# CPTAC/PDC Gene-Set Pipeline Proposal

Full design spec: `docs/superpowers/specs/2026-06-29-cptac-pdc-pipeline-design.md`

## Overview

This pipeline extracts phosphoregulation gene sets from CPTAC discovery cohorts processed by the CDAP (Clinical Data Analysis Pipeline) and publicly hosted on the NCI Proteomic Data Commons (PDC). For each cohort it runs a tumor-vs-adjacent-normal comparison at the phosphosite level, protein-adjusted where proteome data are available, and emits signed gene sets in the format consumed by the DIG portal.

Phase 1 implements a single cohort (`ccrcc`) and a single model (`PT1`). Phase 2 extends the runner to batch over all 10 registered discovery cohorts. Phase 3 adds cluster submit scripts for production-scale runs.

## Four-stage pipeline flow

### 1. Fetch

`fetch_pdc_study.run_fetch` issues GraphQL queries against the PDC public API for each cohort's phospho and proteome study IDs. It selects the canonical CDAP phosphosite and proteome report TSVs from the file listing, downloads them via PDC signed URLs, and cross-checks each file's MD5 against the value returned by the API. It derives two additional files: `sample_annotations.tsv` (aliquot → condition mapping from biospecimen `sample_type`) and `pdc_file_manifest.tsv` (provenance metadata including DRS URIs).

### 2. Prepare

The dig engine's `ptm_prepare_public` workflow ingests the CDAP phosphosite and proteome TSVs plus `sample_annotations.tsv`. It normalises the CDAP `"<aliquot_submitter_id> Log Ratio"` column header format to bare aliquot IDs before joining sample annotations, and accepts lowercase phospho residue characters (`s`, `t`, `y`) from CDAP notation. Output is a standardised PTM matrix, protein matrix, and sample metadata file.

### 3. Extract

The dig engine's `ptm_site_matrix` converter runs the tumor-vs-normal contrast on the prepared matrices. For model `PT1` with `--protein_adjustment_run_mode compare_if_protein`, the extractor emits both a protein-adjusted variant and an unadjusted variant when a protein matrix is present. Gene-level aggregation uses `signed_topk_mean` over the top 200 phosphosites by effect size.

### 4. Provenance overlay

`build_pdc_provenance_overlay.write_overlay` constructs a `provenance_overlay.json` consumed by the extractor. The overlay supplies CRDC-native provenance fields for each prepared input file (phospho matrix, protein matrix, sample metadata). See the CRDC-vs-CFDE section below for the specific values and the rationale for the choices made.

## Model PT1

`PT1` is the Phase-1 tumor-vs-adjacent-normal phosphoregulation model:

- Study contrast: `case` (Primary Tumor) vs `control` (Solid Tissue Normal)
- Enabled for cohorts where `has_adjacent_normal=true` in the study manifest
- Produces two extractor output variants:
  - `protein_adjusted` — phosphosite log-ratio corrected for matched protein abundance change
  - `unadjusted` — raw phosphosite log-ratio without protein correction
- Both variants are emitted by `ptm_site_matrix` when `--protein_adjustment_run_mode compare_if_protein` is set and a protein matrix is available

## CRDC-vs-CFDE provenance decision

PDC belongs to the **NCI Cancer Research Data Commons (CRDC)**. It is not part of the NIH Common Fund Data Ecosystem (CFDE), and PDC studies do not have entries in the CFDE Workbench. Attempting to construct a genuine CFDE `drc_url` for PDC data would be incorrect.

The repository's provenance schema requires `dcc_url` and `drc_url` fields on every node. For CPTAC data the pipeline populates these fields with CRDC-native values by explicit decision:

| Field | Value | Notes |
|---|---|---|
| `persistent_id` / `local_id` | `drs://dg.4DFC/<file_uuid>` | CRDC GA4GH DRS URI; `file_uuid` is the PDC `file_id` returned by `filesPerStudy` |
| `dcc_url` | `https://pdc.cancer.gov/pdc/study/<pdc_study_id>` | Direct link to the PDC study page for the source file |
| `drc_url` | `https://datacommons.cancer.gov/repository/proteomic-data-commons` | CRDC repository landing page; used as a CRDC stand-in for the `drc_url` slot since PDC has no CFDE presence |
| `md5` / `size` | real values | Engine-computed at fetch time; cross-checked against the PDC API `md5sum` field |

This approach keeps provenance honest (no invented CFDE identifiers) while satisfying the schema's structural requirements. The `drc_url` constant (`CRDC_DRC_URL`) is defined once in `build_pdc_provenance_overlay.py` and applied uniformly to all input and operation nodes.

## Phase 2 and Phase 3 (future work)

**Phase 2** will extend `build_cptac_genesets.sh` to loop over all cohorts in `study_manifest.tsv` with `enabled=true`, running the per-study `run_cptac_ptm_model.py` runner for each cohort × model combination.

**Phase 3** will add cluster submit scripts (e.g. UGER/SLURM job arrays) to parallelize per-cohort runs in production, matching the pattern used by other pipelines in this repository.
