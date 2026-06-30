# CPTAC/PDC Pipeline Inputs

The current CPTAC pipeline supports:

- phosphosite-level tumor-vs-adjacent-normal models:
  - `PT1` — phospho tumor-vs-adjacent-normal, protein-adjusted (or unadjusted where protein data are absent)

Phase 1 enables one cohort:

- `ccrcc` — Clear Cell Renal Cell Carcinoma (CPTAC3, proteome PDC000127, phospho PDC000128)

Nine additional discovery cohorts are registered in the study manifest but set `enabled=false` pending Phase 2 batch rollout.

## Planning inputs

- `CPTAC/config/study_manifest.tsv`
- `CPTAC/config/model_list.tsv`
- `CPTAC/config/model_manifest.tsv`
- `CPTAC/config/model_description_templates.tsv`

`study_manifest.tsv` maps each `cohort_id` to a proteome/phospho `pdc_study_id` pair and records whether the cohort has adjacent normal tissue and whether it is enabled. `model_list.tsv` and `model_manifest.tsv` together supply the per-model workflow and extractor flag sets.

## Data source

All CPTAC CDAP processed matrices are open-access (CC-BY). The pipeline fetches them from the **NCI Proteomic Data Commons (PDC)** public GraphQL API:

- Endpoint: `https://pdc.cancer.gov/graphql` (HTTP POST `{"query": ...}`, no authentication required)
- `acceptDUA: true` is passed in every query to assert acceptance of the PDC data-use agreement

Three queries are issued per cohort: two study-resolution queries to obtain `study_id` UUIDs, one file-listing query per study (keyed by UUID), and one biospecimen query (keyed by `pdc_study_id`):

- `study(pdc_study_id: "<id>", acceptDUA: true)` — resolves a `pdc_study_id` (e.g. `PDC000127`) to its `study_id` UUID; the pipeline first calls this for each study because `filesPerStudy(pdc_study_id: ...)` returns null-filled records for some study versions, whereas `filesPerStudy(study_id: <UUID>)` is reliable
- `filesPerStudy(study_id: "<uuid>", acceptDUA: true)` — returns `file_id`, `file_name`, `md5sum`, `file_size`, `signedUrl { url }` for every file in the study
- `biospecimenPerStudy(pdc_study_id: "<id>", acceptDUA: true)` — returns `aliquot_submitter_id`, `sample_submitter_id`, `case_submitter_id`, `sample_type` for every aliquot (keyed by `pdc_study_id`; this query is unaffected by the null-fill issue)

For `ccrcc` the biospecimen query returns 110 Primary Tumor aliquots (mapped to condition `case`) and 84 Solid Tissue Normal aliquots (mapped to condition `control`). Aliquots with `sample_type` values that do not map to either condition (e.g. Not Reported, Cell Line) are written to `sample_annotations.tsv` with an empty condition field and are excluded from contrasts at the prepare step.

## Per-study inputs fetched automatically

The pipeline resolves and downloads two CDAP report TSVs per cohort. For `ccrcc` the concrete file names are:

- phosphosite report: `CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv`
- proteome report: `CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv`

The phosphosite report header is `Phosphosite`, followed by approximately 207 sample columns formatted as `"<aliquot_submitter_id> Log Ratio"`, then `Peptide`, `Gene`, `Organism`. Phosphosite identifiers use lowercase residue letters and versioned RefSeq accessions, e.g. `NP_000005.2:y708`.

Two derived files are written under the cohort fetch directory:

- `sample_annotations.tsv` — one row per aliquot; columns include `sample_id_raw` (aliquot submitter ID), `condition` (`case` / `control` / empty), `case_submitter_id`, `sample_submitter_id`, `sample_type`
- `pdc_file_manifest.tsv` — one row per downloaded report; columns include `local_path`, `file_id`, `md5sum`, `file_size`, `drs_uri`, `role`, `pdc_study_id`

MD5 sums are cross-checked against the PDC API value at fetch time; a mismatch raises an error before any downstream step runs.

## Engine dependency

The pipeline requires a checkout of `dig-gene-set-extractors`, passed via `--dig_dir`. The runner invokes the engine as:

```
python -m geneset_extractors.cli <subcommand> [flags]
```

with `cwd=<dig_dir>` and `PYTHONPATH=<dig_dir>/src`. Two CPTAC-specific fixes must be present in the engine checkout:

- Lowercase phospho residues: the prepare workflow must accept `y`, `s`, `t` residue characters (not only uppercase)
- CDAP column-name join: the `"<aliquot_submitter_id> Log Ratio"` column header pattern must be stripped to a bare aliquot ID before the sample annotation join

Two engine subcommands are used per model run:

- prepare: `workflows ptm_prepare_public --input_mode cdap_files`
- extract: `convert ptm_site_matrix`

## Model PT1

`PT1` is the primary tumor-vs-adjacent-normal phosphoregulation model:

- contrast: `case` (`Primary Tumor`) vs `control` (`Solid Tissue Normal`) within a cohort
- prepare workflow: `geneset_extractors.cli workflows ptm_prepare_public`
  - `--ptm_type phospho`
  - `--assay_type_policy warn`
- extractor: `geneset_extractors.cli convert ptm_site_matrix`
  - `--study_contrast condition_a_vs_b`
  - `--condition_a case --condition_b control`
  - `--protein_adjustment_run_mode compare_if_protein`
  - `--select top_k --top_k 200`
  - `--gene_aggregation signed_topk_mean`
- `--protein_adjustment_run_mode compare_if_protein` produces two output variants: protein-adjusted and unadjusted; both are emitted when a protein matrix is available

## Software inputs

- `dig-gene-set-extractors/` checkout (the engine, Python)
- Python 3.x (stdlib-only runner; no third-party packages required for the runner itself)
- Engine dependencies as declared in the `dig-gene-set-extractors` environment

## Entry points

```
CPTAC/run/build_cptac_genesets.sh  [all runner args passed through]
CPTAC/src/run_cptac_ptm_model.py   [invoked directly for per-study runs]
```

Required arguments:

- `--dig_dir` — path to a `dig-gene-set-extractors` checkout
- `--cohort_id` — a `cohort_id` from `study_manifest.tsv` (e.g. `ccrcc`)
- `--out_root` — root directory for all outputs
- `--model_id` — model to run (default: `PT1`)

Optional:

- `--config_dir` — override the default `CPTAC/config/` directory
- `--python_bin` — Python interpreter to use (default: `python`)
- `--offline`, `--source_dir`, `--api_cache_json` — offline / air-gapped mode (see below)

## Output layout

The main output tree for a per-study, per-model run is:

- `<out_root>/genesets/<cohort_id>/fetch/` — downloaded reports, sample annotations, file manifest
- `<out_root>/genesets/<cohort_id>/models/<model_id>/prepared/` — engine-prepared PTM and protein matrices, sample metadata
- `<out_root>/genesets/<cohort_id>/models/<model_id>/extractor/` — gene sets (GMT), extractor outputs
- `<out_root>/genesets/<cohort_id>/models/<model_id>/provenance_overlay.json` — CRDC-native provenance overlay consumed by the extractor
- `<out_root>/genesets/<cohort_id>/models/<model_id>/run.log` — captured stdout/stderr from both engine invocations
- `<out_root>/genesets/<cohort_id>/models/<model_id>/commands.md` — the exact command lines that were run

## Offline / air-gapped path

Pass `--offline --source_dir <dir> --api_cache_json <file>` to skip all network calls:

- `--source_dir` must contain the phosphosite and proteome report TSVs under their canonical file names
- `--api_cache_json` must be a JSON file with keys `phospho_files`, `proteome_files`, and `biospecimen`, each holding the raw GraphQL response payload that would have been returned by the corresponding PDC query

MD5 verification still runs in offline mode against the values stored in the API cache. This is also the standard test path.

## Canonical checklist

Before running the pipeline end to end, confirm:

- [ ] `dig-gene-set-extractors` checkout is present and its path is known; the two CPTAC engine fixes (lowercase residues; CDAP Log Ratio column join) are included
- [ ] The target cohort is listed in `study_manifest.tsv` with `enabled=true`
- [ ] Network access to `https://pdc.cancer.gov/graphql` is available, OR offline assets (`--source_dir`, `--api_cache_json`) are staged
- [ ] `--out_root` directory exists or is creatable
- [ ] Python interpreter resolves `geneset_extractors.cli` from `<dig_dir>/src` (test with `PYTHONPATH=<dig_dir>/src python -c "import geneset_extractors"`)
