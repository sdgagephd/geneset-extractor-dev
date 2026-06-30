"""Per-study CPTAC PTM gene-set runner: fetch -> prepare -> overlay -> extract (stdlib only)."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cptac_selection_io as sio
import fetch_pdc_study as fetch
import build_pdc_provenance_overlay as ovl


def engine_cmd(dig_dir: Path, python_bin: str, *args: str) -> tuple[list[str], dict]:
    cmd = [python_bin, "-m", "geneset_extractors.cli", *args]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(dig_dir) / "src")
    return cmd, env


def _run(cmd: list[str], env: dict, dig_dir: Path, log_lines: list[str]) -> None:
    log_lines.append("$ " + " ".join(shlex.quote(c) for c in cmd))
    completed = subprocess.run(
        cmd, cwd=str(dig_dir), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False
    )
    log_lines.append(completed.stdout)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, cmd, output=completed.stdout)


def write_model_sidecars(
    *,
    extractor_dir: Path,
    model_id: str,
    cohort_label: str,
    phospho_pdc_study_id: str,
    proteome_pdc_study_id: str,
) -> None:
    """Emit one geneset.model.json per variant carrying the description-template vars.

    The shared refresh tool (`refresh_model_metadata_and_provenance.py`) renders model
    descriptions from `config/model_description_templates.tsv`, resolving `{model.<var>}`
    placeholders against this sidecar. Writing it here keeps CPTAC on the shared metadata
    tooling instead of a bespoke description pathway.
    """
    manifest = extractor_dir / "manifest.tsv"
    if not manifest.exists():
        return
    for row in sio.read_tsv(manifest):
        meta_rel = (row.get("meta_path") or "").strip()
        if not meta_rel:
            continue
        sidecar = (extractor_dir / meta_rel).parent / "geneset.model.json"
        sidecar.write_text(
            json.dumps(
                {
                    "model_id": model_id,
                    "cohort_label": cohort_label,
                    "variant_id": (row.get("variant_id") or "").strip(),
                    "phospho_pdc_study_id": phospho_pdc_study_id,
                    "proteome_pdc_study_id": proteome_pdc_study_id,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )


def run_model(
    *,
    dig_dir: str | Path,
    cohort_id: str,
    model_id: str,
    out_root: str | Path,
    config_dir: str | Path | None = None,
    offline: bool = False,
    source_dir: str | Path | None = None,
    api_cache_json: str | Path | None = None,
    python_bin: str = "python",
) -> Path:
    dig_dir = Path(dig_dir)
    config_dir = Path(config_dir) if config_dir else sio.default_config_dir()
    studies = sio.load_study_manifest(config_dir / "study_manifest.tsv")
    models = sio.load_models(config_dir / "model_list.tsv", config_dir / "model_manifest.tsv")
    study = studies[cohort_id]
    model = models[model_id]

    cohort_out = Path(out_root) / "genesets" / cohort_id
    model_out = cohort_out / "models" / model_id
    fetch_dir = cohort_out / "fetch"
    prepared_dir = model_out / "prepared"
    extractor_dir = model_out / "extractor"
    model_out.mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    # 1. fetch (cohort-level)
    fetched = fetch.run_fetch(
        cohort_id=cohort_id,
        cohort_label=study["cohort_label"],
        proteome_pdc_study_id=study["proteome_pdc_study_id"],
        phospho_pdc_study_id=study["phospho_pdc_study_id"],
        out_dir=fetch_dir,
        offline=offline,
        source_dir=source_dir,
        api_cache_json=api_cache_json,
    )

    try:
        # 2. prepare
        pflags = sio.prepare_flags(model)
        prepare_args = [
            "workflows", "ptm_prepare_public",
            "--input_mode", "cdap_files",
            "--ptm_report_tsv", str(fetched["phospho_report"]),
            "--protein_report_tsv", str(fetched["proteome_report"]),
            "--sample_annotations_tsv", str(fetched["sample_annotations"]),
            "--out_dir", str(prepared_dir),
            "--organism", "human",
            "--ptm_type", pflags.get("ptm_type", "phospho"),
            "--study_id", cohort_id,
            "--study_label", study["cohort_label"],
            "--assay_type_policy", pflags.get("assay_type_policy", "warn"),
        ]
        cmd, env = engine_cmd(dig_dir, python_bin, *prepare_args)
        _run(cmd, env, dig_dir, log_lines)

        # 3. overlay (read the manifest the fetch step wrote)
        manifest_rows = sio.read_tsv(fetched["file_manifest"])
        operation_meta = {
            "script_url": "https://github.com/broadinstitute/geneset-extractor-dev/blob/main/CPTAC/src/run_cptac_ptm_model.py",
        }
        written = ovl.write_overlay(
            manifest_rows=manifest_rows,
            prepared_dir=str(prepared_dir),
            operation_meta=operation_meta,
            out_dir=model_out,
        )

        # 4. extract
        eflags = sio.extractor_flags(model)
        extract_args = [
            "convert", "ptm_site_matrix",
            "--ptm_matrix_tsv", str(prepared_dir / "ptm_matrix.tsv"),
            "--sample_metadata_tsv", str(prepared_dir / "sample_metadata.tsv"),
            "--protein_matrix_tsv", str(prepared_dir / "protein_matrix.tsv"),
            "--study_contrast", eflags.get("study_contrast", "condition_a_vs_b"),
            "--condition_a", eflags.get("condition_a", "case"),
            "--condition_b", eflags.get("condition_b", "control"),
            "--protein_adjustment_run_mode", eflags.get("protein_adjustment_run_mode", "compare_if_protein"),
            # CPTAC CDAP proteome reports are gene-level (no per-row protein accession), so the
            # phospho-site -> proteome subtraction must join on gene symbol. Without this the
            # converter falls back to protein_accession, which is populated (RefSeq) on the phospho
            # side but empty on the gene-level proteome side -> 0 matched sites -> subtract == none.
            "--protein_accession_column", eflags.get("protein_accession_column", "gene_symbol"),
            "--select", eflags.get("select", "top_k"),
            "--top_k", eflags.get("top_k", "200"),
            "--gene_aggregation", eflags.get("gene_aggregation", "signed_topk_mean"),
            "--out_dir", str(extractor_dir),
            "--organism", "human",
            "--genome_build", "human",
            "--ptm_type", pflags.get("ptm_type", "phospho"),
            "--use_reference_bundle", "false",
            "--emit_small_gene_sets", "true",
            "--provenance_overlay_json", str(written["overlay_json"]),
        ]
        cmd, env = engine_cmd(dig_dir, python_bin, *extract_args)
        _run(cmd, env, dig_dir, log_lines)

        # 5. model sidecars (description-template vars for the shared refresh tooling)
        write_model_sidecars(
            extractor_dir=extractor_dir,
            model_id=model_id,
            cohort_label=study["cohort_label"],
            phospho_pdc_study_id=study["phospho_pdc_study_id"],
            proteome_pdc_study_id=study["proteome_pdc_study_id"],
        )

        return model_out
    finally:
        (model_out / "run.log").write_text("\n".join(log_lines), encoding="utf-8")
        (model_out / "commands.md").write_text(
            "# Commands\n\n```\n" + "\n".join(l for l in log_lines if l.startswith("$ ")) + "\n```\n",
            encoding="utf-8",
        )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run one CPTAC PTM gene-set model for one cohort.")
    p.add_argument("--dig_dir", required=True)
    p.add_argument("--cohort_id", required=True)
    p.add_argument("--model_id", default="PT1")
    p.add_argument("--out_root", required=True)
    p.add_argument("--config_dir")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--source_dir")
    p.add_argument("--api_cache_json")
    p.add_argument("--python_bin", default="python")
    args = p.parse_args(argv)
    model_out = run_model(
        dig_dir=args.dig_dir,
        cohort_id=args.cohort_id,
        model_id=args.model_id,
        out_root=args.out_root,
        config_dir=args.config_dir,
        offline=args.offline,
        source_dir=args.source_dir,
        api_cache_json=args.api_cache_json,
        python_bin=args.python_bin,
    )
    print(f"model_out={model_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
