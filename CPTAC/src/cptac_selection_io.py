"""Config IO + selection helpers for the CPTAC/PDC pipeline (stdlib only)."""
from __future__ import annotations

import csv
from pathlib import Path


def default_config_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "config"


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [{k: ("" if v is None else str(v)).strip() for k, v in row.items()} for row in reader]


def load_study_manifest(path: str | Path) -> dict[str, dict]:
    return {row["cohort_id"]: row for row in read_tsv(path)}


def load_models(model_list_path: str | Path, model_manifest_path: str | Path) -> dict[str, dict]:
    manifest = {row["model_id"]: row for row in read_tsv(model_manifest_path)}
    models: dict[str, dict] = {}
    for row in read_tsv(model_list_path):
        mid = row["model_id"]
        merged = dict(row)
        merged.update(manifest.get(mid, {}))
        models[mid] = merged
    return models


def _strip_prefix_flags(model: dict, prefix: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, value in model.items():
        if not key.startswith(prefix):
            continue
        if value == "" or value == "NA":
            continue
        out[key[len(prefix):]] = value
    return out


def extractor_flags(model: dict) -> dict[str, str]:
    return _strip_prefix_flags(model, "extractor_")


def prepare_flags(model: dict) -> dict[str, str]:
    return _strip_prefix_flags(model, "prepare_")


def enabled_ids(rows: dict[str, dict]) -> list[str]:
    return [k for k, v in rows.items() if v.get("enabled") == "true"]
