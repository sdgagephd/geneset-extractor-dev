"""Fetch CPTAC CDAP reports + provenance metadata from the PDC GraphQL API (stdlib only)."""
from __future__ import annotations

DRS_PREFIX = "drs://dg.4DFC/"

_TUMOR_TYPES = {"primary tumor", "tumor", "recurrent tumor", "metastatic"}
_NORMAL_TYPES = {"solid tissue normal", "normal", "blood derived normal", "adjacent normal"}


def classify_condition(sample_type: str) -> str:
    value = (sample_type or "").strip().lower()
    if value in _TUMOR_TYPES:
        return "case"
    if value in _NORMAL_TYPES:
        return "control"
    return ""


def parse_files_per_study(payload: dict, *, role: str, pdc_study_id: str) -> list[dict]:
    files = (payload.get("data") or {}).get("filesPerStudy") or []
    rows: list[dict] = []
    for f in files:
        file_id = str(f.get("file_id") or "")
        signed = f.get("signedUrl") or {}
        rows.append(
            {
                "file_id": file_id,
                "file_name": str(f.get("file_name") or ""),
                "md5sum": str(f.get("md5sum") or ""),
                "file_size": str(f.get("file_size") or ""),
                "data_category": str(f.get("data_category") or ""),
                "signed_url": str(signed.get("url") or ""),
                "role": role,
                "pdc_study_id": pdc_study_id,
                "drs_uri": DRS_PREFIX + file_id if file_id else "",
            }
        )
    return rows


def pick_report_file(rows: list[dict], *, kind: str) -> dict:
    """kind: 'phosphosite' or 'proteome'. Prefer the gene/site-level report TSV."""
    needles = {
        "phosphosite": ["phosphosite"],
        "proteome": ["proteome"],
    }[kind]
    avoid = {
        "phosphosite": ["phosphopeptide", "peptide", "spectral", "precursor", "mzid", "mzml"],
        "proteome": ["phospho", "phosphosite", "peptide", "spectral", "precursor", "mzid", "mzml"],
    }[kind]
    candidates = [
        r
        for r in rows
        if any(n in r["file_name"].lower() for n in needles)
        and not any(a in r["file_name"].lower() for a in avoid)
    ]
    if not candidates:
        raise ValueError(f"No {kind} report file found among {[r['file_name'] for r in rows]}")
    # Prefer the shortest matching name (the canonical '.phosphosite.tmtNN.tsv' / '.tmtNN.tsv').
    return sorted(candidates, key=lambda r: len(r["file_name"]))[0]


def parse_biospecimen(payload: dict) -> list[dict]:
    items = (payload.get("data") or {}).get("biospecimenPerStudy") or []
    rows: list[dict] = []
    for b in items:
        sample_type = str(b.get("sample_type") or "")
        rows.append(
            {
                "aliquot_submitter_id": str(b.get("aliquot_submitter_id") or ""),
                "sample_submitter_id": str(b.get("sample_submitter_id") or ""),
                "case_submitter_id": str(b.get("case_submitter_id") or ""),
                "sample_type": sample_type,
                "condition": classify_condition(sample_type),
            }
        )
    return rows


import argparse
import csv
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

PDC_GRAPHQL_URL = "https://pdc.cancer.gov/graphql"

_FILES_QUERY = """
{ filesPerStudy (study_id: "%s" acceptDUA: true) {
    file_id file_name md5sum file_size data_category file_type
    signedUrl { url } } }
"""

_STUDY_RESOLVE_QUERY = """
{ study (pdc_study_id: "%s" acceptDUA: true) {
    pdc_study_id study_id analytical_fraction } }
"""

_BIOSPEC_QUERY = """
{ biospecimenPerStudy (pdc_study_id: "%s" acceptDUA: true) {
    aliquot_submitter_id sample_submitter_id case_submitter_id sample_type } }
"""

_SAMPLE_ANNO_COLUMNS = [
    "sample_id_raw", "condition", "group", "study_id", "study_label",
    "case_submitter_id", "sample_submitter_id", "aliquot_submitter_id",
    "sample_type", "tissue_type",
]
_MANIFEST_COLUMNS = ["local_path", "file_id", "md5sum", "file_size", "drs_uri", "role", "pdc_study_id"]


def md5_file(path: str | Path) -> str:
    h = hashlib.md5()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def graphql_post(query: str) -> dict:
    body = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(
        PDC_GRAPHQL_URL, data=body, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=120) as resp:  # noqa: S310 (trusted PDC endpoint)
        return json.loads(resp.read().decode("utf-8"))


def parse_study_id(payload: dict) -> str:
    """Extract the (latest) study_id UUID from a PDC `study` query response."""
    studies = (payload.get("data") or {}).get("study") or []
    for s in studies:
        study_id = str(s.get("study_id") or "")
        if study_id:
            return study_id
    raise ValueError("could not resolve study_id from PDC study query response")


def resolve_study_id(pdc_study_id: str) -> str:
    """Resolve a pdc_study_id (e.g. PDC000127) to its current study_id UUID."""
    return parse_study_id(graphql_post(_STUDY_RESOLVE_QUERY % pdc_study_id))


def download_file(url: str, dest: str | Path) -> Path:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=600) as resp, dest.open("wb") as out:  # noqa: S310
        shutil.copyfileobj(resp, out)
    return dest


def write_sample_annotations(rows: list[dict], out_path: str | Path, *, study_id: str, study_label: str) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=_SAMPLE_ANNO_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "sample_id_raw": r["aliquot_submitter_id"],
                    "condition": r["condition"],
                    "group": study_id,
                    "study_id": study_id,
                    "study_label": study_label,
                    "case_submitter_id": r["case_submitter_id"],
                    "sample_submitter_id": r["sample_submitter_id"],
                    "aliquot_submitter_id": r["aliquot_submitter_id"],
                    "sample_type": r["sample_type"],
                    "tissue_type": r["condition"],
                }
            )


def write_file_manifest(rows: list[dict], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=_MANIFEST_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def _load_api_cache(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _acquire_report(
    *, file_rec: dict, out_dir: Path, source_dir: Path | None, offline: bool
) -> Path:
    dest = out_dir / file_rec["file_name"]
    if offline:
        src = Path(source_dir) / file_rec["file_name"]
        if not src.exists():
            raise FileNotFoundError(f"offline: expected report {src}")
        shutil.copyfile(src, dest)
    else:
        download_file(file_rec["signed_url"], dest)
    computed = md5_file(dest)
    if file_rec["md5sum"] and computed != file_rec["md5sum"]:
        raise ValueError(
            f"md5 mismatch for {file_rec['file_name']}: api={file_rec['md5sum']} computed={computed}"
        )
    return dest


def run_fetch(
    *,
    cohort_id: str,
    cohort_label: str,
    proteome_pdc_study_id: str,
    phospho_pdc_study_id: str,
    out_dir: str | Path,
    offline: bool = False,
    source_dir: str | Path | None = None,
    api_cache_json: str | Path | None = None,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if offline:
        if not api_cache_json:
            raise ValueError("offline mode requires api_cache_json")
        cache = _load_api_cache(api_cache_json)
        phospho_payload = cache["phospho_files"]
        proteome_payload = cache["proteome_files"]
        biospec_payload = cache["biospecimen"]
    else:
        phospho_study_uuid = resolve_study_id(phospho_pdc_study_id)
        proteome_study_uuid = resolve_study_id(proteome_pdc_study_id)
        phospho_payload = graphql_post(_FILES_QUERY % phospho_study_uuid)
        proteome_payload = graphql_post(_FILES_QUERY % proteome_study_uuid)
        biospec_payload = graphql_post(_BIOSPEC_QUERY % phospho_pdc_study_id)

    phospho_files = parse_files_per_study(phospho_payload, role="phospho", pdc_study_id=phospho_pdc_study_id)
    proteome_files = parse_files_per_study(proteome_payload, role="proteome", pdc_study_id=proteome_pdc_study_id)
    phospho_rec = pick_report_file(phospho_files, kind="phosphosite")
    proteome_rec = pick_report_file(proteome_files, kind="proteome")

    phospho_path = _acquire_report(file_rec=phospho_rec, out_dir=out_dir, source_dir=source_dir, offline=offline)
    proteome_path = _acquire_report(file_rec=proteome_rec, out_dir=out_dir, source_dir=source_dir, offline=offline)

    biospec_rows = parse_biospecimen(biospec_payload)
    sample_anno_path = out_dir / "sample_annotations.tsv"
    write_sample_annotations(biospec_rows, sample_anno_path, study_id=cohort_id, study_label=cohort_label)

    manifest_rows = [
        {**phospho_rec, "local_path": str(phospho_path)},
        {**proteome_rec, "local_path": str(proteome_path)},
    ]
    manifest_path = out_dir / "pdc_file_manifest.tsv"
    write_file_manifest(manifest_rows, manifest_path)

    return {
        "phospho_report": phospho_path,
        "proteome_report": proteome_path,
        "sample_annotations": sample_anno_path,
        "file_manifest": manifest_path,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Fetch CPTAC CDAP reports + provenance metadata from PDC.")
    p.add_argument("--cohort_id", required=True)
    p.add_argument("--cohort_label", required=True)
    p.add_argument("--proteome_pdc_study_id", required=True)
    p.add_argument("--phospho_pdc_study_id", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--offline", action="store_true")
    p.add_argument("--source_dir")
    p.add_argument("--api_cache_json")
    args = p.parse_args(argv)
    result = run_fetch(
        cohort_id=args.cohort_id,
        cohort_label=args.cohort_label,
        proteome_pdc_study_id=args.proteome_pdc_study_id,
        phospho_pdc_study_id=args.phospho_pdc_study_id,
        out_dir=args.out_dir,
        offline=args.offline,
        source_dir=args.source_dir,
        api_cache_json=args.api_cache_json,
    )
    print(json.dumps({k: str(v) for k, v in result.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
