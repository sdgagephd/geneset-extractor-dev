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
