import hashlib
import json
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))
import fetch_pdc_study as fetch  # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures"
REPORTS = FIX / "reports"


def _md5(path: Path) -> str:
    h = hashlib.md5()
    h.update(path.read_bytes())
    return h.hexdigest()


def _build_api_cache(tmp_path: Path) -> Path:
    """Assemble an api_cache.json whose md5sums match the toy report files."""
    phospho = json.loads((FIX / "files_per_study.ccrcc_phospho.json").read_text())
    phospho["data"]["filesPerStudy"][0]["md5sum"] = _md5(REPORTS / "CCRCC_Phosphoproteome.phosphosite.tmt10.tsv")
    proteome = json.loads((FIX / "files_per_study.ccrcc_proteome.json").read_text())
    proteome["data"]["filesPerStudy"][0]["md5sum"] = _md5(REPORTS / "CCRCC_Proteome.tmt10.tsv")
    biospec = json.loads((FIX / "biospecimen_per_study.ccrcc.json").read_text())
    cache = {
        "phospho_files": phospho,
        "proteome_files": proteome,
        "biospecimen": biospec,
    }
    out = tmp_path / "api_cache.json"
    out.write_text(json.dumps(cache))
    return out


def test_run_fetch_offline(tmp_path):
    api_cache = _build_api_cache(tmp_path)
    out_dir = tmp_path / "fetch"
    result = fetch.run_fetch(
        cohort_id="ccrcc",
        cohort_label="Clear Cell RCC",
        proteome_pdc_study_id="PDC000127",
        phospho_pdc_study_id="PDC000128",
        out_dir=out_dir,
        offline=True,
        source_dir=REPORTS,
        api_cache_json=api_cache,
    )
    # Reports copied into fetch/
    assert result["phospho_report"].exists()
    assert result["proteome_report"].exists()
    # sample_annotations.tsv has 2 case + 2 control keyed by aliquot id
    anno = (out_dir / "sample_annotations.tsv").read_text().splitlines()
    header = anno[0].split("\t")
    assert "sample_id_raw" in header and "condition" in header
    body = anno[1:]
    conds = [line.split("\t")[header.index("condition")] for line in body]
    assert conds.count("case") == 2 and conds.count("control") == 2
    # pdc_file_manifest.tsv has DRS uris for both reports
    man = (out_dir / "pdc_file_manifest.tsv").read_text()
    assert "drs://dg.4DFC/11111111-1111-1111-1111-111111111111" in man
    assert "drs://dg.4DFC/33333333-3333-3333-3333-333333333333" in man


def test_run_fetch_offline_md5_mismatch_fails(tmp_path):
    api_cache = _build_api_cache(tmp_path)
    cache = json.loads(api_cache.read_text())
    cache["phospho_files"]["data"]["filesPerStudy"][0]["md5sum"] = "deadbeef"
    api_cache.write_text(json.dumps(cache))
    import pytest

    with pytest.raises(ValueError, match="md5"):
        fetch.run_fetch(
            cohort_id="ccrcc",
            cohort_label="Clear Cell RCC",
            proteome_pdc_study_id="PDC000127",
            phospho_pdc_study_id="PDC000128",
            out_dir=tmp_path / "fetch2",
            offline=True,
            source_dir=REPORTS,
            api_cache_json=api_cache,
        )
