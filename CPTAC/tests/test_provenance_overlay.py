import json
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))
import build_pdc_provenance_overlay as ovl  # noqa: E402


def _manifest():
    return [
        {"local_path": "/x/CCRCC_Phosphoproteome.phosphosite.tmt10.tsv", "file_id": "111", "md5sum": "aaa", "drs_uri": "drs://dg.4DFC/111", "role": "phospho", "pdc_study_id": "PDC000128"},
        {"local_path": "/x/CCRCC_Proteome.tmt10.tsv", "file_id": "333", "md5sum": "ccc", "drs_uri": "drs://dg.4DFC/333", "role": "proteome", "pdc_study_id": "PDC000127"},
    ]


def test_overlay_keys_prepared_inputs_with_pdc_identity():
    overlay = ovl.build_overlay(
        manifest_rows=_manifest(),
        prepared_dir="/run/prepared",
        operation_meta={"script_url": "https://example/run.py"},
    )
    inputs = overlay["inputs"]
    ptm = inputs["/run/prepared/ptm_matrix.tsv"]
    assert ptm["persistent_id"] == "111"
    assert ptm["local_id"] == "drs://dg.4DFC/111"
    assert ptm["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000128"
    assert ptm["drc_url"] == ovl.CRDC_DRC_URL
    prot = inputs["/run/prepared/protein_matrix.tsv"]
    assert prot["persistent_id"] == "333"
    assert prot["local_id"] == "drs://dg.4DFC/333"
    assert prot["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000127"
    meta = inputs["/run/prepared/sample_metadata.tsv"]
    assert meta["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000128"
    assert overlay["gene_set"]["drc_url"] == ovl.CRDC_DRC_URL
    assert overlay["operation"]["script_url"] == "https://example/run.py"


def test_write_overlay_emits_files(tmp_path):
    res = ovl.write_overlay(
        manifest_rows=_manifest(),
        prepared_dir=str(tmp_path / "prepared"),
        operation_meta={"script_url": "https://example/run.py"},
        out_dir=tmp_path,
    )
    overlay = json.loads(Path(res["overlay_json"]).read_text())
    assert "/prepared/ptm_matrix.tsv" in "".join(overlay["inputs"].keys())
    smap = Path(res["source_map_tsv"]).read_text().splitlines()
    assert smap[0].split("\t") == ["local_path", "source_uri"]
    assert any("drs://dg.4DFC/111" in line for line in smap[1:])
