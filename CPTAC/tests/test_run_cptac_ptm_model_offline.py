import json
import shutil
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))
import run_cptac_ptm_model as runner  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[3]
DIG_DIR = REPO_ROOT / "dig-gene-set-extractors"
FIX = Path(__file__).resolve().parent / "fixtures"
CONFIG = Path(__file__).resolve().parents[1] / "config"


def _api_cache(tmp_path):
    import hashlib

    def md5(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()

    phospho = json.loads((FIX / "files_per_study.ccrcc_phospho.json").read_text())
    phospho["data"]["filesPerStudy"][0]["md5sum"] = md5(FIX / "reports" / "CCRCC_Phosphoproteome.phosphosite.tmt10.tsv")
    proteome = json.loads((FIX / "files_per_study.ccrcc_proteome.json").read_text())
    proteome["data"]["filesPerStudy"][0]["md5sum"] = md5(FIX / "reports" / "CCRCC_Proteome.tmt10.tsv")
    biospec = json.loads((FIX / "biospecimen_per_study.ccrcc.json").read_text())
    out = tmp_path / "api_cache.json"
    out.write_text(json.dumps({"phospho_files": phospho, "proteome_files": proteome, "biospecimen": biospec}))
    return out


@pytest.mark.skipif(not DIG_DIR.exists(), reason="dig-gene-set-extractors checkout not found")
def test_run_model_offline_end_to_end(tmp_path):
    out_root = tmp_path / "cptac_outputs"
    model_dir = runner.run_model(
        dig_dir=DIG_DIR,
        cohort_id="ccrcc",
        model_id="PT1",
        out_root=out_root,
        config_dir=CONFIG,
        offline=True,
        source_dir=FIX / "reports",
        api_cache_json=_api_cache(tmp_path),
        python_bin=sys.executable,
    )
    # At least one variant dir with the 3-file contract.
    provs = list(Path(model_dir).rglob("geneset.provenance.json"))
    assert provs, "no geneset.provenance.json produced"
    for prov in provs:
        assert (prov.parent / "geneset.tsv").exists()
        assert (prov.parent / "geneset.meta.json").exists()
    # The ptm_matrix.tsv File node carries the PDC phosphosite DRS id + PDC dcc_url.
    graph = json.loads(provs[0].read_text())
    graph = list(graph.values())[0] if "nodes" not in graph else graph
    file_nodes = [n for n in graph["nodes"] if n["type"] == "File" and n["name"] == "ptm_matrix.tsv"]
    assert file_nodes, "ptm_matrix.tsv File node missing"
    c = file_nodes[0]["c2m2_properties"]
    assert c["persistent_id"] == "11111111-1111-1111-1111-111111111111"
    assert c["local_id"] == "drs://dg.4DFC/11111111-1111-1111-1111-111111111111"
    assert file_nodes[0]["dcc_url"] == "https://pdc.cancer.gov/pdc/study/PDC000128"
