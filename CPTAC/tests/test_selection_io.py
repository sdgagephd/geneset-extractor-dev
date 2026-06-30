import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

import cptac_selection_io as sio  # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "config"


def test_study_manifest_has_ccrcc():
    studies = sio.load_study_manifest(CONFIG / "study_manifest.tsv")
    assert "ccrcc" in studies
    row = studies["ccrcc"]
    assert row["proteome_pdc_study_id"] == "PDC000127"
    assert row["phospho_pdc_study_id"] == "PDC000128"
    assert row["enabled"] == "true"


def test_models_merge_and_flags():
    models = sio.load_models(CONFIG / "model_list.tsv", CONFIG / "model_manifest.tsv")
    assert "PT1" in models
    pt1 = models["PT1"]
    assert pt1["model_family"] == "tumor_vs_normal"
    ex = sio.extractor_flags(pt1)
    assert ex["study_contrast"] == "condition_a_vs_b"
    assert ex["condition_a"] == "case"
    assert ex["protein_adjustment_run_mode"] == "compare_if_protein"
    pr = sio.prepare_flags(pt1)
    assert pr["ptm_type"] == "phospho"
    assert pr["assay_type_policy"] == "warn"
