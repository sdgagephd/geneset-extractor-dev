import json
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))
import fetch_pdc_study as fetch  # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures"


def test_classify_condition():
    assert fetch.classify_condition("Primary Tumor") == "case"
    assert fetch.classify_condition("Solid Tissue Normal") == "control"
    assert fetch.classify_condition("") == ""


def test_parse_files_per_study():
    payload = json.loads((FIX / "files_per_study.ccrcc_phospho.json").read_text())
    rows = fetch.parse_files_per_study(payload, role="phospho", pdc_study_id="PDC000128")
    assert len(rows) == 2
    site = fetch.pick_report_file(rows, kind="phosphosite")
    assert site["file_name"] == "CCRCC_Phosphoproteome.phosphosite.tmt10.tsv"
    assert site["drs_uri"] == "drs://dg.4DFC/11111111-1111-1111-1111-111111111111"
    assert site["md5sum"] == "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    assert site["signed_url"] == "https://example.test/ccrcc_phospho.tsv"


def test_parse_biospecimen_maps_condition():
    payload = json.loads((FIX / "biospecimen_per_study.ccrcc.json").read_text())
    rows = fetch.parse_biospecimen(payload)
    by_aliquot = {r["aliquot_submitter_id"]: r for r in rows}
    assert by_aliquot["CPT0000010003"]["condition"] == "case"
    assert by_aliquot["CPT0000010004"]["condition"] == "control"
    assert sum(1 for r in rows if r["condition"] == "case") == 2
    assert sum(1 for r in rows if r["condition"] == "control") == 2


def test_parse_study_id_returns_uuid():
    payload = {"data": {"study": [
        {"pdc_study_id": "PDC000127",
         "study_id": "dbe94609-1fb3-11e9-b7f8-0a80fada099c",
         "analytical_fraction": "Proteome"}]}}
    assert fetch.parse_study_id(payload) == "dbe94609-1fb3-11e9-b7f8-0a80fada099c"


def test_parse_study_id_raises_when_absent():
    import pytest
    with pytest.raises(ValueError):
        fetch.parse_study_id({"data": {"study": []}})


def test_pick_report_file_selects_canonical_tmt_over_noise():
    # A realistic proteome file listing: only the .tmt10.tsv is the gene-level report.
    proteome_rows = [
        {"file_name": "Supplementary_Data_Proteome_TMT.tar.gz"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.peptides.tsv"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.qcmetrics.tsv"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.summary.tsv"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv"},
    ]
    assert (
        fetch.pick_report_file(proteome_rows, kind="proteome")["file_name"]
        == "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv"
    )

    # A phospho listing: must pick phosphosite, never phosphopeptide/peptides.
    phospho_rows = [
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.phosphopeptide.tmt10.tsv"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.peptides.tsv"},
        {"file_name": "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv"},
    ]
    assert (
        fetch.pick_report_file(phospho_rows, kind="phosphosite")["file_name"]
        == "CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv"
    )

    # proteome selector must NOT pick a phosphoproteome report.
    only_phospho = [{"file_name": "CPTAC3_X_Phosphoproteome.phosphosite.tmt10.tsv"}]
    import pytest
    with pytest.raises(ValueError):
        fetch.pick_report_file(only_phospho, kind="proteome")


def test_safe_name_strips_path_components():
    assert fetch._safe_name("CPTAC3_X_Proteome.tmt10.tsv") == "CPTAC3_X_Proteome.tmt10.tsv"
    assert fetch._safe_name("../../etc/passwd") == "passwd"
    assert fetch._safe_name("a/b/c.phosphosite.tmt10.tsv") == "c.phosphosite.tmt10.tsv"
