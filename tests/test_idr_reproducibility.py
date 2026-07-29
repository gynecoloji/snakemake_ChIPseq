import importlib.util, pathlib

_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "idr_reproducibility.py"
_spec = importlib.util.spec_from_file_location("idr_reproducibility", _p)
ir = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ir)


def test_ratio_and_flag():
    assert ir.ratio(100, 100) == 1.0
    assert ir.ratio(100, 50) == 2.0
    assert ir.ratio(50, 100) == 2.0
    assert ir.ratio(0, 100) == float("inf")
    assert ir.flag(1.2, 1.5) == "pass"          # both < 2
    assert ir.flag(1.2, 2.5) == "warn"          # one < 2
    assert ir.flag(3.0, 4.0) == "fail"          # neither < 2
    assert ir.flag(float("inf"), 1.0) == "fail"


def test_build_rows_reproducible():
    rows = ir.build_rows(
        groups=["cJUN"],
        members={"cJUN": ["cJUN_1", "cJUN_2"]},
        true_idr_counts={"cJUN": 10000},
        self_counts={"cJUN_1": 9000, "cJUN_2": 9500},
        pool_counts={"cJUN": 10500},
    )
    r = rows[0]
    assert r["condition"] == "cJUN"
    assert r["Nt"] == 10000 and r["Np"] == 10500
    assert r["N1"] == 9000 and r["N2"] == 9500
    assert abs(r["self_consistency_ratio"] - 9500 / 9000) < 1e-9
    assert abs(r["rescue_ratio"] - 10500 / 10000) < 1e-9
    assert r["status"] == "pass"


def test_build_rows_flags_irreproducible():
    rows = ir.build_rows(
        groups=["bad"],
        members={"bad": ["a", "b"]},
        true_idr_counts={"bad": 100},
        self_counts={"a": 5000, "b": 500},       # self-consistency 10x -> fail
        pool_counts={"bad": 5000},               # rescue 50x -> fail
    )
    assert rows[0]["status"] == "fail"


def test_write_tables_roundtrip(tmp_path):
    rows = ir.build_rows(["g"], {"g": ["r1", "r2"]}, {"g": 10}, {"r1": 9, "r2": 8}, {"g": 11})
    tsv = tmp_path / "repro.tsv"
    mqc = tmp_path / "repro_mqc.txt"
    ir.write_tables(rows, str(tsv), str(mqc))
    head, r1 = tsv.read_text().splitlines()[:2]
    assert head.split("\t") == ["condition", "Nt", "N1", "N2", "Np",
                                "self_consistency_ratio", "rescue_ratio", "status"]
    assert r1.startswith("g\t10\t9\t8\t11\t")
    assert "# section_name: 'IDR reproducibility (ENCODE)'" in mqc.read_text()
