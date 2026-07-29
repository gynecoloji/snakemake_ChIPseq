import importlib.util, pathlib

_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "peak_overlap.py"
_spec = importlib.util.spec_from_file_location("peak_overlap", _p)
po = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(po)


def test_load_long(tmp_path):
    f = tmp_path / "long.tsv"
    f.write_text("sample_a\tsample_b\tjaccard\nA\tA\t1.0\nA\tB\t0.25\nB\tA\t0.25\nB\tB\t1.0\n")
    rows = po.load_long(str(f))
    assert ("A", "B", 0.25) in rows
    assert len(rows) == 4


def test_load_long_handles_bad_values(tmp_path):
    f = tmp_path / "long.tsv"
    f.write_text("sample_a\tsample_b\tjaccard\nA\tB\t\nB\tA\tnan\n")
    rows = po.load_long(str(f))
    # unparseable jaccard -> 0.0, never raises
    assert rows == [("A", "B", 0.0), ("B", "A", 0.0)]


def test_to_matrix_symmetric():
    rows = [("A", "A", 1.0), ("A", "B", 0.25), ("B", "A", 0.25), ("B", "B", 1.0)]
    samples, m = po.to_matrix(rows)
    assert samples == ["A", "B"]
    assert m == [[1.0, 0.25], [0.25, 1.0]]


def test_to_matrix_missing_pairs_zero():
    samples, m = po.to_matrix([("A", "A", 1.0), ("B", "B", 1.0)])
    assert samples == ["A", "B"]
    assert m[0][1] == 0.0 and m[1][0] == 0.0


def test_write_matrix_roundtrip(tmp_path):
    samples, m = po.to_matrix([("A", "A", 1.0), ("A", "B", 0.5), ("B", "A", 0.5), ("B", "B", 1.0)])
    out = tmp_path / "matrix.tsv"
    po.write_matrix(samples, m, str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == "\tA\tB"
    assert lines[1] == "A\t1.0000\t0.5000"
