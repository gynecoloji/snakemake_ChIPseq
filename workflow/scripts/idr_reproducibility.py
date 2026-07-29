#!/usr/bin/env python3
"""ENCODE-style reproducibility summary for 2-replicate ChIP-seq conditions.

Combines the true-replicate IDR peak count (Nt), self-pseudoreplicate counts
(N1, N2) and the pooled-pseudoreplicate count (Np) into:
  * self-consistency ratio = max(N1, N2) / min(N1, N2)
  * rescue ratio           = max(Nt, Np) / min(Nt, Np)
ENCODE flags a condition as reproducible when BOTH ratios are < 2.
"""
from pathlib import Path


def _count_lines(path):
    p = Path(path)
    if not p.exists():
        return 0
    return sum(1 for ln in p.read_text().splitlines() if ln.strip())


def _read_count(path):
    """Read a single integer count from a file (0 if missing/unparseable)."""
    try:
        return int(Path(path).read_text().split()[0])
    except (OSError, ValueError, IndexError):
        return 0


def ratio(a, b):
    """max/min ratio; inf if the smaller side is zero."""
    lo, hi = min(a, b), max(a, b)
    return (hi / lo) if lo > 0 else float("inf")


def flag(self_consistency, rescue):
    """ENCODE verdict: pass if both ratios < 2, warn if one is, else fail."""
    if self_consistency == float("inf") or rescue == float("inf"):
        return "fail"
    if self_consistency < 2 and rescue < 2:
        return "pass"
    if self_consistency < 2 or rescue < 2:
        return "warn"
    return "fail"


def build_rows(groups, members, true_idr_counts, self_counts, pool_counts):
    """Combine already-resolved counts into per-condition summary rows."""
    rows = []
    for g in groups:
        reps = members[g]
        n1 = self_counts.get(reps[0], 0)
        n2 = self_counts.get(reps[1], 0) if len(reps) > 1 else 0
        nt = true_idr_counts.get(g, 0)
        npool = pool_counts.get(g, 0)
        sc = ratio(n1, n2)
        rr = ratio(nt, npool)
        rows.append({"condition": g, "Nt": nt, "N1": n1, "N2": n2, "Np": npool,
                     "self_consistency_ratio": sc, "rescue_ratio": rr,
                     "status": flag(sc, rr)})
    return rows


def _fmt(x):
    return "inf" if x == float("inf") else f"{x:.2f}"


def _row_fields(r):
    return [r["condition"], str(r["Nt"]), str(r["N1"]), str(r["N2"]), str(r["Np"]),
            _fmt(r["self_consistency_ratio"]), _fmt(r["rescue_ratio"]), r["status"]]


def write_tables(rows, tsv, mqc):
    tlines = ["condition\tNt\tN1\tN2\tNp\tself_consistency_ratio\trescue_ratio\tstatus"]
    tlines += ["\t".join(_row_fields(r)) for r in rows]
    Path(tsv).write_text("\n".join(tlines) + "\n")

    header = ("# id: idr_reproducibility\n"
              "# section_name: 'IDR reproducibility (ENCODE)'\n"
              "# description: 'Self-consistency and rescue ratios; both < 2 = reproducible.'\n"
              "# plot_type: 'table'\n"
              "Condition\tNt\tN1\tN2\tNp\tSelf-consistency\tRescue\tStatus\n")
    body = "\n".join("\t".join(_row_fields(r)) for r in rows)
    Path(mqc).write_text(header + body + ("\n" if rows else ""))


if "snakemake" in globals():  # pragma: no cover
    sm = snakemake  # noqa: F821
    _groups = list(sm.params.groups)
    _members = {g: list(m) for g, m in dict(sm.params.members).items()}
    _true = {g: _count_lines(p) for g, p in dict(sm.params.true_idr).items()}
    _idir = str(sm.params.repro_idr_dir)
    _self, _pool = {}, {}
    for g in _groups:
        for s in _members[g]:
            _self[s] = _read_count(f"{_idir}/self__{s}.n.txt")
        _pool[g] = _read_count(f"{_idir}/pool__{g}.n.txt")
    _rows = build_rows(_groups, _members, _true, _self, _pool)
    write_tables(_rows, str(sm.output.tsv), str(sm.output.mqc))
