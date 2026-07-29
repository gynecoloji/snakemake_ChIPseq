#!/usr/bin/env python3
"""Assemble a peak-set Jaccard matrix from pairwise bedtools jaccard output and
render a heatmap.

Input is the long-format TSV written by the peak_jaccard rule
(sample_a, sample_b, jaccard); output is a square matrix TSV + a PNG heatmap.
Pure-Python matrix assembly (matplotlib only for the plot) so it's unit-testable.
"""
import argparse
from pathlib import Path


def load_long(path):
    """Read (sample_a, sample_b, jaccard) rows; header line skipped."""
    rows = []
    lines = Path(path).read_text().splitlines()
    for ln in lines[1:]:
        if not ln.strip():
            continue
        f = ln.split("\t")
        try:
            j = float(f[2])
        except (IndexError, ValueError):
            j = 0.0
        if j != j:            # bedtools jaccard emits nan for empty-vs-empty
            j = 0.0
        rows.append((f[0], f[1], j))
    return rows


def to_matrix(rows):
    """Return (samples, matrix) with samples sorted; missing pairs = 0.0."""
    samples = sorted({a for a, _, _ in rows} | {b for _, b, _ in rows})
    idx = {s: i for i, s in enumerate(samples)}
    m = [[0.0] * len(samples) for _ in samples]
    for a, b, j in rows:
        m[idx[a]][idx[b]] = j
    return samples, m


def write_matrix(samples, matrix, path):
    lines = ["\t" + "\t".join(samples)]
    for s, row in zip(samples, matrix):
        lines.append(s + "\t" + "\t".join(f"{v:.4f}" for v in row))
    Path(path).write_text("\n".join(lines) + "\n")


def plot_heatmap(samples, matrix, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(samples)
    fig, ax = plt.subplots(figsize=(max(4, 0.6 * n + 3), max(3, 0.6 * n + 2)))
    im = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    short = [s.replace("GSF", "").split("-")[-1] if "-" in s else s for s in samples]
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(short, fontsize=8)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{matrix[i][j]:.2f}", ha="center", va="center",
                    color="white" if matrix[i][j] < 0.6 else "black", fontsize=7)
    ax.set_title("Peak-set Jaccard similarity")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=200)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Peak-set Jaccard matrix + heatmap.")
    ap.add_argument("long_tsv")
    ap.add_argument("-m", "--matrix", required=True)
    ap.add_argument("-p", "--plot", required=True)
    a = ap.parse_args(argv)
    samples, matrix = to_matrix(load_long(a.long_tsv))
    Path(a.matrix).parent.mkdir(parents=True, exist_ok=True)
    write_matrix(samples, matrix, a.matrix)
    plot_heatmap(samples, matrix, a.plot)


if "snakemake" in globals():  # pragma: no cover
    sm = snakemake  # noqa: F821
    _samples, _matrix = to_matrix(load_long(str(sm.input.long)))
    write_matrix(_samples, _matrix, str(sm.output.matrix))
    plot_heatmap(_samples, _matrix, str(sm.output.plot))
elif __name__ == "__main__":
    main()
