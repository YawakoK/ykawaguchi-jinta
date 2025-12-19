BIN = 2_500_000
MAX_SCAFFOLDS = 51
USE_INFERRED_SCALE = False
USE_ORIENTATION = False
USE_LOG_SCALE = False

import bisect
import argparse
import os

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42  # embed TrueType so text stays text in PDF
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

CONTACT_CMAP = LinearSegmentedColormap.from_list(
    "contact_white_red_black_soft",
    [
        (0.0, "#ffffff"),
        (0.6, "#e82e2e"),
        (1.0, "#990808"),
    ],
)

parser = argparse.ArgumentParser(description="Plot Hi-C contact map with scaffold reordering.")
parser.add_argument(
    "-a",
    "--assembly0",
    required=True,
    help="Original assembly file (required)",
)
parser.add_argument(
    "-c",
    "--assembly-curated",
    required=True,
    help="Curated assembly file (required)",
)
parser.add_argument(
    "-m",
    "--hic-matrix",
    required=True,
    help="Hi-C matrix TSV (required)",
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="Output filename with extension (png or pdf, e.g., plot.png or plot.pdf)",
)
parser.add_argument(
    "-s",
    "--scale-factor",
    type=float,
    default=2,
    help="Scale factor to apply to matrix coordinates before binning (default: %(default)s)",
)
parser.add_argument(
    "--max-scaffolds",
    type=int,
    help="Limit number of curated scaffolds (0 for all)",
)
args = parser.parse_args()

assembly0 = args.assembly0
assembly_curated = args.assembly_curated
hic_matrix = args.hic_matrix
SCALE_FACTOR = args.scale_factor if args.scale_factor is not None else 2
MAX_SCAFFOLDS = None if args.max_scaffolds == 0 else args.max_scaffolds
OUTPUT_BASE = os.path.splitext(args.output)[0]
OUTPUT_EXT = os.path.splitext(args.output)[1].lstrip(".").lower()
OUTPUT_FORMATS = [OUTPUT_EXT] if OUTPUT_EXT else ["png"]


def load_contigs(path):
    contig_lengths = []
    starts_bp = []

    cur_bp = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                break
            _, _, length = line[1:].split()[:3]
            length = int(length)
            contig_lengths.append(length)
            starts_bp.append(cur_bp)
            cur_bp += length

    return contig_lengths, starts_bp, cur_bp


contig_lengths, starts_bp, genome_bp = load_contigs(assembly0)
n_contig = len(contig_lengths)


def read_scaffolds(path):
    scaffolds = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            vals = [int(x) for x in line.split()]
            scaffolds.append(vals)
    return scaffolds


scaffolds0 = read_scaffolds(assembly0)
scaffolds = read_scaffolds(assembly_curated)


def build_starts_from_scaffolds(scaffolds, contig_lengths):
    starts_bp = [-1] * len(contig_lengths)
    cur_bp = 0
    for contig_list in scaffolds:
        for signed_cid in contig_list:
            cid = abs(signed_cid)
            i = cid - 1
            if starts_bp[i] != -1:
                continue
            starts_bp[i] = cur_bp
            cur_bp += contig_lengths[i]
    for i, s in enumerate(starts_bp):
        if s == -1:
            starts_bp[i] = cur_bp
            cur_bp += contig_lengths[i]
    return starts_bp, cur_bp


# define genome bin coordinates based on assembly0 scaffold order
if scaffolds0:
    starts_bp, genome_bp = build_starts_from_scaffolds(scaffolds0, contig_lengths)

df = pd.read_csv(hic_matrix, sep="\t", names=["i", "j", "v"])

max_coord = max(df.i.max(), df.j.max())
matrix_bins_raw = max_coord // BIN + 1
expected_bins = (genome_bp + BIN - 1) // BIN
inferred_scale = max(1, round(expected_bins / matrix_bins_raw))
if SCALE_FACTOR is not None:
    scale_factor = SCALE_FACTOR
elif USE_INFERRED_SCALE:
    scale_factor = inferred_scale
else:
    scale_factor = 1
bin_size = BIN * scale_factor
matrix_bins = max_coord // bin_size + 1

print("genome_bp:", genome_bp)
print("expected bins @BIN:", expected_bins)
print("matrix bins (raw):", matrix_bins_raw)
print("scale factor inferred:", inferred_scale)
print("scale factor used:", scale_factor)
print("matrix bins (scaled):", matrix_bins)


def contig_bin_range(cid, bin_size):
    i = cid - 1
    start_bp = starts_bp[i]
    end_bp = start_bp + contig_lengths[i]
    start_bin = start_bp // bin_size
    end_bin = (end_bp - 1) // bin_size
    if end_bp <= start_bp:
        return start_bin, start_bin - 1
    return start_bin, end_bin


def contig_for_bp(bp):
    idx = bisect.bisect_right(starts_bp, bp) - 1
    if idx < 0:
        return 1
    if idx >= n_contig:
        return n_contig
    return idx + 1


def bins_by_contig(n_bins, bin_size):
    contig_bins = {}
    for b in range(n_bins):
        bp = b * bin_size
        cid = contig_for_bp(bp)
        contig_bins.setdefault(cid, []).append(b)
    return contig_bins


def contigs_in_range(scaff_list, bin_size, max_bin):
    order = []
    seen = set()
    for contig_list in scaff_list:
        for signed_cid in contig_list:
            cid = abs(signed_cid)
            start, end = contig_bin_range(cid, bin_size)
            if end < start or end < 0 or start > max_bin:
                continue
            if cid in seen:
                continue
            seen.add(cid)
            order.append(cid)
    return order


def build_bin_order(scaff_list, contig_bins):
    selected_bins = []
    edges = [0]
    for contig_list in scaff_list:
        for signed_cid in contig_list:
            cid = abs(signed_cid)
            bins = contig_bins.get(cid, [])
            if USE_ORIENTATION and signed_cid < 0:
                bins = list(reversed(bins))
            selected_bins.extend(bins)
        edges.append(len(selected_bins))
    return selected_bins, edges


def build_matrix(selected_bins, df):
    bin_set = set(selected_bins)
    sub = df[(df.bi.isin(bin_set)) & (df.bj.isin(bin_set))]
    bin_index = {b: k for k, b in enumerate(selected_bins)}
    n = len(selected_bins)
    mat = np.zeros((n, n), dtype=np.float32)
    for _, r in sub.iterrows():
        ii = bin_index[r.bi]
        jj = bin_index[r.bj]
        mat[ii, jj] += r.v
        if ii != jj:
            mat[jj, ii] += r.v
    return mat


def plot_matrix(mat, edges, title, bin_size, draw_edges=True, formats=("png",)):
    data = np.log1p(mat) if USE_LOG_SCALE else mat
    n = data.shape[0]
    extent = (0, n * bin_size, 0, n * bin_size)
    plt.figure(figsize=(8, 8))
    plt.imshow(
        data,
        cmap=CONTACT_CMAP,
        origin="lower",
        vmin=0,
        vmax=np.percentile(data, 99),
        interpolation="none",
        aspect="equal"
    )
    plt.colorbar(label="KR contact" if not USE_LOG_SCALE else "log1p(KR contact)")
    plt.title(title)
    ax = plt.gca()
    # tick positions in bin units, labels in Mb
    tick_bins = max(1, n // 10)
    xticks = np.arange(0, n, tick_bins)
    xticklabels = [f"{(t * bin_size) / 1_000_000:.0f}" for t in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yticks(xticks)
    ax.set_yticklabels(xticklabels)
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("Position (Mb)")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-0.5, n - 0.5)
    if draw_edges:
        for pos in edges[1:-1]:
            plt.axhline(pos, color="white", lw=0.6)
            plt.axvline(pos, color="white", lw=0.6)
    plt.tight_layout()
    for ext in formats:
        out = f"{OUTPUT_BASE}.{ext}"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        print("saved:", out)
    plt.close()


# prepare bin columns in matrix (scale coordinates if needed)
df["bi"] = (df.i * scale_factor) // bin_size
df["bj"] = (df.j * scale_factor) // bin_size
n_full = int(max(df.bi.max(), df.bj.max()) + 1)
contig_bins = bins_by_contig(n_full, bin_size)

# raw assembly order (from assembly0 scaffolds)
raw_scaffolds = scaffolds0 if scaffolds0 else [[cid + 1] for cid in range(n_contig)]
raw_bins, raw_edges = build_bin_order(raw_scaffolds, contig_bins)
print("raw bins:", len(raw_bins), "matrix bins:", n_full)
raw_mat = build_matrix(raw_bins, df)
plot_matrix(
    raw_mat,
    raw_edges,
    f"Hi-C: assembly order (bin={bin_size//1_000_000}Mb x{scale_factor})",
    bin_size=bin_size,
    draw_edges=False,
    formats=OUTPUT_FORMATS,
)

print("scaffolds (all):", len(scaffolds))
if MAX_SCAFFOLDS is not None:
    scaffolds = scaffolds[:MAX_SCAFFOLDS]
    print("using first scaffolds:", len(scaffolds))

# contig order sanity for visible bins
raw_contig_order = contigs_in_range(raw_scaffolds, bin_size, n_full - 1)
cur_contig_order = contigs_in_range(scaffolds, bin_size, n_full - 1)
contig_diff = sum(
    1 for a, b in zip(raw_contig_order, cur_contig_order) if a != b
)
print("contig order diff (raw vs curated):", contig_diff, "of", len(raw_contig_order))
print("raw contigs (first 20):", raw_contig_order[:20])
print("curated contigs (first 20):", cur_contig_order[:20])

# show scaffolds where order differs within visible bins
def scaffold_signature(scaff_list, bin_size, max_bin):
    sig = []
    for contig_list in scaff_list:
        contigs = []
        for signed_cid in contig_list:
            cid = abs(signed_cid)
            start, end = contig_bin_range(cid, bin_size)
            if end < start or end < 0 or start > max_bin:
                continue
            contigs.append(cid)
        sig.append(contigs)
    return sig

raw_sig = scaffold_signature(raw_scaffolds, bin_size, n_full - 1)
cur_sig = scaffold_signature(scaffolds, bin_size, n_full - 1)
diff_idx = [
    i for i, (a, b) in enumerate(zip(raw_sig, cur_sig))
    if a != b and (a or b)
]
print("scaffold order diff count:", len(diff_idx))
print("scaffold order diff idx (first 10):", diff_idx[:10])
for i in diff_idx[:5]:
    print("raw scaffold", i + 1, raw_sig[i])
    print("cur scaffold", i + 1, cur_sig[i])

# curated scaffold order
cur_bins, cur_edges = build_bin_order(scaffolds, contig_bins)
print("curated bins:", len(cur_bins))
if len(raw_bins) == len(cur_bins):
    diff = sum(1 for a, b in zip(raw_bins, cur_bins) if a != b)
    print("bin order diff (raw vs curated):", diff, "of", len(raw_bins))
cur_mat = build_matrix(cur_bins, df)
plot_matrix(
    cur_mat,
    cur_edges,
    f"Hi-C: curated scaffolds (bin={bin_size//1_000_000}Mb x{scale_factor})",
    bin_size=bin_size,
    draw_edges=False,
    formats=OUTPUT_FORMATS,
)
