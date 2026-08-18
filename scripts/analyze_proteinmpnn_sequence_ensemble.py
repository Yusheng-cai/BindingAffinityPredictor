#!/usr/bin/env python3
"""Analyze identity, clustering, and position entropy for exp006."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

from affinity_benchmark.metrics.sequence_ensemble import (
    STANDARD_AMINO_ACIDS,
    pairwise_identity_matrix,
    position_frequencies,
    shannon_entropy,
)


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 14,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    }
)


IDENTITY_CUTS = (0.50, 0.60, 0.70, 0.80, 0.90)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences",
        type=Path,
        default=root
        / "runs/exp006_proteinmpnn_fixed_backbone_ensemble"
        / "proteinmpnn/fixed_backbone_cpu/sequences.csv",
    )
    parser.add_argument(
        "--backbone",
        type=Path,
        default=root
        / "runs/exp003_rfdiffusion_insr_binder_smoke"
        / "rfdiffusion/seed42/raw/design_ppi_42.pdb",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "reports/exp006_proteinmpnn_fixed_backbone_ensemble",
    )
    return parser.parse_args()


def load_sequences(path: Path) -> tuple[list[dict[str, str]], list[str], list[str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No sequence records in {path}")
    sequences = [row["sequence"] for row in rows]
    labels = [
        f"s{int(row['batch_seed'])}_{int(row['sample']):02d}" for row in rows
    ]
    if len(labels) != len(set(labels)):
        raise ValueError("Sample labels are not unique")
    return rows, sequences, labels


def parse_ca_coordinates(path: Path) -> dict[str, list[dict[str, object]]]:
    chains: dict[str, list[dict[str, object]]] = {}
    seen: set[tuple[str, int, str]] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            chain = line[21].strip()
            residue_number = int(line[22:26])
            insertion_code = line[26].strip()
            key = (chain, residue_number, insertion_code)
            if key in seen:
                continue
            seen.add(key)
            chains.setdefault(chain, []).append(
                {
                    "residue_number": residue_number,
                    "insertion_code": insertion_code,
                    "residue_name": line[17:20].strip(),
                    "xyz": np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    ),
                }
            )
    return chains


def descriptive(values: np.ndarray) -> dict[str, float]:
    return {
        "minimum": float(np.min(values)),
        "median": float(np.median(values)),
        "maximum": float(np.max(values)),
    }


def write_identity_csv(path: Path, labels: list[str], identity: np.ndarray) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *labels])
        for label, row in zip(labels, identity):
            writer.writerow([label, *[f"{value:.10f}" for value in row]])


def make_identity_figure(
    path: Path, identity: np.ndarray, tree: np.ndarray, order: np.ndarray
) -> None:
    figure = plt.figure(figsize=(10, 11))
    grid = figure.add_gridspec(2, 1, height_ratios=(1.2, 8), hspace=0.03)
    dendrogram_axis = figure.add_subplot(grid[0])
    dendrogram(
        tree,
        ax=dendrogram_axis,
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#4b5563",
    )
    dendrogram_axis.set_ylabel("Hamming\ndistance")
    dendrogram_axis.set_xticks([])
    dendrogram_axis.spines[["top", "right", "bottom"]].set_visible(False)

    heatmap_axis = figure.add_subplot(grid[1])
    ordered = identity[np.ix_(order, order)]
    image = heatmap_axis.imshow(
        ordered, cmap="viridis", vmin=0.35, vmax=1.0, interpolation="nearest"
    )
    heatmap_axis.set_xlabel("sequences in average-linkage order")
    heatmap_axis.set_ylabel("sequences in average-linkage order")
    heatmap_axis.set_xticks([])
    heatmap_axis.set_yticks([])
    colorbar = figure.colorbar(image, ax=heatmap_axis, fraction=0.046, pad=0.04)
    colorbar.set_label("position-wise sequence identity")
    figure.suptitle("ProteinMPNN fixed-backbone sequence identity", y=0.995)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def make_entropy_figure(
    path: Path, entropy: np.ndarray, consensus: list[str]
) -> None:
    positions = np.arange(1, len(entropy) + 1)
    figure, axis = plt.subplots(figsize=(13, 4.8))
    axis.plot(positions, entropy, color="#31688e", linewidth=1.8)
    axis.scatter(
        positions,
        entropy,
        c=entropy,
        cmap="viridis",
        vmin=0,
        vmax=max(entropy),
        s=28,
        zorder=3,
    )
    invariant = entropy < 1e-12
    axis.scatter(
        positions[invariant],
        entropy[invariant],
        color="#d1495b",
        marker="|",
        s=140,
        label=f"invariant in sample (n={int(invariant.sum())})",
        zorder=4,
    )
    for position, residue in zip(positions[invariant], np.asarray(consensus)[invariant]):
        axis.annotate(
            residue,
            (position, 0),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            fontsize=7,
        )
    axis.set(
        xlabel="binder backbone position",
        ylabel="empirical Shannon entropy (nats)",
        title="Position-wise variability across 100 ProteinMPNN samples",
        xlim=(0, len(entropy) + 1),
        ylim=(-0.04, max(entropy) * 1.08),
    )
    axis.legend(frameon=False, loc="upper right")
    axis.grid(axis="y", alpha=0.2)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def make_entropy_backbone_html(
    path: Path,
    chains: dict[str, list[dict[str, object]]],
    entropy: np.ndarray,
    consensus: list[str],
    consensus_frequency: np.ndarray,
) -> None:
    if "A" not in chains or "B" not in chains:
        raise ValueError("Expected target chain A and binder chain B")
    target = chains["A"]
    binder = chains["B"]
    if len(binder) != len(entropy):
        raise ValueError(
            f"Binder has {len(binder)} C-alpha atoms but {len(entropy)} entropy values"
        )
    target_xyz = np.stack([item["xyz"] for item in target])
    binder_xyz = np.stack([item["xyz"] for item in binder])
    hover = [
        (
            f"binder position {i + 1}<br>"
            f"consensus {consensus[i]} ({consensus_frequency[i]:.0%})<br>"
            f"entropy {entropy[i]:.3f} nats"
        )
        for i in range(len(entropy))
    ]
    figure = go.Figure()
    figure.add_trace(
        go.Scatter3d(
            x=target_xyz[:, 0],
            y=target_xyz[:, 1],
            z=target_xyz[:, 2],
            mode="lines",
            line={"color": "#b6c1bd", "width": 5},
            name="fixed receptor chain A",
            hoverinfo="skip",
        )
    )
    figure.add_trace(
        go.Scatter3d(
            x=binder_xyz[:, 0],
            y=binder_xyz[:, 1],
            z=binder_xyz[:, 2],
            mode="lines+markers",
            line={"color": "#334155", "width": 5},
            marker={
                "size": 6,
                "color": entropy,
                "colorscale": "Viridis",
                "cmin": 0,
                "cmax": float(max(entropy)),
                "colorbar": {"title": "entropy<br>(nats)"},
            },
            text=hover,
            hovertemplate="%{text}<extra></extra>",
            name="binder chain B",
        )
    )
    figure.update_layout(
        title=(
            "ProteinMPNN sequence variability on the fixed RFdiffusion backbone"
            "<br><sup>low entropy = constrained; high entropy = variable</sup>"
        ),
        scene={"aspectmode": "data"},
        height=760,
        margin={"l": 0, "r": 0, "t": 90, "b": 0},
    )
    figure.write_html(path, include_plotlyjs="cdn", full_html=True)


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows, sequences, labels = load_sequences(args.sequences)
    identity = pairwise_identity_matrix(sequences)
    off_diagonal = identity[np.triu_indices(len(sequences), k=1)]
    distance = 1.0 - identity
    tree = linkage(
        squareform(distance, checks=True), method="average", optimal_ordering=True
    )
    order = leaves_list(tree)

    frequencies = position_frequencies(sequences)
    entropy = shannon_entropy(frequencies)
    consensus_indices = np.argmax(frequencies, axis=1)
    consensus = [STANDARD_AMINO_ACIDS[index] for index in consensus_indices]
    consensus_frequency = frequencies[
        np.arange(len(consensus_indices)), consensus_indices
    ]

    cluster_assignments: dict[float, np.ndarray] = {}
    cluster_summary: dict[str, dict[str, object]] = {}
    for identity_cut in IDENTITY_CUTS:
        assignments = fcluster(tree, t=1.0 - identity_cut, criterion="distance")
        cluster_assignments[identity_cut] = assignments
        sizes = sorted(Counter(assignments).values(), reverse=True)
        cluster_summary[f"{identity_cut:.2f}"] = {
            "average_linkage_distance_cut": 1.0 - identity_cut,
            "cluster_count": len(sizes),
            "cluster_sizes_descending": sizes,
        }

    write_identity_csv(args.output / "pairwise_identity.csv", labels, identity)
    with (args.output / "sequence_clusters.csv").open("w", newline="") as handle:
        fields = [
            "sample_id",
            "batch_seed",
            "sample",
            "designed_chain_score",
            "dendrogram_order",
            *[f"cluster_at_{int(cut * 100)}pct_identity_cut" for cut in IDENTITY_CUTS],
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        rank = {int(index): position + 1 for position, index in enumerate(order)}
        for index, (label, row) in enumerate(zip(labels, rows)):
            output_row: dict[str, object] = {
                "sample_id": label,
                "batch_seed": row["batch_seed"],
                "sample": row["sample"],
                "designed_chain_score": row["designed_chain_score"],
                "dendrogram_order": rank[index],
            }
            for cut, assignments in cluster_assignments.items():
                output_row[f"cluster_at_{int(cut * 100)}pct_identity_cut"] = int(
                    assignments[index]
                )
            writer.writerow(output_row)

    entropy_fields = [
        "binder_position",
        "consensus_amino_acid",
        "consensus_frequency",
        "entropy_nats",
        "normalized_entropy_log20",
        "effective_amino_acids_exp_entropy",
        "invariant_in_100_samples",
        *[f"frequency_{residue}" for residue in STANDARD_AMINO_ACIDS],
    ]
    with (args.output / "position_entropy.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=entropy_fields)
        writer.writeheader()
        for index in range(len(entropy)):
            row: dict[str, object] = {
                "binder_position": index + 1,
                "consensus_amino_acid": consensus[index],
                "consensus_frequency": f"{consensus_frequency[index]:.10f}",
                "entropy_nats": f"{entropy[index]:.10f}",
                "normalized_entropy_log20": f"{entropy[index] / math.log(20):.10f}",
                "effective_amino_acids_exp_entropy": f"{math.exp(entropy[index]):.10f}",
                "invariant_in_100_samples": bool(entropy[index] < 1e-12),
            }
            row.update(
                {
                    f"frequency_{residue}": f"{frequencies[index, aa_index]:.10f}"
                    for aa_index, residue in enumerate(STANDARD_AMINO_ACIDS)
                }
            )
            writer.writerow(row)

    make_identity_figure(
        args.output / "identity_heatmap.png", identity, tree, order
    )
    make_entropy_figure(
        args.output / "position_entropy.png", entropy, consensus
    )
    make_entropy_backbone_html(
        args.output / "entropy_backbone.html",
        parse_ca_coordinates(args.backbone),
        entropy,
        consensus,
        consensus_frequency,
    )

    summary = {
        "samples": len(sequences),
        "sequence_length": len(sequences[0]),
        "identity_definition": (
            "Fraction of identical residues at corresponding positions; "
            "all sequences use the same fixed 90-residue backbone."
        ),
        "clustering": {
            "distance": "1 - position-wise identity",
            "linkage": "average",
            "hard_cut_sensitivity": cluster_summary,
        },
        "pairwise_identity_off_diagonal": descriptive(off_diagonal),
        "position_entropy_nats": descriptive(entropy),
        "entropy_estimator": (
            "Plug-in Shannon entropy from 100 sampled sequences; no finite-sample "
            "bias correction."
        ),
        "invariant_position_count": int(np.sum(entropy < 1e-12)),
        "invariant_positions": [
            index + 1 for index, value in enumerate(entropy) if value < 1e-12
        ],
        "dendrogram_sample_order": [labels[index] for index in order],
    }
    (args.output / "sequence_analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
