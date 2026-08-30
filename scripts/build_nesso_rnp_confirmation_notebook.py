#!/usr/bin/env python3
"""Build a concise notebook comparing the Nesso discovery and confirmation cohorts."""

from __future__ import annotations

import argparse
from pathlib import Path

import nbformat as nbf


def md(text: str):
    return nbf.v4.new_markdown_cell(text.strip() + "\n")


def code(text: str):
    return nbf.v4.new_code_cell(text.strip() + "\n")


def build() -> nbf.NotebookNode:
    notebook = nbf.v4.new_notebook()
    notebook.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3"},
    }
    notebook.cells = [
        md(
            r"""
# Nesso-1 on 50 additional Runs N' Poses systems

## Question

The original 100-system experiment showed that Nesso's structural predictions
were better for complexes resembling pre-September-2021 PDB complexes. It also
looked as though performance changed sharply near a familiarity score of 50.

We therefore selected **50 additional Runs N' Poses records**, excluding the
exact original 100 records, and repeated the same Nesso-1 inference and
coordinate-free structural scoring protocol.

The new cohort has unique clusters, PDB entries, and ligand SMILES internally.
Across cohorts, 14 structural clusters recur, so pooled uncertainty intervals
resample complete clusters rather than treating all 150 records as independent.

**Result:** the broad familiarity trend reproduces. A universal threshold at
50 does not.
"""
        ),
        code(
            r"""
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display


def repository_root(start: Path) -> Path:
    for candidate in (start.resolve(), *start.resolve().parents):
        if (candidate / "reports").is_dir() and (candidate / "configs").is_dir():
            return candidate
    raise FileNotFoundError("Run this notebook from inside BindingAffinityPredictor")


ROOT = repository_root(Path.cwd())
RESULTS = ROOT / "reports/exp013_nesso1_rnp_confirmation50/results"
comparison = json.loads((RESULTS / "comparison.json").read_text())
systems = pd.read_csv(RESULTS / "combined_per_system_metrics.csv")

assert comparison["cohort_overlap_count"] == 0
assert len(systems) == 150
systems.groupby("cohort").size()
"""
        ),
        md(
            r"""
## 1. What does Spearman $\rho$ tell us?

Spearman $\rho$ asks whether two quantities have a consistent rank ordering.
For errors such as distogram NLL and distance MAE, a negative value means
errors decrease as familiarity increases. For contact AP and pocket F1, a
positive value means accuracy increases with familiarity.

The confirmation cohort should be read independently of the original cohort.
If both show the same direction and their uncertainty intervals exclude zero,
the qualitative trend has replicated.
"""
        ),
        code(
            r"""
metrics = {
    "Distogram NLL ↓": "interface_distogram_nll",
    "Distance MAE ↓": "interface_expected_distance_mae_angstrom",
    "Contact AP ↑": "token_contact_average_precision_at_6a",
    "Pocket F1 ↑": "physical_pocket_f1_at_6a",
}
cohort_labels = {
    "discovery_100": "Original 100",
    "confirmation_50": "New 50",
    "combined_150": "Combined 150",
}
rows = []
for cohort, label in cohort_labels.items():
    aggregate = comparison["aggregates"][cohort]
    for metric_label, metric in metrics.items():
        result = aggregate["spearman_vs_familiarity"][metric]
        rows.append({
            "Cohort": label,
            "Metric": metric_label,
            "rho": result["rho"],
            "CI lower": result["lower_95"],
            "CI upper": result["upper_95"],
        })
rho_table = pd.DataFrame(rows)
display(rho_table.pivot(index="Metric", columns="Cohort", values="rho").round(3))
"""
        ),
        code(
            r"""
fig, ax = plt.subplots(figsize=(9.2, 4.6), constrained_layout=True)
x = np.arange(len(metrics))
width = 0.34
for offset, cohort, color in [
    (-width / 2, "discovery_100", "#4267B2"),
    (width / 2, "confirmation_50", "#E07A38"),
]:
    selected = rho_table[rho_table["Cohort"] == cohort_labels[cohort]].set_index("Metric")
    selected = selected.loc[list(metrics)]
    values = selected["rho"].to_numpy()
    low = values - selected["CI lower"].to_numpy()
    high = selected["CI upper"].to_numpy() - values
    ax.bar(x + offset, values, width, label=cohort_labels[cohort], color=color, alpha=0.9)
    ax.errorbar(x + offset, values, yerr=np.vstack([low, high]), fmt="none", color="#222", capsize=3)
ax.axhline(0, color="#333", linewidth=0.8)
ax.set_xticks(x, list(metrics), rotation=18, ha="right")
ax.set_ylabel("Spearman rho with familiarity")
ax.set_title("The familiarity trend reproduces in 50 new records", loc="left", weight="bold")
ax.legend(frameon=False)
ax.grid(axis="y", alpha=0.2)
ax.spines[["top", "right"]].set_visible(False)
plt.show()
"""
        ),
        md(
            r"""
All four confirmation-cohort correlations point in the expected direction,
and all four 95% bootstrap intervals exclude zero:

- distance MAE: $\rho=-0.630$;
- contact AP: $\rho=+0.529$;
- pocket F1: $\rho=+0.625$.

Nesso's interface entropy also continues to track distance error in the new 50
($\rho=+0.701$). This is evidence that the broad familiarity dependence was
not peculiar to the first selected sample.
"""
        ),
        md(
            r"""
## 2. Does performance jump specifically at 50?

For this diagnostic, we compare the median metric in the 40--50 bin with the
median in the 50--60 bin. This is an exploratory comparison, not a fitted or
pre-established biological threshold.
"""
        ),
        code(
            r"""
boundary_rows = []
for cohort in ("discovery_100", "confirmation_50", "combined_150"):
    for metric_label, metric in metrics.items():
        result = comparison["boundary_40_50_vs_50_60"][cohort][metric]
        boundary_rows.append({
            "Cohort": cohort_labels[cohort],
            "Metric": metric_label,
            "40–50 median": result["median_40_50"],
            "50–60 median": result["median_50_60"],
            "Difference": result["observed_difference"],
            "95% lower": result["bootstrap_lower_95"],
            "95% upper": result["bootstrap_upper_95"],
        })
boundary = pd.DataFrame(boundary_rows)
display(boundary[boundary["Cohort"] == "New 50"].set_index("Metric").round(3))
"""
        ),
        code(
            r"""
bins = ["0-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100"]
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0), constrained_layout=True)
panels = [
    ("interface_expected_distance_mae_angstrom", "Distance MAE (Å)", "lower is better"),
    ("token_contact_average_precision_at_6a", "Contact AP", "higher is better"),
    ("physical_pocket_f1_at_6a", "Pocket F1", "higher is better"),
]
for ax, (metric, title, direction) in zip(axes, panels):
    for cohort, color, marker in [
        ("discovery_100", "#4267B2", "o"),
        ("confirmation_50", "#E07A38", "s"),
    ]:
        aggregate = comparison["aggregates"][cohort]
        values = [aggregate["by_familiarity_bin"][b]["metrics"][metric]["median"] for b in bins]
        ax.plot(bins, values, marker=marker, color=color, label=cohort_labels[cohort])
    ax.axvline(3.5, color="#777", linestyle="--", linewidth=1)
    ax.set_title(title, loc="left", weight="bold")
    ax.text(0.98, 0.96, direction, transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].legend(frameon=False)
plt.show()
"""
        ),
        md(
            r"""
The exact jump did **not** replicate for distogram NLL, distance MAE, or
contact AP. In the new cohort, the 40--50 bin already performed well, and the
50--60 bin was similar or slightly worse. Pocket F1 did improve in both
cohorts.

## Conclusion

1. **Robust result:** Nesso's structural representation is substantially more
   accurate for systems structurally familiar to its training era.
2. **Robust uncertainty result:** higher Nesso interface entropy tends to mean
   larger structural error.
3. **Not supported:** a universal performance transition exactly at score 50.
4. **Possible narrower hypothesis:** pocket recovery may change strongly near
   this region, but testing that requires denser sampling and a predeclared
   changepoint analysis.
"""
        ),
    ]
    return notebook


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(build(), args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
