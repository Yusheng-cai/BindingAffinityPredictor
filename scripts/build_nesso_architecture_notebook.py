#!/usr/bin/env python3
"""Build the executed-ready Nesso-1 architecture companion notebook."""

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "weeks" / "2026-W35" / "notebooks" / "nesso_architecture_visualized.ipynb"


def markdown(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


notebook = nbf.v4.new_notebook()
notebook["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3"},
}

notebook["cells"] = [
    markdown(
        r"""
# Nesso-1 architecture, made visible

This companion notebook turns the architecture discussion in the Week 35 report into inspectable tensors and interactive plots. It follows the path

$$
\text{protein sequence + ligand graph}
\longrightarrow s_i
\longrightarrow z_{ij}^{(0)}
\longrightarrow \text{Pairformer}
\longrightarrow \text{distogram/crop}
\longrightarrow \text{affinity ensemble}.
$$

There are two deliberately separated parts:

1. a **tiny educational calculation** with transparent, randomly initialized projections, used only to explain the mathematics; and
2. a **released-checkpoint trace** loaded from one saved Nesso-1 prediction, used to show the actual tensor dimensions and outputs.

> **Scientific boundary.** The miniature network has no Nesso-1 weights and produces no scientifically meaningful affinity. It mirrors information flow, not checkpoint performance. The final section does not rerun inference; it visualizes raw tensors already saved from the released Nesso-1 checkpoint.

Primary sources: [Nesso repository](https://github.com/Genentech/nesso), the Nesso-1 paper included in this project's literature collection, and the Pairformer description in [AlphaFold 3](https://doi.org/10.1038/s41586-024-07487-w).
"""
    ),
    markdown(
        r"""
## 1. The architecture at one glance

Nesso represents each protein residue and each ligand atom as a **token**. A token vector $s_i$ describes item $i$; a pair vector $z_{ij}$ describes the ordered relationship from token $i$ to token $j$.

| stage | released Nesso-1 width / setting | interpretation |
|---|---:|---|
| ESM-2 residue feature | 1280 | sequence context learned before Nesso training |
| Nesso input token $s_i^{\rm input}$ | 384 | atom/residue chemistry encoded by Nesso |
| initial pair vector $z_{ij}$ | 128 | relationship for every ordered token pair |
| trunk | 48 Pairformer blocks | repeatedly exchange information among pair entries |
| distance head | 64 bins | probability distribution over inter-token distance |
| refinement crop | within predicted 22 Å; at most 256 tokens | focus the expensive all-atom refinement |
| affinity crop | predicted 15 Å neighborhood | retain the putative binding region |
| affinity ensemble | two independent 8-block Pairformer modules | two continuous estimates, then their mean |

The widths are **not atom coordinates**. They are learned feature channels. Nesso does not replace its 384-channel input token with ESM-2. It projects the 1280-channel ESM-2 protein context to 384, combines that result with Nesso's own input token inside a separate ESM module, and converts the combination into a pair update for $z$. The affinity modules also receive a separately projected copy of the ESM context.
"""
    ),
    code(
        """
from pathlib import Path
import json
import math
import sys

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import torch
from safetensors.torch import load_file

REPO_ROOT = Path.cwd().resolve()
while REPO_ROOT != REPO_ROOT.parent and not (REPO_ROOT / "src").exists():
    REPO_ROOT = REPO_ROOT.parent
if not (REPO_ROOT / "src").exists():
    raise RuntimeError("Run this notebook from inside the BindingAffinityPredictor repository")
sys.path.insert(0, str(REPO_ROOT / "src"))

from affinity_benchmark.educational.mini_nesso import (
    expected_distance_and_entropy,
    initialize_pair_representation,
    one_mini_pairformer_step,
    pool_affinity_pairs,
)

torch.set_num_threads(1)
torch.manual_seed(12)
COLORS = {"protein": "#3978c5", "ligand": "#df7b2c"}
print(f"PyTorch {torch.__version__}; CPU teaching calculation; repository={REPO_ROOT}")
"""
    ),
    markdown(
        r"""
## 2. Build a nine-token protein–ligand example

We use six imaginary protein residues and a three-atom ligand. The ligand tokens are atoms because small-molecule chemistry is naturally described by an atom/bond graph; the protein tokens are residues because ESM-2 and the main trunk operate at residue resolution.

First, the Nesso input embedder creates $s_i^{\rm input}\in\mathbb{R}^{384}$ from atom-level and residue-identity features. Separately, for a protein residue, the ESM module computes

$$
\widetilde s_i = W_{\rm input}s_i^{\rm input}
+ \operatorname{MLP}_{\rm ESM}(e_i^{\rm ESM}),
\qquad e_i^{\rm ESM}\in\mathbb{R}^{1280},\quad \widetilde s_i\in\mathbb{R}^{384}.
$$

For a ligand atom the ESM term is zero, so ligand chemistry comes from $s_i^{\rm input}$. The combined $\widetilde s_i$ is projected into a 128-channel update for the pair tensor. This separation matters: ESM supplies sequence context, while Nesso's input embedder supplies explicit atom/residue chemistry.
"""
    ),
    code(
        """
labels = ["P1 ALA", "P2 LYS", "P3 TYR", "P4 SER", "P5 LEU", "P6 ASP",
          "L1 C", "L2 N", "L3 O"]
n_protein, n_ligand = 6, 3
n_tokens = len(labels)
is_protein = torch.tensor([True] * n_protein + [False] * n_ligand)
is_ligand = ~is_protein

g = torch.Generator().manual_seed(5)
esm_1280 = torch.randn(n_protein, 1280, generator=g)
w_esm = torch.randn(1280, 384, generator=g) / math.sqrt(1280)
input_tokens = torch.randn(n_tokens, 384, generator=g) * 0.35
projected_esm = torch.zeros(n_tokens, 384)
projected_esm[:n_protein] = esm_1280 @ w_esm
combined_tokens = input_tokens + projected_esm

fig = go.Figure()

def add_flow_box(x, y, text, color, width=1.25):
    fig.add_shape(type="rect", x0=x-width/2, x1=x+width/2, y0=y-0.34, y1=y+0.34,
                  line=dict(color=color, width=2), fillcolor=color, opacity=0.20,
                  layer="below")
    fig.add_annotation(x=x, y=y, text=text, showarrow=False, align="center",
                       font=dict(size=13, color="#20242a"))

def add_flow_arrow(x0, y0, x1, y1, color):
    fig.add_annotation(x=x1, y=y1, ax=x0, ay=y0,
                       xref="x", yref="y", axref="x", ayref="y",
                       text="", showarrow=True, arrowhead=3, arrowsize=1.2,
                       arrowwidth=2.2, arrowcolor=color)

# Upper lane: protein sequence context from ESM-2.
add_flow_box(0.7, 2.35, "protein<br>sequence", COLORS["protein"])
add_flow_box(2.25, 2.35, "ESM-2<br>1280 / residue", "#719bd0", width=1.35)
add_flow_box(3.9, 2.35, "ESM MLP<br>1280 → 384", "#777777", width=1.35)
add_flow_arrow(1.35, 2.35, 1.55, 2.35, COLORS["protein"])
add_flow_arrow(2.95, 2.35, 3.20, 2.35, COLORS["protein"])

# Lower lane: chemistry encoded by Nesso itself.
add_flow_box(0.85, 0.45, "atom + residue<br>chemistry", COLORS["ligand"], width=1.55)
add_flow_box(2.85, 0.45, "Nesso input sᵢ<br>384 / token", "#df9d60", width=1.55)
add_flow_arrow(1.65, 0.45, 2.05, 0.45, COLORS["ligand"])

# The two lanes meet once, then update the pair tensor.
add_flow_box(5.55, 1.40, "combined s̃ᵢ<br>384 / token", "#6d55a3", width=1.55)
add_flow_box(7.35, 1.40, "left + right<br>projections", "#4f9b75", width=1.35)
add_flow_box(9.05, 1.40, "pair update Δz<br>128 channels", "#4f9b75", width=1.55)
add_flow_arrow(4.60, 2.18, 4.78, 1.75, COLORS["protein"])
add_flow_arrow(3.68, 0.62, 4.76, 1.10, COLORS["ligand"])
add_flow_arrow(6.35, 1.40, 6.65, 1.40, "#4f9b75")
add_flow_arrow(8.05, 1.40, 8.25, 1.40, "#4f9b75")

fig.add_annotation(x=0.0, y=2.95, text="PROTEIN SEQUENCE CONTEXT", showarrow=False,
                   xanchor="left", font=dict(size=12, color=COLORS["protein"]))
fig.add_annotation(x=0.0, y=-0.12, text="NESSO CHEMISTRY FEATURES (PROTEIN + LIGAND)", showarrow=False,
                   xanchor="left", font=dict(size=12, color=COLORS["ligand"]))
fig.update_xaxes(range=[-0.15, 10.0], visible=False, fixedrange=True)
fig.update_yaxes(range=[-0.35, 3.15], visible=False, fixedrange=True)
fig.update_layout(
    title="Two separate feature routes meet once before updating the pair tensor",
    height=500, margin=dict(l=25, r=25, t=70, b=25), plot_bgcolor="white",
)
fig.show()

fig = go.Figure(go.Bar(
    x=labels, y=combined_tokens.norm(dim=-1).numpy(),
    marker_color=[COLORS["protein"]] * n_protein + [COLORS["ligand"]] * n_ligand,
    customdata=np.array(["protein residue"] * n_protein + ["ligand atom"] * n_ligand),
    hovertemplate="%{x}<br>%{customdata}<br>||sᵢ||=%{y:.2f}<extra></extra>",
))
fig.update_layout(title="Toy combined ESM-module token magnitudes (random teaching values)",
                  xaxis_title="token", yaxis_title="Euclidean norm of 384 channels", height=420)
fig.show()
print("Nesso input tokens:", tuple(input_tokens.shape))
print("Protein ESM input:", tuple(esm_1280.shape), "→ combined ESM-module tokens:", tuple(combined_tokens.shape))
"""
    ),
    markdown(
        r"""
## 3. Initialize the pair tensor $z_{ij}^{(0)}$

The initial pair representation is more than a dot product. The released source first constructs

$$
z_{ij}^{\rm init} = W_Ls_i^{\rm input} + W_Rs_j^{\rm input}
 + W_{\rm rel}r_{ij} + W_{\rm bond}b_{ij}.
$$

$r_{ij}$ records chain/sequence relationships and $b_{ij}$ records ligand bonds and bond types. During every recycling pass, the ESM module then adds

$$
\Delta z_{ij}^{\rm ESM}=U_L\widetilde s_i+U_R\widetilde s_j,
\qquad z_{ij}\leftarrow z_{ij}^{\rm init}+\Delta z_{ij}^{\rm ESM}
$$

(plus the recycled previous pair state after the first pass). Because $z_{ij}$ is ordered, $z_{ij}$ need not equal $z_{ji}$.
"""
    ),
    code(
        """
relative = torch.zeros(n_tokens, n_tokens, 4)
bond = torch.zeros(n_tokens, n_tokens, 2)
for i in range(n_tokens):
    for j in range(n_tokens):
        relative[i, j, 0] = float(is_protein[i] and is_protein[j])
        relative[i, j, 1] = float(is_ligand[i] and is_ligand[j])
        relative[i, j, 2] = (j - i) / max(n_tokens - 1, 1)
        relative[i, j, 3] = float(i == j)
bond[6, 7, 0] = bond[7, 6, 0] = 1.0       # toy single bond C--N
bond[7, 8, 1] = bond[8, 7, 1] = 1.0       # toy double bond N==O

z_base = initialize_pair_representation(input_tokens, relative, bond, pair_dim=128, seed=9)
g = torch.Generator().manual_seed(10)
w_esm_left = torch.randn(384, 128, generator=g) / math.sqrt(384)
w_esm_right = torch.randn(384, 128, generator=g) / math.sqrt(384)
delta_z_esm = (combined_tokens @ w_esm_left)[:, None, :] + (combined_tokens @ w_esm_right)[None, :, :]
z0 = z_base + delta_z_esm
z0_norm = z0.norm(dim=-1)
fig = go.Figure(go.Heatmap(
    z=z0_norm.numpy(), x=labels, y=labels, colorscale="Viridis",
    colorbar=dict(title="||zᵢⱼ||"), hovertemplate="i=%{y}<br>j=%{x}<br>norm=%{z:.2f}<extra></extra>",
))
fig.add_vline(x=n_protein - 0.5, line_color="white", line_width=3)
fig.add_hline(y=n_protein - 0.5, line_color="white", line_width=3)
fig.update_layout(title="Initial ordered-pair tensor: protein/protein, protein/ligand, ligand/protein, ligand/ligand blocks",
                  height=620, xaxis_title="destination token j", yaxis_title="source token i",
                  yaxis_autorange="reversed")
fig.show()
print("z initial:", tuple(z_base.shape), "+ ESM pair update:", tuple(delta_z_esm.shape))
print("Pairformer input z shape:", tuple(z0.shape), "= [token i, token j, pair channel]")
"""
    ),
    markdown(
        r"""
## 4. Why outgoing and incoming triangle multiplication are different

For a target pair $(i,j)$, Pairformer asks what every third token $k$ says about that relationship. The two directions use different edges:

$$
\Delta z_{ij}^{\rm out}=\sum_k A(z_{ik})\odot B(z_{jk}),
\qquad
\Delta z_{ij}^{\rm in}=\sum_k A(z_{ki})\odot B(z_{kj}).
$$

- **Outgoing:** $i$ and $j$ both point toward $k$.
- **Incoming:** $k$ points toward both $i$ and $j$.

The learned projections $A$ and $B$ decide what kind of evidence can combine. These updates are not literal forces or measured distances; they are feature updates that let one pair learn from triangles in the interaction graph.

### What happens to the two updates?

They are added through **sequential residual connections**. Let $z^{[b]}$ mean the pair tensor entering Pairformer block $b$. At inference, where dropout is inactive, the released block follows

$$
\begin{aligned}
z^{\rm out} &= z^{[b]} + \Delta z^{\rm out}\!\left(z^{[b]}\right),\\
z^{\rm in}  &= z^{\rm out} + \Delta z^{\rm in}\!\left(z^{\rm out}\right),\\
z^{\rm start} &= z^{\rm in} + \Delta z^{\rm attn,start}\!\left(z^{\rm in}\right),\\
z^{\rm end} &= z^{\rm start} + \Delta z^{\rm attn,end}\!\left(z^{\rm start}\right),\\
z^{[b+1]} &= z^{\rm end} + \operatorname{Transition}\!\left(z^{\rm end}\right).
\end{aligned}
$$

So the short answer is **yes, the original pair information remains through the residual path**, but the incoming update is not independently calculated from the untouched $z^{[b]}$. It is recalculated from $z^{\rm out}$, which already contains the outgoing correction:

$$
z^{\rm in}=z^{[b]}+\Delta z^{\rm out}(z^{[b]})
+\Delta z^{\rm in}\!\left(z^{[b]}+\Delta z^{\rm out}(z^{[b]})\right).
$$

### What do “starting-node” and “ending-node” attention mean?

The names refer to the two ends of the ordered pair $i\rightarrow j$.

- **Starting-node attention** fixes the first endpoint $i$. The query pair $z_{ij}$ attends across the row $\{z_{ik}:k=1,\ldots,N\}$: all relationships that leave $i$.
- **Ending-node attention** fixes the second endpoint $j$. The query pair $z_{ij}$ attends down the column $\{z_{kj}:k=1,\ldots,N\}$: all relationships that enter $j$.

Ignoring multi-head indexing, pair biases, masks, and learned gates, the two operations can be written schematically as

$$
\alpha_{ijk}=\operatorname{softmax}_k
\left(\frac{Q(z_{ij})\cdot K(z_{ik})}{\sqrt d}\right),
\qquad
\Delta z_{ij}^{\rm attn,start}
=W_O\sum_k\alpha_{ijk}V(z_{ik}),
$$

$$
\beta_{ijk}=\operatorname{softmax}_k
\left(\frac{Q(z_{ij})\cdot K(z_{kj})}{\sqrt d}\right),
\qquad
\Delta z_{ij}^{\rm attn,end}
=W_O\sum_k\beta_{ijk}V(z_{kj}).
$$

Thus, starting attention asks what the neighborhood around $i$ says about $i\rightarrow j$, while ending attention asks what the neighborhood around $j$ says about it. The softmax weights determine internal information routing; they are not physical contact or binding probabilities.

The resulting $z^{[b+1]}$ becomes the input to the next Pairformer block. The released Nesso trunk repeats this across 48 blocks. The miniature below keeps outgoing, incoming, starting-node attention, and the transition; it omits ending-node attention to keep the code short.
"""
    ),
    code(
        """
trace = one_mini_pairformer_step(z0, seed=13)
i, j = 2, 7  # P3 TYR -> L2 N
updates = [
    ("initial ||z⁽⁰⁾||", trace.z_initial.norm(dim=-1)),
    ("||Δ outgoing||", trace.delta_outgoing.norm(dim=-1)),
    ("||Δ incoming||", trace.delta_incoming.norm(dim=-1)),
    ("||Δ attention||", trace.delta_attention.norm(dim=-1)),
    ("total ||z final - z⁽⁰⁾||", (trace.z_final - trace.z_initial).norm(dim=-1)),
]
fig = make_subplots(rows=2, cols=3, subplot_titles=[u[0] for u in updates])
for panel, (_, values) in enumerate(updates):
    row, col = panel // 3 + 1, panel % 3 + 1
    fig.add_trace(go.Heatmap(z=values.numpy(), x=labels, y=labels,
                             colorscale="Magma", showscale=(panel == 4)), row=row, col=col)
fig.update_yaxes(autorange="reversed")
fig.update_layout(title="One miniature Pairformer-style block changes every ordered pair", height=790)
fig.show()

stage_names = ["block input", "+ outgoing", "+ incoming", "+ attention", "+ transition"]
stage_tensors = [trace.z_initial, trace.z_after_outgoing, trace.z_after_incoming,
                 trace.z_after_attention, trace.z_final]
cumulative_change = [
    float((stage[i, j] - trace.z_initial[i, j]).norm()) for stage in stage_tensors
]
fig = go.Figure(go.Bar(
    x=stage_names, y=cumulative_change, marker_color=["#cccccc", "#3978c5", "#df7b2c", "#6d55a3", "#4f9b75"],
    text=[f"{value:.3f}" for value in cumulative_change], textposition="outside",
))
fig.update_layout(
    title=f"Residual corrections accumulate for {labels[i]} → {labels[j]}",
    xaxis_title="sequential state inside one miniature block",
    yaxis_title="||current zᵢⱼ − block-input zᵢⱼ||", height=430,
)
fig.show()

# A matrix view makes “starting” versus “ending” precise without a crowded graph.
start_view = np.zeros((n_tokens, n_tokens))
end_view = np.zeros((n_tokens, n_tokens))
start_view[i, :] = 1                       # z[i,k]: hold the starting node i fixed
end_view[:, j] = 1                         # z[k,j]: hold the ending node j fixed
start_view[i, j] = end_view[i, j] = 2      # the query pair z[i,j]
attention_colors = [
    [0.00, "#f1f1f1"], [0.32, "#f1f1f1"],
    [0.34, "#719bd0"], [0.65, "#719bd0"],
    [0.67, "#df7b2c"], [1.00, "#df7b2c"],
]
fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=[
        f"Starting-node: fix i = {labels[i]}; read z[i,k]",
        f"Ending-node: fix j = {labels[j]}; read z[k,j]",
    ],
)
fig.add_trace(go.Heatmap(z=start_view, x=labels, y=labels, colorscale=attention_colors,
                         zmin=0, zmax=2, showscale=False,
                         hovertemplate="source=%{y}<br>destination=%{x}<extra></extra>"), row=1, col=1)
fig.add_trace(go.Heatmap(z=end_view, x=labels, y=labels, colorscale=attention_colors,
                         zmin=0, zmax=2, showscale=False,
                         hovertemplate="source=%{y}<br>destination=%{x}<extra></extra>"), row=1, col=2)
fig.update_yaxes(autorange="reversed", title="source token")
fig.update_xaxes(title="destination token")
fig.update_layout(
    title="Blue entries provide context; orange is the query pair z[i,j]",
    height=570, width=1100,
)
fig.show()

weights = trace.attention_weights[i, j].numpy()
fig = go.Figure(go.Bar(x=labels, y=weights,
                       marker_color=[COLORS["protein"]] * n_protein + [COLORS["ligand"]] * n_ligand))
fig.update_layout(title=f"Starting-node attention for target pair {labels[i]} → {labels[j]}",
                  xaxis_title="third token k", yaxis_title="attention weight αᵢⱼₖ", height=420)
fig.show()
print("Outgoing residual identity error:",
      float((trace.z_after_outgoing - (trace.z_initial + trace.delta_outgoing)).abs().max()))
print("Incoming residual identity error:",
      float((trace.z_after_incoming - (trace.z_after_outgoing + trace.delta_incoming)).abs().max()))
print("Attention softmax check (must equal 1):", weights.sum())
"""
    ),
    markdown(
        r"""
## 5. From pair features to a distogram and a pocket crop

A learned linear head converts each final pair vector into 64 logits. Softmax turns those logits into a distance distribution:

$$
p_{ij}(b)=\operatorname{softmax}_b(W_dz_{ij}),\qquad
\mathbb{E}[d_{ij}]=\sum_b p_{ij}(b)c_b.
$$

Here $c_b$ is the center of distance bin $b$. Nesso uses the predicted distances to identify a provisional pocket. The released configuration refines tokens within **22 Å**, capped at **256 tokens**, and constructs the affinity region using a tighter predicted **15 Å** neighborhood. Because the crop comes from predicted distances, an early pocket error can remove useful context from later stages.

The probabilities below come from a random teaching head. They demonstrate the conversion only; their distances are not predictions.
"""
    ),
    code(
        """
g = torch.Generator().manual_seed(19)
w_dist = torch.randn(128, 64, generator=g) / math.sqrt(128)
dist_logits = trace.z_final @ w_dist
expected_distance, normalized_entropy = expected_distance_and_entropy(dist_logits)
probability = dist_logits[i, j].softmax(dim=-1)
centers = torch.empty(64)
boundaries = torch.linspace(2, 22, 63)
centers[0], centers[-1] = 1.5, 24.5
centers[1:-1] = (boundaries[:-1] + boundaries[1:]) / 2

fig = go.Figure(go.Bar(x=centers.numpy(), y=probability.numpy(), marker_color="#6d55a3"))
fig.add_vline(x=float(expected_distance[i, j]), line_dash="dash", line_color="black",
              annotation_text=f"mean {expected_distance[i,j]:.1f} Å")
fig.update_layout(title=f"Toy 64-bin distogram for {labels[i]} → {labels[j]}",
                  xaxis_title="distance-bin center (Å)", yaxis_title="probability", height=430)
fig.show()

print(f"Expected distance: {expected_distance[i,j]:.2f} Å")
print(f"Normalized entropy: {normalized_entropy[i,j]:.3f} (0=concentrated, 1=uniform)")
"""
    ),
    markdown(
        r"""
## 6. Which pair entries reach the affinity head?

The affinity module does not average the whole protein–protein block. Its pooling mask keeps:

- protein $\rightarrow$ ligand pairs;
- ligand $\rightarrow$ protein pairs; and
- off-diagonal ligand $\rightarrow$ ligand pairs.

Those selected vectors are pooled into a global representation, processed by an affinity-specific network, and mapped to a scalar. The release stores two independently predicted continuous values, `pred_value1` and `pred_value2`; `pred_value` is their mean. It separately reports a binary binder probability. These are different outputs and must not be mixed.
"""
    ),
    code(
        """
pooled, affinity_mask = pool_affinity_pairs(trace.z_final, is_protein, is_ligand)
mask_labels = np.where(affinity_mask.numpy(), "included", "excluded")
fig = go.Figure(go.Heatmap(
    z=affinity_mask.int().numpy(), x=labels, y=labels,
    colorscale=[[0, "#eeeeee"], [1, "#9b6cc1"]], zmin=0, zmax=1, showscale=False,
    text=mask_labels, hovertemplate="i=%{y}<br>j=%{x}<br>%{text}<extra></extra>",
))
fig.add_vline(x=n_protein - 0.5, line_color="white", line_width=3)
fig.add_hline(y=n_protein - 0.5, line_color="white", line_width=3)
fig.update_layout(title="Purple pair entries are pooled for the toy affinity representation",
                  height=600, yaxis_autorange="reversed")
fig.show()

g = torch.Generator().manual_seed(23)
w1, w2 = torch.randn(128, generator=g) / 20, torch.randn(128, generator=g) / 20
toy_members = torch.stack([pooled @ w1, pooled @ w2])
toy_mean = toy_members.mean()
print("Selected pairs:", int(affinity_mask.sum()), "of", n_tokens * n_tokens)
print("Pooled vector shape:", tuple(pooled.shape))
print("Random teaching head values:", toy_members.tolist(), "mean:", toy_mean.item())
"""
    ),
    markdown(
        r"""
## 7. The same objects in one saved released-checkpoint prediction

We now leave the toy network. The following cell loads raw output from one completed Nesso-1 run in experiment `exp012`. This system was given a protein sequence and ligand chemistry; experimental coordinates were used later for scoring, not as Nesso input.

The saved tensor `z` is the released model's refined $256\times256\times128$ pair representation. `pdistogram` was saved on a padded $576\times576$ grid; `refine_mask` selects the same 256 refined tokens. The plots therefore align the refined distogram with `z` and `mol_type` before comparing protein–ligand entries.
"""
    ),
    code(
        """
REAL_DIR = REPO_ROOT / "runs/exp012_nesso1_rnp_distogram_generalization/nesso1/seed42/full100/raw/predictions/7fd6__1__1.A__1.B__1.B"
tensor_path = REAL_DIR / "predictions.safetensors"
affinity_path = REAL_DIR / "affinity.json"
if not tensor_path.exists():
    raise FileNotFoundError(f"Saved Nesso tensors are unavailable: {tensor_path}")

real = load_file(tensor_path)
affinity = json.loads(affinity_path.read_text())
refine_mask = real["refine_mask"].bool()
mol_type = real["mol_type"].long()
real_z = real["z"].float()
real_distogram = real["pdistogram"][refine_mask][:, refine_mask].float()
real_expected, real_entropy = expected_distance_and_entropy(real_distogram)
real_is_protein = mol_type == 0
real_is_ligand = mol_type == 3
n_real_protein = int(real_is_protein.sum())
n_real_ligand = int(real_is_ligand.sum())

print("System: 7fd6__1__1.A__1.B__1.B")
print("mol_type:", tuple(mol_type.shape), f"({n_real_protein} protein residues, {n_real_ligand} ligand atoms)")
print("z:", tuple(real_z.shape))
print("full padded pdistogram:", tuple(real["pdistogram"].shape))
print("refined aligned pdistogram:", tuple(real_distogram.shape))
print("refinement crop tokens:", int(refine_mask.sum()), "/", len(refine_mask))
print("22 Å pocket-mask tokens before final selection:", int(real["pocket_mask"].sum()))
"""
    ),
    code(
        """
boundary = n_real_protein - 0.5
fig = go.Figure(go.Heatmap(
    z=real_z.norm(dim=-1).numpy(), colorscale="Viridis", colorbar=dict(title="||zᵢⱼ||"),
    hovertemplate="i=%{y}<br>j=%{x}<br>norm=%{z:.2f}<extra></extra>",
))
fig.add_vline(x=boundary, line_color="white", line_width=2,
              annotation_text="protein | ligand")
fig.add_hline(y=boundary, line_color="white", line_width=2)
fig.update_layout(title="Released Nesso-1 refined pair tensor (actual saved checkpoint output)",
                  xaxis_title="token j", yaxis_title="token i", yaxis_autorange="reversed", height=680)
fig.show()

pl_expected = real_expected[real_is_protein][:, real_is_ligand]
pl_entropy = real_entropy[real_is_protein][:, real_is_ligand]
fig = make_subplots(rows=1, cols=2, subplot_titles=["Expected protein–ligand distance", "Normalized distogram entropy"])
fig.add_trace(go.Heatmap(z=pl_expected.numpy(), colorscale="Turbo", colorbar=dict(title="Å", x=0.44)), row=1, col=1)
fig.add_trace(go.Heatmap(z=pl_entropy.numpy(), colorscale="Magma", zmin=0, zmax=1,
                         colorbar=dict(title="entropy", x=1.01)), row=1, col=2)
fig.update_yaxes(autorange="reversed", title="protein residue index")
fig.update_xaxes(title="ligand atom index")
fig.update_layout(title=f"Actual protein–ligand distograms: {n_real_protein} residues × {n_real_ligand} atoms",
                  height=620, width=1050)
fig.show()
"""
    ),
    code(
        """
flat_index = int(pl_expected.argmin())
pi, li = np.unravel_index(flat_index, tuple(pl_expected.shape))
protein_indices = torch.where(real_is_protein)[0]
ligand_indices = torch.where(real_is_ligand)[0]
ri, rj = int(protein_indices[pi]), int(ligand_indices[li])
real_prob = real_distogram[ri, rj].softmax(dim=-1)

fig = go.Figure(go.Bar(x=centers.numpy(), y=real_prob.numpy(), marker_color="#3978c5"))
fig.add_vline(x=float(real_expected[ri, rj]), line_dash="dash", line_color="black",
              annotation_text=f"mean {real_expected[ri,rj]:.2f} Å")
fig.update_layout(title=f"Actual 64-bin distogram for closest predicted pair: protein token {ri}, ligand token {rj}",
                  xaxis_title="distance-bin center (Å)", yaxis_title="probability", height=430)
fig.show()

members = [affinity["affinity_pred_value1"], affinity["affinity_pred_value2"]]
mean_value = affinity["affinity_pred_value"]
fig = go.Figure(go.Bar(
    x=["member 1", "member 2", "reported mean"], y=members + [mean_value],
    marker_color=["#719bd0", "#9b6cc1", "#df7b2c"],
    text=[f"{v:.3f}" for v in members + [mean_value]], textposition="outside",
))
fig.update_layout(title="Actual saved Nesso-1 continuous affinity ensemble output",
                  yaxis_title="released continuous score (lower = stronger)", height=430)
fig.show()

print(f"Reported mean score: {mean_value:.6f} = ({members[0]:.6f} + {members[1]:.6f}) / 2")
print(f"Nominal concentration from the released convention: 10^score = {10**mean_value:.3f} µM")
print(f"Separate binary binder probability: {affinity['affinity_probability_binary']:.3f}")
"""
    ),
    markdown(
        r"""
## 8. What this notebook establishes—and what it does not

**We can now see:**

1. Nesso first creates its own 384-channel atom/residue input tokens; ESM-2 separately supplies 1280-channel protein sequence context that is projected to 384 channels.
2. $z_{ij}^{\rm init}$ is built from endpoint input tokens, relative position, and ligand bonds; an ESM-derived endpoint update is then added before Pairformer—not a simple dot product.
3. outgoing/incoming triangle multiplication and pair attention propagate evidence through third tokens $k$.
4. a 64-bin distogram supports expected-distance, uncertainty, and predicted-pocket calculations.
5. selected protein–ligand and ligand–ligand pair features feed an affinity-specific two-member ensemble.
6. the saved checkpoint tensors have the expected released dimensions: $z\in\mathbb{R}^{256\times256\times128}$ and a 64-bin distogram.

**This notebook does not establish:** that an attention weight is a physical interaction, that a low-entropy distance is correct, that the two ensemble values are experimental replicates, or that the continuous score is a universal equilibrium constant. Those require validation against held-out experimental measurements and structures.

The reusable teaching operations are in `src/affinity_benchmark/educational/mini_nesso.py`; focused numerical tests are in `tests/test_mini_nesso.py`.
"""
    ),
]

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(notebook, OUTPUT)
print(f"Wrote {OUTPUT}")
