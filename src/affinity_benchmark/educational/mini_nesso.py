"""Small, transparent operations for teaching the Nesso-1 pair pathway.

These functions mirror the *information flow* of selected released Nesso-1
operations while remaining deliberately small and inspectable.  They do not
contain Nesso weights and must not be used for scientific prediction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import torch
from torch import Tensor


@dataclass(frozen=True)
class MiniPairTrace:
    """Intermediate tensors from one educational pair-update block."""

    z_initial: Tensor
    delta_outgoing: Tensor
    z_after_outgoing: Tensor
    delta_incoming: Tensor
    z_after_incoming: Tensor
    delta_attention: Tensor
    attention_weights: Tensor
    z_after_attention: Tensor
    delta_transition: Tensor
    z_final: Tensor


def _random_matrix(
    rows: int,
    columns: int,
    *,
    generator: torch.Generator,
    dtype: torch.dtype,
    device: torch.device,
    scale: float | None = None,
) -> Tensor:
    """Create a deterministic teaching weight without changing global RNG state."""

    if scale is None:
        scale = 1.0 / math.sqrt(max(rows, 1))
    return torch.randn(
        rows,
        columns,
        generator=generator,
        dtype=dtype,
        device=device,
    ) * scale


def initialize_pair_representation(
    token_features: Tensor,
    relative_features: Tensor,
    bond_features: Tensor,
    *,
    pair_dim: int = 128,
    seed: int = 0,
) -> Tensor:
    """Construct a Nesso-inspired initial pair tensor from transparent inputs.

    Parameters
    ----------
    token_features
        ``[N, C_s]`` combined token features.
    relative_features
        ``[N, N, C_r]`` sequence/chain relationship features.
    bond_features
        ``[N, N, C_b]`` ligand bond features.
    pair_dim
        Width of each pair entry.
    seed
        Deterministic seed for educational projection weights.
    """

    if token_features.ndim != 2:
        raise ValueError("token_features must have shape [N, C_s]")
    n_tokens, token_dim = token_features.shape
    if relative_features.shape[:2] != (n_tokens, n_tokens):
        raise ValueError("relative_features must have shape [N, N, C_r]")
    if bond_features.shape[:2] != (n_tokens, n_tokens):
        raise ValueError("bond_features must have shape [N, N, C_b]")

    generator = torch.Generator(device=token_features.device).manual_seed(seed)
    kwargs = {
        "generator": generator,
        "dtype": token_features.dtype,
        "device": token_features.device,
    }
    w_left = _random_matrix(token_dim, pair_dim, **kwargs)
    w_right = _random_matrix(token_dim, pair_dim, **kwargs)
    w_relative = _random_matrix(relative_features.shape[-1], pair_dim, **kwargs)
    w_bond = _random_matrix(bond_features.shape[-1], pair_dim, **kwargs)

    left = token_features @ w_left
    right = token_features @ w_right
    return (
        left[:, None, :]
        + right[None, :, :]
        + relative_features @ w_relative
        + bond_features @ w_bond
    )


def triangle_multiplication_outgoing(z: Tensor, a_weight: Tensor, b_weight: Tensor) -> Tensor:
    """Educational outgoing update ``sum_k A(z_ik) * B(z_jk)``."""

    a = torch.einsum("ijd,df->ijf", z, a_weight)
    b = torch.einsum("ijd,df->ijf", z, b_weight)
    return torch.einsum("ikd,jkd->ijd", a, b) / math.sqrt(z.shape[0])


def triangle_multiplication_incoming(z: Tensor, a_weight: Tensor, b_weight: Tensor) -> Tensor:
    """Educational incoming update ``sum_k A(z_ki) * B(z_kj)``."""

    a = torch.einsum("ijd,df->ijf", z, a_weight)
    b = torch.einsum("ijd,df->ijf", z, b_weight)
    return torch.einsum("kid,kjd->ijd", a, b) / math.sqrt(z.shape[0])


def starting_node_attention(
    z: Tensor,
    q_weight: Tensor,
    k_weight: Tensor,
    v_weight: Tensor,
) -> tuple[Tensor, Tensor]:
    """One-head educational attention over third tokens while fixing node ``i``.

    Returns an update ``[N, N, D]`` and weights ``[N, N, N]`` indexed by
    ``[i, j, k]``.
    """

    q = torch.einsum("ijd,df->ijf", z, q_weight)
    k = torch.einsum("ijd,df->ijf", z, k_weight)
    v = torch.einsum("ijd,df->ijf", z, v_weight)
    scores = torch.einsum("ijd,ikd->ijk", q, k) / math.sqrt(q.shape[-1])
    weights = scores.softmax(dim=-1)
    update = torch.einsum("ijk,ikd->ijd", weights, v)
    return update, weights


def one_mini_pairformer_step(z: Tensor, *, seed: int = 1) -> MiniPairTrace:
    """Apply outgoing, incoming, attention, and transition residual updates."""

    if z.ndim != 3 or z.shape[0] != z.shape[1]:
        raise ValueError("z must have shape [N, N, D]")
    pair_dim = z.shape[-1]
    generator = torch.Generator(device=z.device).manual_seed(seed)
    kwargs = {
        "generator": generator,
        "dtype": z.dtype,
        "device": z.device,
        "scale": 0.25 / math.sqrt(pair_dim),
    }

    w_out_a = _random_matrix(pair_dim, pair_dim, **kwargs)
    w_out_b = _random_matrix(pair_dim, pair_dim, **kwargs)
    delta_out = triangle_multiplication_outgoing(z, w_out_a, w_out_b)
    z_out = z + delta_out

    w_in_a = _random_matrix(pair_dim, pair_dim, **kwargs)
    w_in_b = _random_matrix(pair_dim, pair_dim, **kwargs)
    delta_in = triangle_multiplication_incoming(z_out, w_in_a, w_in_b)
    z_in = z_out + delta_in

    w_q = _random_matrix(pair_dim, pair_dim, **kwargs)
    w_k = _random_matrix(pair_dim, pair_dim, **kwargs)
    w_v = _random_matrix(pair_dim, pair_dim, **kwargs)
    delta_attention, attention = starting_node_attention(z_in, w_q, w_k, w_v)
    delta_attention = 0.25 * delta_attention
    z_attention = z_in + delta_attention

    w_hidden = _random_matrix(pair_dim, 2 * pair_dim, **kwargs)
    w_final = _random_matrix(2 * pair_dim, pair_dim, **kwargs)
    delta_transition = torch.relu(z_attention @ w_hidden) @ w_final
    z_final = z_attention + delta_transition

    return MiniPairTrace(
        z_initial=z,
        delta_outgoing=delta_out,
        z_after_outgoing=z_out,
        delta_incoming=delta_in,
        z_after_incoming=z_in,
        delta_attention=delta_attention,
        attention_weights=attention,
        z_after_attention=z_attention,
        delta_transition=delta_transition,
        z_final=z_final,
    )


def distogram_bin_centers(
    num_bins: int = 64,
    min_dist: float = 2.0,
    max_dist: float = 22.0,
    *,
    device: torch.device | None = None,
) -> Tensor:
    """Return the bin centers used by the released Nesso inference utility."""

    boundaries = torch.linspace(min_dist, max_dist, num_bins - 1, device=device)
    centers = torch.empty(num_bins, device=device)
    centers[0] = 1.5
    centers[-1] = 24.5
    centers[1:-1] = (boundaries[:-1] + boundaries[1:]) * 0.5
    return centers


def expected_distance_and_entropy(logits: Tensor) -> tuple[Tensor, Tensor]:
    """Convert distogram logits to expected distance and normalized entropy."""

    probabilities = logits.float().softmax(dim=-1)
    centers = distogram_bin_centers(logits.shape[-1], device=logits.device)
    expected = torch.einsum("...b,b->...", probabilities, centers)
    entropy = -(probabilities * probabilities.clamp_min(1e-12).log()).sum(dim=-1)
    entropy = entropy / math.log(logits.shape[-1])
    return expected, entropy


def affinity_pair_mask(is_protein: Tensor, is_ligand: Tensor) -> Tensor:
    """Return the PL/LP/LL off-diagonal mask used for educational pooling."""

    is_protein = is_protein.bool()
    is_ligand = is_ligand.bool()
    mask = (
        is_ligand[:, None] * is_protein[None, :]
        + is_protein[:, None] * is_ligand[None, :]
        + is_ligand[:, None] * is_ligand[None, :]
    ).bool()
    mask.fill_diagonal_(False)
    return mask


def pool_affinity_pairs(z: Tensor, is_protein: Tensor, is_ligand: Tensor) -> tuple[Tensor, Tensor]:
    """Average interface and ligand-pair entries into one global vector."""

    mask = affinity_pair_mask(is_protein, is_ligand)
    if not mask.any():
        raise ValueError("affinity mask is empty")
    return z[mask].mean(dim=0), mask
