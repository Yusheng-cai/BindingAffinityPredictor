"""Minimal PyTorch models illustrating RFdiffusion and ProteinMPNN ideas.

The design goal is inspectability, not biological realism.  The diffusion
model operates on C-alpha coordinates rather than full residue frames.  The
sequence model uses six coarse chemical classes rather than twenty amino
acids.  Both choices are explicit simplifications for teaching.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import torch
from torch import Tensor, nn


ALPHABET = ("H", "D", "A", "+", "-", "N")
TOKEN_NAMES = (
    "hydrophobic",
    "donor",
    "acceptor",
    "positive",
    "negative",
    "neutral",
)
UNKNOWN_TOKEN = len(ALPHABET)


def make_helix(
    length: int,
    radius: float = 2.3,
    rise: float = 1.5,
    degrees_per_residue: float = 100.0,
    device: Optional[torch.device] = None,
) -> Tensor:
    """Return an idealized C-alpha helix centered at its coordinate mean."""

    index = torch.arange(length, dtype=torch.float32, device=device)
    angle = index * math.radians(degrees_per_residue)
    xyz = torch.stack(
        (radius * torch.cos(angle), radius * torch.sin(angle), rise * index),
        dim=-1,
    )
    return xyz - xyz.mean(dim=0, keepdim=True)


def random_rotation(
    batch_size: int,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float32,
) -> Tensor:
    """Sample proper 3D rotation matrices using QR decomposition."""

    matrix = torch.randn(batch_size, 3, 3, device=device, dtype=dtype)
    q, r = torch.linalg.qr(matrix)
    signs = torch.sign(torch.diagonal(r, dim1=-2, dim2=-1))
    signs = torch.where(signs == 0, torch.ones_like(signs), signs)
    q = q * signs[:, None, :]
    determinant = torch.linalg.det(q)
    q[:, :, -1] *= torch.where(determinant < 0, -1.0, 1.0)[:, None]
    return q


def apply_rotation(coordinates: Tensor, rotation: Tensor) -> Tensor:
    """Apply batched rotation matrices to row-vector coordinates."""

    return torch.einsum("bij,bnj->bni", rotation, coordinates)


@dataclass
class SyntheticComplex:
    """A small, procedural target-binder complex and coarse chemistry labels."""

    target: Tensor
    binder: Tensor
    target_tokens: Tensor
    binder_tokens: Tensor


def make_synthetic_complex(
    batch_size: int = 1,
    target_length: int = 14,
    binder_length: int = 11,
    rotate: bool = True,
    coordinate_noise: float = 0.0,
    device: Optional[torch.device] = None,
) -> SyntheticComplex:
    """Construct parallel idealized helices with complementary class labels.

    The complex is deliberately procedural.  Binder labels are assigned as the
    complement of the nearest target label, so the rule is known exactly.
    """

    target = make_helix(target_length, device=device)
    binder = make_helix(binder_length, radius=2.1, rise=1.45, device=device)
    binder = binder + torch.tensor([6.0, 0.0, 0.5], device=device)

    target = target.unsqueeze(0).expand(batch_size, -1, -1).clone()
    binder = binder.unsqueeze(0).expand(batch_size, -1, -1).clone()
    if coordinate_noise:
        target = target + coordinate_noise * torch.randn_like(target)
        binder = binder + coordinate_noise * torch.randn_like(binder)

    if rotate:
        rotation = random_rotation(batch_size, device=device, dtype=target.dtype)
        target = apply_rotation(target, rotation)
        binder = apply_rotation(binder, rotation)

    target_tokens = (
        torch.arange(target_length, device=device) % len(ALPHABET)
    ).unsqueeze(0).expand(batch_size, -1).clone()
    distances = torch.cdist(binder, target)
    nearest_target = distances.argmin(dim=-1)
    nearest_tokens = torch.gather(target_tokens, 1, nearest_target)
    complement = torch.tensor([0, 2, 1, 4, 3, 5], device=device)
    binder_tokens = complement[nearest_tokens]
    return SyntheticComplex(target, binder, target_tokens, binder_tokens)


def cosine_schedule(
    steps: int,
    s: float = 0.008,
    device: Optional[torch.device] = None,
) -> Dict[str, Tensor]:
    """Construct a cosine DDPM schedule with steps indexed from 0 to T-1."""

    grid = torch.linspace(0, steps, steps + 1, device=device)
    alpha_bar = torch.cos(((grid / steps + s) / (1 + s)) * math.pi / 2).square()
    alpha_bar = alpha_bar / alpha_bar[0]
    betas = 1 - alpha_bar[1:] / alpha_bar[:-1]
    betas = betas.clamp(1e-5, 0.999)
    alphas = 1 - betas
    alpha_bar = torch.cumprod(alphas, dim=0)
    return {"betas": betas, "alphas": alphas, "alpha_bar": alpha_bar}


def q_sample(x0: Tensor, timestep: Tensor, alpha_bar: Tensor, noise: Optional[Tensor] = None) -> Tuple[Tensor, Tensor]:
    """Sample x_t from q(x_t | x_0)."""

    if noise is None:
        noise = torch.randn_like(x0)
    a = alpha_bar[timestep].view(-1, 1, 1)
    xt = a.sqrt() * x0 + (1 - a).sqrt() * noise
    return xt, noise


def _rbf(distance: Tensor, bins: int = 16, maximum: float = 20.0) -> Tensor:
    centers = torch.linspace(0.0, maximum, bins, device=distance.device, dtype=distance.dtype)
    width = maximum / bins
    return torch.exp(-((distance[..., None] - centers) / width).square())


def _index_features(length: int, batch_size: int, device: torch.device, dtype: torch.dtype) -> Tensor:
    if length == 1:
        fraction = torch.zeros(1, device=device, dtype=dtype)
    else:
        fraction = torch.linspace(0, 1, length, device=device, dtype=dtype)
    features = torch.stack(
        (fraction, torch.sin(math.pi * fraction), torch.cos(math.pi * fraction)),
        dim=-1,
    )
    return features.unsqueeze(0).expand(batch_size, -1, -1)


class EquivariantCoordinateBlock(nn.Module):
    """Scalar message passing plus an equivariant relative-vector update."""

    def __init__(self, hidden_dim: int, rbf_bins: int = 16) -> None:
        super().__init__()
        self.rbf_bins = rbf_bins
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim + rbf_bins, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.coordinate_mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, 1),
            nn.Tanh(),
        )

    def forward(self, coordinates: Tensor, node_features: Tensor, movable: Tensor) -> Tuple[Tensor, Tensor]:
        batch_size, nodes, _ = coordinates.shape
        relative = coordinates[:, :, None, :] - coordinates[:, None, :, :]
        distance = relative.square().sum(dim=-1).add(1e-8).sqrt()
        radial = _rbf(distance, self.rbf_bins)
        hi = node_features[:, :, None, :].expand(-1, -1, nodes, -1)
        hj = node_features[:, None, :, :].expand(-1, nodes, -1, -1)
        messages = self.edge_mlp(torch.cat((hi, hj, radial), dim=-1))

        pair_mask = 1 - torch.eye(nodes, device=coordinates.device, dtype=coordinates.dtype)
        pair_mask = pair_mask.view(1, nodes, nodes, 1)
        messages = messages * pair_mask
        aggregated = messages.sum(dim=2) / max(nodes - 1, 1)
        node_features = node_features + self.node_mlp(torch.cat((node_features, aggregated), dim=-1))

        scalar = self.coordinate_mlp(messages) * pair_mask
        direction = relative / distance[..., None]
        # A square-root normalization keeps updates usable in this tiny model
        # while preventing their scale from growing linearly with graph size.
        delta = (direction * scalar).sum(dim=2) / math.sqrt(max(nodes - 1, 1))
        coordinates = coordinates + delta * movable[..., None]
        return coordinates, node_features


class MiniEquivariantDenoiser(nn.Module):
    """C-alpha-only, target-conditioned, equivariant x0 predictor."""

    def __init__(self, hidden_dim: int = 64, layers: int = 4, max_steps: int = 50) -> None:
        super().__init__()
        self.max_steps = max_steps
        self.input_projection = nn.Sequential(nn.Linear(7, hidden_dim), nn.SiLU())
        self.blocks = nn.ModuleList(
            EquivariantCoordinateBlock(hidden_dim) for _ in range(layers)
        )

    def forward(self, noisy_binder: Tensor, target: Tensor, timestep: Tensor) -> Tensor:
        batch_size, binder_length, _ = noisy_binder.shape
        target_length = target.shape[1]
        coordinates = torch.cat((noisy_binder, target), dim=1)

        binder_type = torch.tensor([1.0, 0.0], device=coordinates.device, dtype=coordinates.dtype)
        target_type = torch.tensor([0.0, 1.0], device=coordinates.device, dtype=coordinates.dtype)
        types = torch.cat(
            (
                binder_type.view(1, 1, 2).expand(batch_size, binder_length, 2),
                target_type.view(1, 1, 2).expand(batch_size, target_length, 2),
            ),
            dim=1,
        )
        indices = torch.cat(
            (
                _index_features(binder_length, batch_size, coordinates.device, coordinates.dtype),
                _index_features(target_length, batch_size, coordinates.device, coordinates.dtype),
            ),
            dim=1,
        )
        time = timestep.to(coordinates.dtype) / max(self.max_steps - 1, 1)
        time_features = torch.stack((torch.sin(math.pi * time), torch.cos(math.pi * time)), dim=-1)
        time_features = time_features[:, None, :].expand(-1, binder_length + target_length, -1)
        node_features = self.input_projection(torch.cat((types, indices, time_features), dim=-1))

        movable = torch.cat(
            (
                torch.ones(batch_size, binder_length, device=coordinates.device, dtype=coordinates.dtype),
                torch.zeros(batch_size, target_length, device=coordinates.device, dtype=coordinates.dtype),
            ),
            dim=1,
        )
        for block in self.blocks:
            coordinates, node_features = block(coordinates, node_features, movable)
        return coordinates[:, :binder_length]


@torch.no_grad()
def sample_ddim(
    model: MiniEquivariantDenoiser,
    target: Tensor,
    binder_length: int,
    schedule: Dict[str, Tensor],
    initial_noise: Optional[Tensor] = None,
    start_step: Optional[int] = None,
) -> Tuple[Tensor, List[Tensor]]:
    """Deterministic DDIM-style sampling for visualizing a reverse trajectory."""

    alpha_bar = schedule["alpha_bar"]
    batch_size = target.shape[0]
    if initial_noise is None:
        current = torch.randn(batch_size, binder_length, 3, device=target.device)
    else:
        current = initial_noise.clone()
    if start_step is None:
        start_step = len(alpha_bar) - 1
    if not 0 <= start_step < len(alpha_bar):
        raise ValueError(f"start_step must be in [0, {len(alpha_bar) - 1}]")
    trajectory = [current.detach().cpu()]
    model.eval()
    for step in reversed(range(start_step + 1)):
        timestep = torch.full((batch_size,), step, device=target.device, dtype=torch.long)
        predicted_x0 = model(current, target, timestep)
        if step == 0:
            current = predicted_x0
        else:
            a_t = alpha_bar[step]
            a_previous = alpha_bar[step - 1]
            predicted_noise = (current - a_t.sqrt() * predicted_x0) / (1 - a_t).sqrt().clamp_min(1e-6)
            current = a_previous.sqrt() * predicted_x0 + (1 - a_previous).sqrt() * predicted_noise
        trajectory.append(current.detach().cpu())
    return current, trajectory


class ScalarMessageBlock(nn.Module):
    """Invariant message-passing block for the miniature sequence model."""

    def __init__(self, hidden_dim: int, rbf_bins: int = 16) -> None:
        super().__init__()
        self.rbf_bins = rbf_bins
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim + rbf_bins + 1, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(self, coordinates: Tensor, node_features: Tensor, chain_id: Tensor) -> Tensor:
        nodes = coordinates.shape[1]
        distance = torch.cdist(coordinates, coordinates)
        radial = _rbf(distance, self.rbf_bins)
        same_chain = (chain_id[:, :, None] == chain_id[:, None, :]).to(coordinates.dtype)[..., None]
        hi = node_features[:, :, None, :].expand(-1, -1, nodes, -1)
        hj = node_features[:, None, :, :].expand(-1, nodes, -1, -1)
        messages = self.edge_mlp(torch.cat((hi, hj, radial, same_chain), dim=-1))
        pair_mask = 1 - torch.eye(nodes, device=coordinates.device, dtype=coordinates.dtype)
        messages = messages * pair_mask.view(1, nodes, nodes, 1)
        aggregated = messages.sum(dim=2) / max(nodes - 1, 1)
        return node_features + self.node_mlp(torch.cat((node_features, aggregated), dim=-1))


class MiniProteinMPNN(nn.Module):
    """Random-order masked sequence decoder on a fixed residue graph.

    This intentionally compact decoder is ProteinMPNN-like rather than a
    line-for-line reproduction.  Previously decoded tokens are supplied as
    context, and unknown binder positions use a dedicated mask token.
    """

    def __init__(self, hidden_dim: int = 64, layers: int = 3, classes: int = len(ALPHABET)) -> None:
        super().__init__()
        self.classes = classes
        self.geometry_projection = nn.Sequential(nn.Linear(5, hidden_dim), nn.SiLU())
        self.sequence_embedding = nn.Embedding(classes + 1, hidden_dim)
        self.blocks = nn.ModuleList(ScalarMessageBlock(hidden_dim) for _ in range(layers))
        self.output = nn.Sequential(nn.LayerNorm(hidden_dim), nn.Linear(hidden_dim, classes))

    def forward(
        self,
        binder: Tensor,
        target: Tensor,
        target_tokens: Tensor,
        partial_binder_tokens: Tensor,
    ) -> Tensor:
        batch_size, binder_length, _ = binder.shape
        target_length = target.shape[1]
        coordinates = torch.cat((binder, target), dim=1)
        chain_id = torch.cat(
            (
                torch.ones(batch_size, binder_length, device=binder.device, dtype=torch.long),
                torch.zeros(batch_size, target_length, device=binder.device, dtype=torch.long),
            ),
            dim=1,
        )
        types = torch.nn.functional.one_hot(chain_id, num_classes=2).to(binder.dtype)
        indices = torch.cat(
            (
                _index_features(binder_length, batch_size, binder.device, binder.dtype),
                _index_features(target_length, batch_size, binder.device, binder.dtype),
            ),
            dim=1,
        )
        tokens = torch.cat((partial_binder_tokens, target_tokens), dim=1)
        node_features = self.geometry_projection(torch.cat((types, indices), dim=-1))
        node_features = node_features + self.sequence_embedding(tokens)
        for block in self.blocks:
            node_features = block(coordinates, node_features, chain_id)
        return self.output(node_features[:, :binder_length])

    @torch.no_grad()
    def sample(
        self,
        binder: Tensor,
        target: Tensor,
        target_tokens: Tensor,
        temperature: float = 1.0,
        order: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Tensor]:
        """Sample binder classes one position at a time in an arbitrary order."""

        batch_size, binder_length, _ = binder.shape
        if order is None:
            order = torch.stack(
                [torch.randperm(binder_length, device=binder.device) for _ in range(batch_size)]
            )
        sequence = torch.full(
            (batch_size, binder_length), UNKNOWN_TOKEN, device=binder.device, dtype=torch.long
        )
        log_probability = torch.zeros(batch_size, device=binder.device)
        self.eval()
        for rank in range(binder_length):
            logits = self(binder, target, target_tokens, sequence) / temperature
            for batch in range(batch_size):
                position = order[batch, rank]
                distribution = torch.distributions.Categorical(logits=logits[batch, position])
                token = distribution.sample()
                sequence[batch, position] = token
                log_probability[batch] += distribution.log_prob(token)
        return sequence, log_probability


def random_partial_sequence(labels: Tensor) -> Tuple[Tensor, Tensor]:
    """Reveal a random prefix of a random decoding order for training."""

    batch_size, length = labels.shape
    partial = torch.full_like(labels, UNKNOWN_TOKEN)
    predict_mask = torch.ones_like(labels, dtype=torch.bool)
    for batch in range(batch_size):
        order = torch.randperm(length, device=labels.device)
        revealed = torch.randint(0, length, (1,), device=labels.device).item()
        known_positions = order[:revealed]
        partial[batch, known_positions] = labels[batch, known_positions]
        predict_mask[batch, known_positions] = False
    return partial, predict_mask


def compatibility_matrix(device: Optional[torch.device] = None) -> Tensor:
    """Return the explicit toy binder-row/target-column chemistry table."""

    matrix = torch.full((len(ALPHABET), len(ALPHABET)), -0.25, device=device)
    matrix[0, 0] = 1.0  # hydrophobic burial
    matrix[1, 2] = 1.0  # donor -> acceptor
    matrix[2, 1] = 1.0  # acceptor -> donor
    matrix[3, 4] = 1.0  # positive -> negative
    matrix[4, 3] = 1.0  # negative -> positive
    matrix[5, 5] = 0.2
    matrix[1, 1] = matrix[2, 2] = -0.6
    matrix[3, 3] = matrix[4, 4] = -1.0
    return matrix


def interface_score(
    binder: Tensor,
    target: Tensor,
    binder_tokens: Tensor,
    target_tokens: Tensor,
    contact_cutoff: float = 7.0,
    contact_width: float = 0.8,
) -> Tensor:
    """Evaluate a transparent toy interface-quality score.

    This score is not a binding free energy.  It combines smooth pair contacts,
    a steric-overlap penalty, and an exposed-hydrophobe penalty.
    """

    distances = torch.cdist(binder, target)
    contacts = torch.sigmoid((contact_cutoff - distances) / contact_width)
    matrix = compatibility_matrix(binder.device)
    pair_compatibility = matrix[binder_tokens[:, :, None], target_tokens[:, None, :]]
    chemical = (contacts * pair_compatibility).sum(dim=(1, 2))

    clash = torch.relu(2.8 - distances).square().sum(dim=(1, 2))
    contact_fraction = contacts.max(dim=-1).values
    exposed_hydrophobe = ((binder_tokens == 0).to(binder.dtype) * (1 - contact_fraction)).sum(dim=1)
    return chemical - 2.0 * clash - 0.5 * exposed_hydrophobe


@torch.no_grad()
def sample_with_score_guidance(
    model: MiniProteinMPNN,
    binder: Tensor,
    target: Tensor,
    target_tokens: Tensor,
    beta: float = 0.0,
    temperature: float = 1.0,
    order: Optional[Tensor] = None,
) -> Tuple[Tensor, Tensor]:
    """Bias sequence logits using changes in the transparent interface score.

    Unknown positions are provisionally treated as neutral while evaluating
    each candidate token. This is an intentionally simple illustration of
    logit guidance, not a production sequence-design algorithm.
    """

    batch_size, binder_length, _ = binder.shape
    if order is None:
        order = torch.stack(
            [torch.randperm(binder_length, device=binder.device) for _ in range(batch_size)]
        )
    sequence = torch.full(
        (batch_size, binder_length), UNKNOWN_TOKEN, device=binder.device, dtype=torch.long
    )
    log_probability = torch.zeros(batch_size, device=binder.device)
    model.eval()
    for rank in range(binder_length):
        logits = model(binder, target, target_tokens, sequence) / temperature
        for batch in range(batch_size):
            position = order[batch, rank]
            provisional = sequence[batch : batch + 1].clone()
            provisional[provisional == UNKNOWN_TOKEN] = 5
            candidate_scores = []
            for token in range(model.classes):
                candidate = provisional.clone()
                candidate[0, position] = token
                candidate_scores.append(
                    interface_score(
                        binder[batch : batch + 1],
                        target[batch : batch + 1],
                        candidate,
                        target_tokens[batch : batch + 1],
                    )[0]
                )
            guidance = torch.stack(candidate_scores)
            guidance = guidance - guidance.mean()
            guided_logits = logits[batch, position] + beta * guidance
            distribution = torch.distributions.Categorical(logits=guided_logits)
            token = distribution.sample()
            sequence[batch, position] = token
            log_probability[batch] += torch.log_softmax(logits[batch, position], dim=-1)[token]
    return sequence, log_probability


def tokens_to_string(tokens: Tensor) -> List[str]:
    """Convert a batch of coarse token sequences to compact strings."""

    return ["".join(ALPHABET[token] for token in row.tolist()) for row in tokens]
