import torch

from affinity_benchmark.educational.mini_nesso import (
    affinity_pair_mask,
    expected_distance_and_entropy,
    initialize_pair_representation,
    one_mini_pairformer_step,
    pool_affinity_pairs,
    triangle_multiplication_incoming,
    triangle_multiplication_outgoing,
)


def test_initial_pair_representation_shape_and_seed():
    token = torch.arange(20, dtype=torch.float32).reshape(5, 4) / 20
    relative = torch.zeros(5, 5, 2)
    bond = torch.zeros(5, 5, 1)
    first = initialize_pair_representation(token, relative, bond, pair_dim=7, seed=3)
    second = initialize_pair_representation(token, relative, bond, pair_dim=7, seed=3)
    assert first.shape == (5, 5, 7)
    assert torch.equal(first, second)


def test_triangle_multiplication_matches_explicit_einsum():
    z = torch.arange(27, dtype=torch.float32).reshape(3, 3, 3) / 10
    identity = torch.eye(3)
    outgoing = triangle_multiplication_outgoing(z, identity, identity)
    incoming = triangle_multiplication_incoming(z, identity, identity)
    scale = 3**0.5
    assert torch.allclose(outgoing, torch.einsum("ikd,jkd->ijd", z, z) / scale)
    assert torch.allclose(incoming, torch.einsum("kid,kjd->ijd", z, z) / scale)


def test_attention_weights_are_probabilities_and_trace_shapes_match():
    z = torch.randn(6, 6, 8, generator=torch.Generator().manual_seed(2))
    trace = one_mini_pairformer_step(z, seed=4)
    assert trace.z_final.shape == z.shape
    assert trace.attention_weights.shape == (6, 6, 6)
    assert torch.allclose(
        trace.attention_weights.sum(dim=-1),
        torch.ones(6, 6),
        atol=1e-6,
    )


def test_expected_distance_and_entropy_for_certain_bin():
    logits = torch.full((2, 64), -30.0)
    logits[0, 0] = 30.0
    logits[1, -1] = 30.0
    expected, entropy = expected_distance_and_entropy(logits)
    assert torch.allclose(expected, torch.tensor([1.5, 24.5]), atol=1e-5)
    assert torch.all(entropy < 1e-5)


def test_affinity_pool_uses_pl_lp_and_off_diagonal_ll_pairs():
    z = torch.arange(4 * 4 * 2, dtype=torch.float32).reshape(4, 4, 2)
    protein = torch.tensor([True, True, False, False])
    ligand = ~protein
    pooled, mask = pool_affinity_pairs(z, protein, ligand)
    expected_mask = affinity_pair_mask(protein, ligand)
    assert torch.equal(mask, expected_mask)
    assert mask.sum().item() == 10
    assert torch.allclose(pooled, z[mask].mean(dim=0))
