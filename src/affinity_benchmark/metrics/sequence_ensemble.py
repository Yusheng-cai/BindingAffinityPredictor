"""Metrics for equal-length, fixed-backbone sequence ensembles."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def sequence_array(
    sequences: Sequence[str], alphabet: str = STANDARD_AMINO_ACIDS
) -> np.ndarray:
    """Return an N x L character array after validating an aligned ensemble."""
    if not sequences:
        raise ValueError("At least one sequence is required")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError(f"Sequences must have one common length, found {lengths}")
    unknown = sorted(set("".join(sequences)) - set(alphabet))
    if unknown:
        raise ValueError(f"Sequences contain residues outside the alphabet: {unknown}")
    return np.asarray([list(sequence) for sequence in sequences], dtype="U1")


def pairwise_identity_matrix(sequences: Sequence[str]) -> np.ndarray:
    """Calculate ungapped position-wise identity for equal-length sequences."""
    array = sequence_array(sequences)
    return np.mean(array[:, None, :] == array[None, :, :], axis=2)


def position_frequencies(
    sequences: Sequence[str], alphabet: str = STANDARD_AMINO_ACIDS
) -> np.ndarray:
    """Return empirical amino-acid frequencies with shape L x len(alphabet)."""
    array = sequence_array(sequences, alphabet=alphabet)
    return np.stack([np.mean(array == residue, axis=0) for residue in alphabet], axis=1)


def shannon_entropy(frequencies: np.ndarray) -> np.ndarray:
    """Calculate plug-in Shannon entropy in nats along the last axis."""
    frequencies = np.asarray(frequencies, dtype=float)
    if frequencies.ndim < 1:
        raise ValueError("Frequencies must have at least one dimension")
    if np.any(frequencies < 0):
        raise ValueError("Frequencies cannot be negative")
    totals = frequencies.sum(axis=-1)
    if not np.allclose(totals, 1.0):
        raise ValueError("Frequencies must sum to one along the last axis")
    terms = np.zeros_like(frequencies)
    positive = frequencies > 0
    terms[positive] = frequencies[positive] * np.log(frequencies[positive])
    return -terms.sum(axis=-1)
