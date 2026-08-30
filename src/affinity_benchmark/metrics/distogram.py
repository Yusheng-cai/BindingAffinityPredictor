"""Dependency-light metrics for categorical inter-token distance predictions."""

from __future__ import annotations

import numpy as np


def distogram_boundaries(
    num_bins: int = 64, min_distance: float = 2.0, max_distance: float = 22.0
) -> np.ndarray:
    """Return the ordered upper boundaries used by the Nesso-1 distogram."""

    if num_bins < 2:
        raise ValueError("num_bins must be at least 2")
    if min_distance >= max_distance:
        raise ValueError("min_distance must be smaller than max_distance")
    return np.linspace(min_distance, max_distance, num_bins - 1, dtype=float)


def distogram_centers(
    num_bins: int = 64, min_distance: float = 2.0, max_distance: float = 22.0
) -> np.ndarray:
    """Return Nesso-1's representative distance for every distogram bin."""

    boundaries = distogram_boundaries(num_bins, min_distance, max_distance)
    centers = np.empty(num_bins, dtype=float)
    centers[0] = 1.5
    centers[-1] = 24.5
    centers[1:-1] = (boundaries[:-1] + boundaries[1:]) / 2.0
    return centers


def distance_bin_indices(
    distances: np.ndarray,
    num_bins: int = 64,
    min_distance: float = 2.0,
    max_distance: float = 22.0,
) -> np.ndarray:
    """Map distances to categorical bins, including the two open end bins."""

    values = np.asarray(distances, dtype=float)
    boundaries = distogram_boundaries(num_bins, min_distance, max_distance)
    return np.searchsorted(boundaries, values, side="left")


def _log_softmax(logits: np.ndarray) -> np.ndarray:
    values = np.asarray(logits, dtype=float)
    maximum = np.max(values, axis=-1, keepdims=True)
    shifted = values - maximum
    return shifted - np.log(np.sum(np.exp(shifted), axis=-1, keepdims=True))


def distogram_probabilities(logits: np.ndarray) -> np.ndarray:
    """Convert distogram logits to normalized probabilities."""

    return np.exp(_log_softmax(logits))


def expected_distances(
    logits: np.ndarray,
    min_distance: float = 2.0,
    max_distance: float = 22.0,
) -> np.ndarray:
    """Calculate expected distances using the centers defined by Nesso-1."""

    values = np.asarray(logits)
    centers = distogram_centers(values.shape[-1], min_distance, max_distance)
    return np.sum(distogram_probabilities(values) * centers, axis=-1)


def distogram_negative_log_likelihood(
    logits: np.ndarray,
    distances: np.ndarray,
    mask: np.ndarray | None = None,
    min_distance: float = 2.0,
    max_distance: float = 22.0,
) -> float:
    """Mean negative log probability assigned to the observed distance bins."""

    values = np.asarray(logits)
    observed = np.asarray(distances, dtype=float)
    if values.shape[:-1] != observed.shape:
        raise ValueError("logits and distances have incompatible shapes")
    selected = np.isfinite(observed)
    if mask is not None:
        supplied = np.asarray(mask, dtype=bool)
        if supplied.shape != observed.shape:
            raise ValueError("mask and distances must have the same shape")
        selected &= supplied
    if not np.any(selected):
        raise ValueError("no finite distances were selected")

    bins = distance_bin_indices(
        observed, values.shape[-1], min_distance, max_distance
    )
    log_probs = _log_softmax(values)
    chosen = np.take_along_axis(log_probs, bins[..., None], axis=-1)[..., 0]
    return float(-np.mean(chosen[selected]))


def contact_probabilities(
    logits: np.ndarray,
    cutoff: float = 6.0,
    min_distance: float = 2.0,
    max_distance: float = 22.0,
) -> np.ndarray:
    """Approximate P(distance <= cutoff) by summing bins centered below cutoff."""

    values = np.asarray(logits)
    centers = distogram_centers(values.shape[-1], min_distance, max_distance)
    return np.sum(distogram_probabilities(values)[..., centers <= cutoff], axis=-1)


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    """Return ranking average precision for binary labels."""

    y = np.asarray(labels, dtype=bool).ravel()
    s = np.asarray(scores, dtype=float).ravel()
    if y.shape != s.shape:
        raise ValueError("labels and scores must have the same shape")
    positives = int(np.sum(y))
    if positives == 0:
        raise ValueError("average precision requires at least one positive")
    order = np.argsort(-s, kind="mergesort")
    ranked = y[order]
    precision = np.cumsum(ranked) / np.arange(1, len(ranked) + 1)
    return float(np.sum(precision[ranked]) / positives)


def binary_auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Return AUROC from all positive-negative score comparisons, ties = 0.5."""

    y = np.asarray(labels, dtype=bool).ravel()
    s = np.asarray(scores, dtype=float).ravel()
    if y.shape != s.shape:
        raise ValueError("labels and scores must have the same shape")
    positive = s[y]
    negative = s[~y]
    if len(positive) == 0 or len(negative) == 0:
        raise ValueError("AUROC requires both positive and negative labels")
    comparisons = positive[:, None] - negative[None, :]
    return float(np.mean((comparisons > 0) + 0.5 * (comparisons == 0)))


def precision_recall_f1(
    true_labels: np.ndarray, predicted_labels: np.ndarray
) -> tuple[float, float, float]:
    """Return binary precision, recall, and F1."""

    truth = np.asarray(true_labels, dtype=bool).ravel()
    pred = np.asarray(predicted_labels, dtype=bool).ravel()
    if truth.shape != pred.shape:
        raise ValueError("true and predicted labels must have the same shape")
    tp = int(np.sum(truth & pred))
    fp = int(np.sum(~truth & pred))
    fn = int(np.sum(truth & ~pred))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return float(precision), float(recall), float(f1)
