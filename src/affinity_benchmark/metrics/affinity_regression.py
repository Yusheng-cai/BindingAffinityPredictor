"""Regression metrics for within-assay protein--ligand affinity benchmarks."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy.stats import kendalltau, pearsonr, spearmanr


R_KCAL_PER_MOL_K = 0.00198720425864083
DEFAULT_TEMPERATURE_K = 298.15


def log10_micromolar_to_kcal_per_mol(
    value: np.ndarray | Iterable[float] | float,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
) -> np.ndarray:
    """Convert log10(concentration / micromolar) to RT ln(concentration / M)."""

    values = np.asarray(value, dtype=float)
    return R_KCAL_PER_MOL_K * temperature_k * np.log(10.0) * (values - 6.0)


def pairwise_mae(observed: Iterable[float], predicted: Iterable[float]) -> float:
    """Return MAE over all unordered pairwise affinity differences."""

    observed_array = np.asarray(list(observed), dtype=float)
    predicted_array = np.asarray(list(predicted), dtype=float)
    if observed_array.shape != predicted_array.shape or observed_array.size < 2:
        raise ValueError("observed and predicted must have one common length >= 2")
    upper = np.triu_indices(observed_array.size, k=1)
    observed_delta = observed_array[:, None] - observed_array[None, :]
    predicted_delta = predicted_array[:, None] - predicted_array[None, :]
    return float(np.mean(np.abs(observed_delta[upper] - predicted_delta[upper])))


def assay_metrics(observed: Iterable[float], predicted: Iterable[float]) -> dict[str, float]:
    """Calculate the paper-aligned metrics for one assay or target."""

    observed_array = np.asarray(list(observed), dtype=float)
    predicted_array = np.asarray(list(predicted), dtype=float)
    if observed_array.shape != predicted_array.shape or observed_array.size < 2:
        raise ValueError("observed and predicted must have one common length >= 2")
    centered = predicted_array - np.mean(predicted_array) + np.mean(observed_array)
    return {
        "n": int(observed_array.size),
        "pearson_r": float(pearsonr(observed_array, predicted_array).statistic),
        "spearman_rho": float(spearmanr(observed_array, predicted_array).statistic),
        "kendall_tau": float(kendalltau(observed_array, predicted_array).statistic),
        "pairwise_mae_kcal_mol": pairwise_mae(observed_array, predicted_array),
        "mae_kcal_mol": float(np.mean(np.abs(predicted_array - observed_array))),
        "centered_mae_kcal_mol": float(np.mean(np.abs(centered - observed_array))),
    }


def compound_weighted_average(
    metrics_by_target: dict[str, dict[str, float]], metric: str
) -> float:
    """Average a per-target metric with weights proportional to compound count."""

    weights = np.asarray([result["n"] for result in metrics_by_target.values()], dtype=float)
    values = np.asarray([result[metric] for result in metrics_by_target.values()], dtype=float)
    return float(np.average(values, weights=weights))


def bootstrap_weighted_metric(
    observations_by_target: dict[str, tuple[np.ndarray, np.ndarray]],
    metric: str,
    iterations: int = 2000,
    seed: int = 20260817,
) -> dict[str, float]:
    """Bootstrap compounds within each target and return a percentile interval."""

    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    for _ in range(iterations):
        resampled: dict[str, dict[str, float]] = {}
        for target, (observed, predicted) in observations_by_target.items():
            indices = rng.integers(0, len(observed), size=len(observed))
            result = assay_metrics(observed[indices], predicted[indices])
            if not np.isfinite(result[metric]):
                break
            resampled[target] = result
        if len(resampled) == len(observations_by_target):
            estimates.append(compound_weighted_average(resampled, metric))
    if not estimates:
        raise ValueError(f"no finite bootstrap estimates for {metric}")
    array = np.asarray(estimates)
    return {
        "iterations_requested": iterations,
        "iterations_finite": int(array.size),
        "seed": seed,
        "lower_95": float(np.quantile(array, 0.025)),
        "median": float(np.quantile(array, 0.5)),
        "upper_95": float(np.quantile(array, 0.975)),
    }


def paired_bootstrap_weighted_metric_difference(
    observations_by_target: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    metric: str,
    iterations: int = 2000,
    seed: int = 20260817,
) -> dict[str, float]:
    """Bootstrap a paired model-A minus model-B weighted metric difference."""

    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    for _ in range(iterations):
        resampled_a: dict[str, dict[str, float]] = {}
        resampled_b: dict[str, dict[str, float]] = {}
        for target, (observed, predicted_a, predicted_b) in observations_by_target.items():
            if not (len(observed) == len(predicted_a) == len(predicted_b)):
                raise ValueError(f"paired observations differ in length for {target}")
            indices = rng.integers(0, len(observed), size=len(observed))
            result_a = assay_metrics(observed[indices], predicted_a[indices])
            result_b = assay_metrics(observed[indices], predicted_b[indices])
            if not (np.isfinite(result_a[metric]) and np.isfinite(result_b[metric])):
                break
            resampled_a[target] = result_a
            resampled_b[target] = result_b
        if len(resampled_a) == len(observations_by_target):
            estimates.append(
                compound_weighted_average(resampled_a, metric)
                - compound_weighted_average(resampled_b, metric)
            )
    if not estimates:
        raise ValueError(f"no finite paired bootstrap estimates for {metric}")
    array = np.asarray(estimates)
    return {
        "iterations_requested": iterations,
        "iterations_finite": int(array.size),
        "seed": seed,
        "lower_95": float(np.quantile(array, 0.025)),
        "median": float(np.quantile(array, 0.5)),
        "upper_95": float(np.quantile(array, 0.975)),
    }
