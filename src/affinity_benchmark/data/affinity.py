"""Explicit concentration and logarithmic-affinity unit conversions."""

from __future__ import annotations

import math
from numbers import Real


_UNIT_TO_MOLAR = {
    "M": 1.0,
    "mM": 1e-3,
    "uM": 1e-6,
    "µM": 1e-6,
    "μM": 1e-6,
    "nM": 1e-9,
    "pM": 1e-12,
}


def concentration_to_molar(value: Real, unit: str) -> float:
    """Convert a positive finite concentration to molar units."""

    concentration = _validate_concentration(value)
    try:
        factor = _UNIT_TO_MOLAR[unit]
    except KeyError as error:
        supported = ", ".join(_UNIT_TO_MOLAR)
        raise ValueError(f"unsupported concentration unit {unit!r}; use {supported}") from error
    return concentration * factor


def concentration_to_log10_micromolar(value: Real, unit: str) -> float:
    """Return log10(concentration / 1 micromolar)."""

    molar = concentration_to_molar(value, unit)
    return math.log10(molar / 1e-6)


def concentration_to_px(value: Real, unit: str) -> float:
    """Return pX = -log10(concentration / 1 molar)."""

    return -math.log10(concentration_to_molar(value, unit))


def log10_micromolar_to_px(log_value: Real) -> float:
    """Convert log10(X/micromolar) to pX without changing endpoint identity."""

    numeric = _validate_finite_number(log_value, "log_value")
    return 6.0 - numeric


def px_to_log10_micromolar(px: Real) -> float:
    """Convert pX to log10(X/micromolar) without changing endpoint identity."""

    numeric = _validate_finite_number(px, "px")
    return 6.0 - numeric


def _validate_concentration(value: Real) -> float:
    numeric = _validate_finite_number(value, "concentration")
    if numeric <= 0:
        raise ValueError("concentration must be greater than zero")
    return numeric


def _validate_finite_number(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be finite")
    return numeric
