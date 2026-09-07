"""Coulomb barrier estimates for charged-particle emission channels."""

from __future__ import annotations


def coulomb_barrier_mev(
    z_residual: int,
    a_residual: int,
    y_protons: int,
    r0_fm: float = 1.45,
) -> float:
    """Approximate Coulomb barrier scale for emitting y protons.

    Uses a simple touching-sphere estimate:
        Vc ≈ 1.44 * Z_r * y / (r0 (A_r^{1/3} + y^{1/3})) MeV
    Neutrons (y=0) return 0.

    This is a barrier *scale*, not a sharp reaction threshold. Tunneling
    allows yield below Eth+Vc; prefer TENDL/TALYS σ when available.
    """
    if y_protons <= 0:
        return 0.0
    if z_residual < 0 or a_residual <= 0:
        return 0.0
    radius = r0_fm * (a_residual ** (1.0 / 3.0) + (y_protons ** (1.0 / 3.0)))
    if radius <= 0:
        return 0.0
    return 1.44 * z_residual * y_protons / radius


def heuristic_barrier_indicator_mev(eth_mev: float, barrier_mev: float) -> float:
    """Eth + Vc heuristic indicator (not a sharp second threshold)."""
    return eth_mev + max(0.0, barrier_mev)


# Backward-compatible alias
def effective_threshold_mev(eth_mev: float, barrier_mev: float) -> float:
    return heuristic_barrier_indicator_mev(eth_mev, barrier_mev)
