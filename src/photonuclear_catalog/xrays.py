"""Atomic X-ray energies for distinguishing X-rays from nuclear gammas."""

from __future__ import annotations

# Approximate elemental Kα / Kβ energies in keV.
# Used as a secondary filter when LiveChart marks a line under rad_types=g
# but the energy matches a known atomic X-ray of the decaying element's Z
# (or nearby daughters after β±/EC). Sources: typical characteristic X-ray tables.

_K_ALPHA_KEV: dict[int, float] = {
    26: 6.40,
    27: 6.93,
    28: 7.48,
    29: 8.05,
    30: 8.64,
    46: 21.18,
    47: 22.16,
    48: 23.17,
    49: 24.21,
    50: 25.27,
    51: 26.36,
    52: 27.47,
    53: 28.61,
    54: 29.78,
    55: 30.97,
    56: 32.19,
    57: 33.44,
    58: 34.72,
    74: 59.32,
    75: 61.14,
    76: 63.00,
    77: 64.90,
    78: 66.83,
    79: 68.80,
    80: 70.82,
    81: 72.87,
    82: 74.97,
    83: 77.11,
}

_K_BETA_KEV: dict[int, float] = {
    26: 7.06,
    27: 7.65,
    28: 8.26,
    29: 8.91,
    30: 9.57,
    74: 67.24,
    75: 69.31,
    76: 71.41,
    77: 73.56,
    78: 75.75,
    79: 77.98,
    80: 80.25,
    81: 82.57,
    82: 84.94,
    83: 87.34,
}


def is_atomic_xray(z: int, energy_kev: float, tol_kev: float = 1.25) -> bool:
    """Return True if energy matches a known K X-ray near element Z.

    Checks Z-1, Z, Z+1 because EC/β± daughters often dominate X-ray lines
    attributed to the parent decay dataset.
    Soft L-band for heavy nuclei (Z≥70, Eγ≲14 keV) is treated as atomic.
    """
    for zz in (z - 1, z, z + 1):
        if zz <= 0:
            continue
        for table in (_K_ALPHA_KEV, _K_BETA_KEV):
            ref = table.get(zz)
            if ref is not None and abs(energy_kev - ref) <= tol_kev:
                return True
    # Heavy-element L X-rays cluster well below nuclear IT lines of interest
    # (e.g. 84.66 keV 196Aum1, 98.9 keV 195Ptm1 stay).
    if z >= 70 and energy_kev < 14.0:
        return True
    return False
