"""Photonuclear (γ, xn yp) threshold and product enumeration."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from photonuclear_catalog.coulomb import coulomb_barrier_mev, heuristic_barrier_indicator_mev
from photonuclear_catalog.masses import ME_HYDROGEN_KEV, ME_NEUTRON_KEV, MassTable
from photonuclear_catalog.nubase import (
    NubaseRecord,
    NubaseTable,
    format_half_life,
    pretty_name,
    pretty_name_nuclear,
)
from photonuclear_catalog.nuclide import Nuclide

DEFAULT_SIGMA_NEGLIGIBLE_MB = 1e-3


def _round_mb(value: float | None) -> float | None:
    """Round mb values; keep scientific fidelity for ≪ 1e-3 mb."""
    if value is None:
        return None
    if abs(value) >= 1e-3:
        return round(value, 6)
    return float(f"{value:.6g}")


@dataclass(frozen=True)
class Product:
    target: str
    product: str
    product_pretty: str
    z: int
    a: int
    n: int
    isomer: str  # internal NUBASE label (m/n/...)
    isomer_display: str  # m1/m2/...
    x_neutrons: int
    y_protons: int
    channel: str
    channel_note: str
    eth_mev: float
    coulomb_barrier_mev: float
    heuristic_barrier_indicator_mev: float
    excitation_kev: float
    half_life: str
    is_metastable: bool
    is_stable: bool
    jp: str
    decay: str
    mass_estimated: bool
    half_life_estimated: bool = False
    excitation_estimated: bool = False
    sigma_max_mb: float | None = None
    e_at_sigma_max_mev: float | None = None
    sigma_status: str = "not_requested"
    sigma_scope: str = "missing"  # state_specific | ground_state | missing
    sigma_mt: int | None = None  # free-nucleon ENDF MT expected/used
    sigma_ground_state_reference_mb: float | None = None
    sigma_channel_total_mb: float | None = None  # MF3(MT)
    e_at_channel_total_max_mev: float | None = None
    sigma_sum_mf10_lfs_mb: float | None = None  # Σ MF10 LFS for MT+IZAP
    mf10_coverage_ratio: float | None = None  # ΣMF10 / MF3
    e_at_total_max_mev: float | None = None
    sigma_lfs: int | None = None
    sigma_lfs_excitation_kev: float | None = None
    sigma_qi_mev: float | None = None
    sigma_source_variant: str | None = None
    sigma_source_url: str | None = None
    sigma_source_sha256: str | None = None
    viability: str = "kinematic_open"
    # viability: kinematic_open | kinematic_only | sigma_negligible |
    #            sigma_missing | barrier_uncertain | state_specific_missing |
    #            channel_supported_state_split_missing

    def to_dict(self) -> dict:
        d = asdict(self)
        for key, val in list(d.items()):
            if val is None:
                d[key] = ""
        return d


def threshold_mev(
    masses: MassTable,
    z_t: int,
    a_t: int,
    x: int,
    y: int,
    excitation_kev: float = 0.0,
) -> float | None:
    """Return Eth in MeV for (γ, xn yp), or None if masses missing."""
    z_r = z_t - y
    a_r = a_t - x - y
    if z_r < 0 or a_r < 1 or a_r < z_r:
        return None
    mt = masses.get(z_t, a_t)
    mr = masses.get(z_r, a_r)
    if mt is None or mr is None:
        return None
    eth_kev = (
        mr.mass_excess_kev
        + x * ME_NEUTRON_KEV
        + y * ME_HYDROGEN_KEV
        - mt.mass_excess_kev
        + excitation_kev
    )
    return eth_kev / 1000.0


def classify_viability(
    *,
    eth_mev: float,
    heuristic_barrier_mev: float,
    emax_mev: float,
    sigma_max_mb: float | None,
    sigma_status: str,
    sigma_scope: str,
    sigma_negligible_mb: float = DEFAULT_SIGMA_NEGLIGIBLE_MB,
) -> str:
    """Classify kinematic / TENDL / heuristic-barrier viability.

    TENDL residual σ (when present) takes precedence over the Coulomb heuristic.
    Eth + Vc is only a heuristic indicator, not a sharp second threshold.
    """
    if eth_mev > emax_mev:
        return "kinematic_closed"

    if sigma_status == "channel_total_only":
        return "channel_supported_state_split_missing"

    if sigma_status == "state_specific_missing":
        return "kinematic_only"

    tendl_ok = sigma_status == "ok" and sigma_max_mb is not None
    if tendl_ok:
        if sigma_max_mb < sigma_negligible_mb:
            return "sigma_negligible"
        return "kinematic_open"

    if sigma_status in {
        "residual_not_in_tendl",
        "mt_not_in_tendl",
        "mt_unmapped",
        "mt_unspecified",
        "no_sigma_in_window",
        "tendl_missing",
    }:
        if heuristic_barrier_mev > emax_mev + 1e-9:
            return "barrier_uncertain"
        return "sigma_missing"

    if heuristic_barrier_mev > emax_mev + 1e-9:
        return "barrier_uncertain"
    return "kinematic_open"


def enumerate_products(
    target: Nuclide,
    masses: MassTable,
    nubase: NubaseTable,
    emax_mev: float = 45.0,
    include_stable_products: bool = True,
    tendl: object | None = None,
    sigma_negligible_mb: float = DEFAULT_SIGMA_NEGLIGIBLE_MB,
) -> list[Product]:
    """List residuals reachable with Eth <= emax.

    Includes radioactive gs, metastable isomers, and (by default) stable
    residuals such as 196Pt. Eth is purely kinematic.
    """
    if target.isomer:
        raise ValueError("Target must be a ground-state nuclide (no isomer suffix).")
    if not masses.has(target.z, target.a):
        raise ValueError(f"Target {target.symbol} not found in AME2020 mass table.")

    max_x = min(target.n, 40)
    max_y = min(target.z, 20)

    products: list[Product] = []
    seen: set[tuple[int, int, str, int, int]] = set()

    if tendl is not None:
        try:
            tendl.ensure_target(target.z, target.a)  # type: ignore[attr-defined]
        except Exception:  # noqa: BLE001
            pass

    for y in range(0, max_y + 1):
        for x in range(0, max_x + 1):
            if x == 0 and y == 0:
                continue
            eth_gs = threshold_mev(masses, target.z, target.a, x, y, 0.0)
            if eth_gs is None:
                continue
            if eth_gs > emax_mev:
                if x > 2 and eth_gs > emax_mev + 10:
                    break
                continue

            z_r = target.z - y
            a_r = target.a - x - y
            mt = masses.get(target.z, target.a)
            mr = masses.get(z_r, a_r)
            mass_est = bool((mt and mt.estimated) or (mr and mr.estimated))
            channel = _channel_name(x, y)
            channel_note = _channel_note(x, y)
            barrier = coulomb_barrier_mev(z_r, a_r, y)
            states = nubase.states_for(z_r, a_r)
            if not states:
                states = [
                    NubaseRecord(
                        z=z_r,
                        a=a_r,
                        isomer_index=0,
                        isomer_label="",
                        half_life_raw="?",
                        half_life_unit="",
                        is_stable=False,
                        is_particle_unstable=False,
                        excitation_kev=0.0,
                        abundance=None,
                        jp="",
                        decay="",
                    )
                ]

            for state in states:
                if state.is_particle_unstable:
                    continue
                if state.is_stable and not state.is_metastable:
                    if not include_stable_products:
                        continue
                elif not state.include_in_catalog:
                    continue

                exc = float(state.excitation_kev or 0.0)
                eth = eth_gs if not state.is_metastable else threshold_mev(
                    masses, target.z, target.a, x, y, exc
                )
                if eth is None or eth > emax_mev:
                    continue

                key = (z_r, a_r, state.isomer_label, x, y)
                if key in seen:
                    continue
                seen.add(key)

                heur = heuristic_barrier_indicator_mev(eth, barrier)
                sigma_mb: float | None = None
                e_at: float | None = None
                sigma_status = "not_requested"
                sigma_scope = "missing"
                sigma_mt: int | None = None
                sigma_gs: float | None = None
                sigma_channel: float | None = None
                e_at_channel: float | None = None
                sigma_sum: float | None = None
                coverage: float | None = None
                e_at_total: float | None = None
                sigma_lfs: int | None = None
                lfs_ex: float | None = None
                sigma_qi: float | None = None
                src_variant: str | None = None
                src_url: str | None = None
                src_sha: str | None = None

                if tendl is not None:
                    try:
                        from photonuclear_catalog.tendl import free_nucleon_mt

                        sigma_mt = free_nucleon_mt(x, y)
                        ordinal = 0
                        if state.is_metastable:
                            ordinal = state.isomer_index if state.isomer_index > 0 else 1
                        lookup = tendl.residual_lookup(  # type: ignore[attr-defined]
                            target.z,
                            target.a,
                            z_r,
                            a_r,
                            eth,
                            emax_mev,
                            mt=sigma_mt,
                            x_neutrons=x,
                            y_protons=y,
                            is_metastable=state.is_metastable,
                            isomer_ordinal=ordinal,
                            excitation_kev=exc if state.is_metastable else 0.0,
                            # Channel totals once per (Z,A,MT) in [Eth_gs, Emax];
                            # isomer Eth only clips the state-specific curve.
                            channel_e_lo_mev=eth_gs,
                        )
                        sigma_mb = lookup.sigma_max_mb
                        e_at = lookup.e_at_sigma_max_mev
                        sigma_status = lookup.status
                        sigma_scope = lookup.scope
                        if lookup.mt is not None:
                            sigma_mt = lookup.mt
                        sigma_lfs = lookup.lfs
                        sigma_gs = lookup.sigma_ground_state_mb
                        sigma_channel = lookup.sigma_channel_total_mb
                        e_at_channel = lookup.e_at_channel_total_max_mev
                        sigma_sum = lookup.sigma_sum_mf10_lfs_mb
                        coverage = lookup.mf10_coverage_ratio
                        e_at_total = lookup.e_at_total_max_mev
                        lfs_ex = lookup.lfs_excitation_kev
                        sigma_qi = lookup.sigma_qi_mev
                        src_variant = lookup.source_variant
                        src_url = lookup.source_url
                        src_sha = lookup.source_sha256
                    except Exception as exc:  # noqa: BLE001
                        sigma_mb, sigma_status = None, f"tendl_error:{exc}"
                        sigma_scope = "missing"

                viability = classify_viability(
                    eth_mev=eth,
                    heuristic_barrier_mev=heur,
                    emax_mev=emax_mev,
                    sigma_max_mb=sigma_mb,
                    sigma_status=sigma_status,
                    sigma_scope=sigma_scope,
                    sigma_negligible_mb=sigma_negligible_mb,
                )

                products.append(
                    Product(
                        target=target.symbol,
                        product=pretty_name(
                            z_r,
                            a_r,
                            state.isomer_label,
                            isomer_index=state.isomer_index,
                            display=True,
                        ),
                        product_pretty=pretty_name_nuclear(
                            z_r,
                            a_r,
                            state.isomer_label,
                            isomer_index=state.isomer_index,
                        ),
                        z=z_r,
                        a=a_r,
                        n=a_r - z_r,
                        isomer=state.isomer_label,
                        isomer_display=state.display_isomer,
                        x_neutrons=x,
                        y_protons=y,
                        channel=channel,
                        channel_note=channel_note,
                        eth_mev=round(eth, 6),
                        coulomb_barrier_mev=round(barrier, 4),
                        heuristic_barrier_indicator_mev=round(heur, 6),
                        excitation_kev=exc,
                        half_life=format_half_life(state),
                        is_metastable=state.is_metastable,
                        is_stable=bool(state.is_stable and not state.is_metastable),
                        jp=state.jp,
                        decay=state.decay,
                        mass_estimated=mass_est,
                        half_life_estimated=state.half_life_estimated,
                        excitation_estimated=state.excitation_estimated,
                        sigma_max_mb=_round_mb(sigma_mb),
                        e_at_sigma_max_mev=(
                            None if e_at is None else round(e_at, 4)
                        ),
                        sigma_status=sigma_status,
                        sigma_scope=sigma_scope,
                        sigma_mt=sigma_mt,
                        sigma_ground_state_reference_mb=_round_mb(sigma_gs),
                        sigma_channel_total_mb=_round_mb(sigma_channel),
                        e_at_channel_total_max_mev=(
                            None if e_at_channel is None else round(e_at_channel, 4)
                        ),
                        sigma_sum_mf10_lfs_mb=_round_mb(sigma_sum),
                        mf10_coverage_ratio=coverage,
                        e_at_total_max_mev=(
                            None if e_at_total is None else round(e_at_total, 4)
                        ),
                        sigma_lfs=sigma_lfs,
                        sigma_lfs_excitation_kev=(
                            None if lfs_ex is None else round(lfs_ex, 3)
                        ),
                        sigma_qi_mev=(
                            None if sigma_qi is None else round(sigma_qi, 6)
                        ),
                        sigma_source_variant=src_variant,
                        sigma_source_url=src_url,
                        sigma_source_sha256=src_sha,
                        viability=viability,
                    )
                )

    products.sort(
        key=lambda p: (
            _viability_rank(p.viability),
            0 if not p.is_stable else 1,
            -(p.sigma_max_mb or -1.0),
            p.heuristic_barrier_indicator_mev,
            p.eth_mev,
            p.y_protons,
            p.x_neutrons,
            p.a,
            p.z,
            p.isomer,
        )
    )
    return products


def filter_by_emax(products: list[Product], emax_mev: float) -> list[Product]:
    return [p for p in products if p.eth_mev <= emax_mev]


def _viability_rank(viability: str) -> int:
    order = {
        "kinematic_open": 0,
        "channel_supported_state_split_missing": 1,
        "kinematic_only": 2,
        "sigma_missing": 3,
        "barrier_uncertain": 4,
        "sigma_negligible": 5,
        "kinematic_closed": 6,
    }
    return order.get(viability, 9)


def _channel_name(x: int, y: int) -> str:
    parts: list[str] = []
    if y == 1:
        parts.append("p")
    elif y > 1:
        parts.append(f"{y}p")
    if x == 1:
        parts.append("n")
    elif x > 1:
        parts.append(f"{x}n")
    body = "".join(parts) if parts else "0"
    return f"(γ,{body})"


def _channel_note(x: int, y: int) -> str:
    """Clarify free-nucleon channel bookkeeping vs clustered emission."""
    base = "Free-nucleon balance (n+p count); not a unique reaction mechanism."
    if y >= 2 and x >= 2 and y == x:
        return (
            f"{base} Same residual Z,A can also form via clustered channels "
            f"e.g. (γ,α) when mass equivalent; those are not expanded here."
        )
    if y >= 2 or (y >= 1 and x >= 1):
        return (
            f"{base} Clustered emission (α, d, t, …) is out of scope and may "
            f"populate the same residual."
        )
    return base
