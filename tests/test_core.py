"""Unit tests for parsers and threshold physics."""

from __future__ import annotations

from pathlib import Path

import pytest

from photonuclear_catalog.masses import MassTable, default_mass_path
from photonuclear_catalog.nubase import NubaseTable, default_nubase_path
from photonuclear_catalog.nuclide import parse_nuclide
from photonuclear_catalog.reactions import enumerate_products, filter_by_emax, threshold_mev


@pytest.fixture(scope="module")
def masses() -> MassTable:
    return MassTable(default_mass_path())


@pytest.fixture(scope="module")
def nubase() -> NubaseTable:
    return NubaseTable(default_nubase_path())


def test_parse_nuclide_forms() -> None:
    assert parse_nuclide("208Pb").symbol == "208Pb"
    assert parse_nuclide("Pb-208").symbol == "208Pb"
    assert parse_nuclide("63Cu").z == 29
    assert parse_nuclide("99mTc").isomer == "m"
    assert parse_nuclide("207Pbm").symbol == "207Pbm"
    assert parse_nuclide("207Pbm").isomer == "m"
    assert parse_nuclide("206Pbn").isomer == "n"
    # Catalog product ids use m1/m2 display suffixes (filter_catalog parses these).
    assert parse_nuclide("207Pbm1").z == 82
    assert parse_nuclide("207Pbm1").a == 207
    assert parse_nuclide("207Pbm1").isomer == "m"
    assert parse_nuclide("196Aum1").isomer == "m"
    assert parse_nuclide("195Ptm1").z == 78 and parse_nuclide("195Ptm1").isomer == "m"
    assert parse_nuclide("152Eum3").isomer == "p"
    assert parse_nuclide("192Irm1").isomer == "m"
    # Two-letter symbols must not be eaten by prefix-isomer "n"/"m".
    assert parse_nuclide("59Ni").z == 28 and parse_nuclide("59Ni").isomer == ""
    assert parse_nuclide("93Nb").z == 41 and parse_nuclide("93Nb").isomer == ""
    assert parse_nuclide("14N").z == 7 and parse_nuclide("14N").isomer == ""
    assert parse_nuclide("23Na").z == 11
    assert parse_nuclide("110nAg").z == 47 and parse_nuclide("110nAg").isomer == "n"
    assert parse_nuclide("99mTc").z == 43 and parse_nuclide("99mTc").isomer == "m"


def test_mass_table_basics(masses: MassTable) -> None:
    assert len(masses) > 3000
    pb = masses.get(82, 208)
    assert pb is not None
    assert pb.element.strip() == "Pb"
    n = masses.get(0, 1)
    assert n is not None


def test_nubase_stable_and_isomer(nubase: NubaseTable) -> None:
    assert nubase.is_stable_ground(82, 208)
    tc = nubase.get(43, 99, "m")
    assert tc is not None
    assert tc.is_metastable
    assert not tc.is_stable


def test_threshold_pb208_gn(masses: MassTable) -> None:
    # 208Pb(γ,n)207Pb threshold is around 7.4 MeV.
    eth = threshold_mev(masses, 82, 208, x=1, y=0)
    assert eth is not None
    assert 7.0 < eth < 8.0


def test_enumerate_filters_energy_and_radioactive(
    masses: MassTable, nubase: NubaseTable
) -> None:
    target = parse_nuclide("208Pb")
    products_45 = enumerate_products(target, masses, nubase, emax_mev=45.0)
    assert products_45
    assert all(p.eth_mev <= 45.0 for p in products_45)
    assert all(
        p.is_stable or p.half_life != "stable" or p.is_metastable for p in products_45
    )
    assert all(hasattr(p, "heuristic_barrier_indicator_mev") for p in products_45)
    assert any(p.isomer_display.startswith("m") for p in products_45 if p.is_metastable)

    products_10 = filter_by_emax(products_45, 10.0)
    assert products_10
    assert len(products_10) < len(products_45)
    assert all(p.eth_mev <= 10.0 for p in products_10)
    assert any(p.channel == "(γ,n)" for p in products_10)


def test_channels_include_mixed(masses: MassTable, nubase: NubaseTable) -> None:
    target = parse_nuclide("63Cu")
    products = enumerate_products(target, masses, nubase, emax_mev=45.0)
    channels = {p.channel for p in products}
    assert "(γ,n)" in channels or "(γ,2n)" in channels
    assert any(p.y_protons > 0 for p in products)


def test_coulomb_and_tendl_annotation(masses: MassTable, nubase: NubaseTable) -> None:
    from photonuclear_catalog.coulomb import coulomb_barrier_mev
    from photonuclear_catalog.tendl import TendlPhotonuclear

    target = parse_nuclide("197Au")
    tendl = TendlPhotonuclear(cache_dir="data/tendl_cache")
    products = enumerate_products(
        target, masses, nubase, emax_mev=44.0, tendl=tendl
    )
    assert products
    # Stable residuals included (user request).
    assert any(p.product == "196Pt" and p.is_stable for p in products)
    charged = [p for p in products if p.y_protons > 0]
    assert charged
    assert all(p.coulomb_barrier_mev > 0 for p in charged)
    assert all(
        p.heuristic_barrier_indicator_mev >= p.eth_mev for p in products
    )
    au196 = [p for p in products if p.z == 79 and p.a == 196 and not p.isomer]
    assert au196
    assert au196[0].sigma_status == "ok"
    assert au196[0].sigma_scope == "ground_state"
    assert au196[0].sigma_mt == 4
    assert au196[0].sigma_max_mb is not None and au196[0].sigma_max_mb > 1.0
    assert au196[0].e_at_sigma_max_mev is not None
    assert au196[0].sigma_sum_mf10_lfs_mb is not None
    assert au196[0].sigma_channel_total_mb is not None
    # Total Σ_LFS must be >= ground-state LFS=0; coverage ≤ ~1 vs MF3.
    assert au196[0].sigma_sum_mf10_lfs_mb + 1e-9 >= au196[0].sigma_max_mb
    assert au196[0].mf10_coverage_ratio is not None
    assert au196[0].mf10_coverage_ratio <= 1.0 + 1e-3
    # Isomer m1 should be state_specific (LFS>0), not a copy of LFS=0.
    aum = [p for p in products if p.product == "196Aum1"]
    assert aum
    assert aum[0].sigma_scope == "state_specific"
    assert aum[0].sigma_max_mb is not None
    assert aum[0].sigma_max_mb != au196[0].sigma_max_mb
    # Isomer without LFS match must not inherit LFS=0 as sigma_max.
    aum2 = [p for p in products if p.product == "196Aum2"]
    assert aum2
    assert aum2[0].sigma_status == "state_specific_missing"
    assert aum2[0].sigma_max_mb is None
    assert aum2[0].viability == "kinematic_only"
    assert aum2[0].sigma_ground_state_reference_mb == au196[0].sigma_max_mb
    # Charged free-nucleon MT must not pick clustered (γ,d)/(γ,α) curves.
    pt195 = [p for p in products if p.product == "195Pt" and p.channel == "(γ,pn)"]
    assert pt195
    assert pt195[0].sigma_mt == 28
    assert pt195[0].sigma_status == "ok"
    assert abs(pt195[0].sigma_max_mb - 0.283225) < 1e-4
    assert abs((pt195[0].e_at_sigma_max_mev or 0) - 35.0) < 1e-6
    assert pt195[0].sigma_source_variant == "s60"
    ptm = [p for p in products if p.product == "195Ptm1" and p.channel == "(γ,pn)"]
    assert ptm
    assert ptm[0].sigma_mt == 28
    # Emax=44 interpolation raises max above the 40 MeV tabulated point.
    assert ptm[0].sigma_max_mb is not None
    assert abs(ptm[0].sigma_max_mb - 0.233263) < 1e-4
    # s60 recovers MT152 / MT42 that s30/mt200 falsely marked missing.
    au192 = [p for p in products if p.product == "192Au"]
    assert au192 and au192[0].sigma_mt == 152
    assert au192[0].sigma_status == "ok"
    assert abs((au192[0].sigma_max_mb or 0) - 0.853975) < 1e-4
    aum2 = [p for p in products if p.product == "192Aum2"]
    assert aum2 and aum2[0].sigma_lfs == 15
    assert abs((aum2[0].sigma_max_mb or 0) - 0.0240126) < 1e-5
    pt193 = [p for p in products if p.product == "193Pt"]
    assert pt193 and pt193[0].sigma_mt == 42
    assert abs((pt193[0].sigma_max_mb or 0) - 0.0436598) < 1e-5
    # Tiny but present free-nucleon σ → sigma_negligible (not mt_not_in_tendl).
    ir194 = [p for p in products if p.z == 77 and p.a == 194 and not p.isomer and p.y_protons == 2]
    assert ir194
    assert ir194[0].sigma_mt == 44
    assert ir194[0].sigma_status == "ok"
    assert ir194[0].viability == "sigma_negligible"
    assert ir194[0].sigma_max_mb is not None
    assert ir194[0].sigma_max_mb < 1e-3
    ir192 = [p for p in products if p.z == 77 and p.a == 192 and not p.isomer and p.y_protons == 2]
    assert ir192
    assert ir192[0].sigma_mt == 179
    assert ir192[0].sigma_status == "ok"
    assert ir192[0].viability == "sigma_negligible"
    # MF3 exists, MF10 state split missing → channel_total_only.
    pt196 = [p for p in products if p.product == "196Pt" and p.y_protons == 1 and p.x_neutrons == 0]
    assert pt196
    assert pt196[0].sigma_mt == 103
    assert pt196[0].sigma_status == "channel_total_only"
    assert pt196[0].sigma_max_mb is None
    assert abs((pt196[0].sigma_channel_total_mb or 0) - 0.351179) < 1e-4
    assert pt196[0].viability == "channel_supported_state_split_missing"
    # 195Aum2 must not clip channel totals with isomer Eth.
    au195 = [p for p in products if p.product == "195Au"]
    au195m2 = [p for p in products if p.product == "195Aum2"]
    assert au195 and au195m2
    assert au195[0].sigma_mt == 16
    assert abs((au195[0].sigma_ground_state_reference_mb or 0) - 109.4047) < 1e-3
    assert abs((au195[0].sigma_sum_mf10_lfs_mb or 0) - 129.6231) < 1e-3
    assert au195m2[0].sigma_ground_state_reference_mb == au195[0].sigma_ground_state_reference_mb
    assert au195m2[0].sigma_sum_mf10_lfs_mb == au195[0].sigma_sum_mf10_lfs_mb
    assert coulomb_barrier_mev(77, 192, 2) > 10.0


def test_free_nucleon_mt_map() -> None:
    from photonuclear_catalog.tendl import free_nucleon_mt

    assert free_nucleon_mt(1, 0) == 4
    assert free_nucleon_mt(1, 1) == 28
    assert free_nucleon_mt(1, 2) == 44
    assert free_nucleon_mt(3, 2) == 179
    assert free_nucleon_mt(0, 1) == 103
    # Unmapped / clustered-equivalent balances stay unmapped.
    assert free_nucleon_mt(0, 0) is None
    assert free_nucleon_mt(9, 9) is None


def test_filter_catalog_accepts_m1_product_ids(
    masses: MassTable, nubase: NubaseTable
) -> None:
    """Filtering must parse catalog ids like 207Pbm1 (not crash)."""
    from photonuclear_catalog.catalog import CatalogResult, filter_catalog
    from photonuclear_catalog.gammas import GammaLine

    products = enumerate_products(
        parse_nuclide("208Pb"), masses, nubase, emax_mev=20.0, include_stable_products=True
    )
    result = CatalogResult(
        target="208Pb",
        emax_mev=20.0,
        products=products,
        gammas=[
            GammaLine(
                product="207Pbm1",
                energy_kev=569.698,
                intensity=97.9,
                intensity_unc=None,
                decay_mode="IT",
                parent_energy_kev=1633.356,
                half_life="806 ms",
                rank=1,
            )
        ],
        use_tendl=False,
        meta={"include_stable_products": True, "exclude_xrays": True},
    )
    filtered = filter_catalog(result, 10.0, masses=masses, nubase=nubase)
    assert filtered.emax_mev == 10.0
    assert all(p.eth_mev <= 10.0 for p in filtered.products)
    assert filtered.meta.get("filtered_from_emax_mev") == 20.0
    # 207Pbm Eth is ~9 MeV-class; gamma row should survive when isomer remains.
    remaining = {(p.z, p.a, p.isomer) for p in filtered.products}
    if (82, 207, "m") in remaining:
        assert any(g.product == "207Pbm1" for g in filtered.gammas)


def test_filter_catalog_refreshes_tendl_to_new_emax(
    masses: MassTable, nubase: NubaseTable
) -> None:
    """Filter 44→20 MeV must recompute σ (not keep the 44 MeV window)."""
    from photonuclear_catalog.catalog import build_catalog, filter_catalog

    cache = Path("data/tendl_cache")
    wide = build_catalog(
        "197Au",
        emax_mev=44.0,
        fetch_gammas=False,
        masses=masses,
        nubase=nubase,
        use_tendl=True,
        tendl_cache_dir=cache,
    )
    direct = build_catalog(
        "197Au",
        emax_mev=20.0,
        fetch_gammas=False,
        masses=masses,
        nubase=nubase,
        use_tendl=True,
        tendl_cache_dir=cache,
    )
    filtered = filter_catalog(
        wide, 20.0, masses=masses, nubase=nubase, tendl_cache_dir=cache
    )

    def _ptm(rows):
        return [p for p in rows if p.product == "195Ptm1" and p.channel == "(γ,pn)"]

    wide_pt = _ptm(wide.products)
    direct_pt = _ptm(direct.products)
    filt_pt = _ptm(filtered.products)
    assert wide_pt and direct_pt and filt_pt

    # Stale-bug signature: filtered would keep the wide σ (~0.23 mb) instead of ~3e-5.
    assert wide_pt[0].sigma_max_mb is not None
    assert direct_pt[0].sigma_max_mb is not None
    assert abs(wide_pt[0].sigma_max_mb - direct_pt[0].sigma_max_mb) > 1e-4

    assert filt_pt[0].sigma_max_mb == direct_pt[0].sigma_max_mb
    assert filt_pt[0].e_at_sigma_max_mev == direct_pt[0].e_at_sigma_max_mev
    assert filt_pt[0].viability == direct_pt[0].viability
    assert filtered.emax_mev == 20.0

