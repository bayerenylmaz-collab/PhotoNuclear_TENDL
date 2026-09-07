"""Tests for human-readable reports and nuclear-gamma filters."""

from __future__ import annotations

import json
from pathlib import Path

from photonuclear_catalog.catalog import CatalogResult, export_result
from photonuclear_catalog.gammas import GammaClient, GammaLine
from photonuclear_catalog.nuclide import parse_nuclide
from photonuclear_catalog.reactions import Product
from photonuclear_catalog.report import render_markdown


def _sample_product(**kwargs) -> Product:
    base = dict(
        target="208Pb",
        product="207Pbm1",
        product_pretty="^{207m1}Pb",
        z=82,
        a=207,
        n=125,
        isomer="m",
        isomer_display="m1",
        x_neutrons=1,
        y_protons=0,
        channel="(γ,n)",
        channel_note="Free-nucleon balance",
        eth_mev=9.001,
        coulomb_barrier_mev=0.0,
        heuristic_barrier_indicator_mev=9.001,
        excitation_kev=1633.356,
        half_life="806 ms",
        is_metastable=True,
        is_stable=False,
        jp="13/2+",
        decay="IT=100",
        mass_estimated=False,
        sigma_max_mb=100.0,
        e_at_sigma_max_mev=14.0,
        sigma_status="ok",
        sigma_scope="state_specific",
        viability="kinematic_open",
    )
    base.update(kwargs)
    return Product(**base)  # type: ignore[arg-type]


def _sample_result() -> CatalogResult:
    products = [
        _sample_product(),
        _sample_product(
            product="206Hgm1",
            product_pretty="^{206m1}Hg",
            z=80,
            a=206,
            n=126,
            x_neutrons=0,
            y_protons=2,
            channel="(γ,2p)",
            eth_mev=17.483,
            coulomb_barrier_mev=12.5,
            heuristic_barrier_indicator_mev=29.983,
            excitation_kev=2102.4,
            half_life="2.088 us",
            jp="",
            decay="",
            sigma_max_mb=None,
            e_at_sigma_max_mev=None,
            sigma_status="residual_not_in_tendl",
            sigma_scope="missing",
            viability="barrier_uncertain",
        ),
    ]
    gammas = [
        GammaLine("207Pbm1", 569.698, 97.9, 1.4, "IT", 1633.356, "806 ms", 1),
        GammaLine("207Pbm1", 1063.656, 88.8, 1.3, "IT", 1633.356, "806 ms", 2),
    ]
    return CatalogResult(
        target="208Pb",
        emax_mev=20.0,
        products=products,
        gammas=gammas,
        min_gamma_kev=0.0,
        gamma_limit=10,
        use_tendl=True,
    )


def test_markdown_report_has_per_isotope_blocks() -> None:
    md = render_markdown(_sample_result(), min_gamma_kev=0.0, gamma_limit=10)
    assert "207Pbm1" in md or "207m1" in md
    assert "569.7" in md or "569.698" in md
    assert "Bariyer göstergesi" in md
    assert "barrier_uncertain" in md or "Uygunluk" in md


def test_export_writes_html_and_md(tmp_path: Path) -> None:
    paths = export_result(_sample_result(), tmp_path)
    assert paths["report_md"].exists()
    assert paths["report_html"].exists()
    html = paths["report_html"].read_text(encoding="utf-8")
    assert "207Pbm1" in html or "207m1" in html
    assert "<table>" in html
    # Values stay left-aligned under headers (right-align looked "shifted")
    assert "text-align: left" in html
    assert "text-align: right" not in html
    assert "font-variant-numeric: tabular-nums;" in html


def test_report_prose_has_no_changelog_language(tmp_path: Path) -> None:
    paths = export_result(_sample_result(), tmp_path)
    html = paths["report_html"].read_text(encoding="utf-8")
    md = paths["report_md"].read_text(encoding="utf-8")
    for text in (html, md):
        low = text.lower()
        assert "eski etiket" not in low
        assert "decay_rads" not in low
        assert "rad_types=x" not in low
        assert "heuristic_barrier_indicator" not in text
        assert "parent_match_conflict:" not in text
        assert "Bariyer göstergesi" in text or "bariyer göstergesi" in low



def test_xray_removed_via_livechart_x_set(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    g_rows = [
        {
            "energy": "65.122",
            "intensity": "21.9",
            "unc_i": "",
            "decay": "EC+B+",
            "decay_%": "93",
            "p_energy": "0",
            "p_z": "79",
            "p_n": "117",
            "p_symbol": "Au",
            "half_life": "6.1669",
            "unit_hl": "d",
            "half_life_sec": "532800",
        },
        {
            "energy": "355.73",
            "intensity": "87",
            "unc_i": "",
            "decay": "EC+B+",
            "decay_%": "93",
            "p_energy": "0",
            "p_z": "79",
            "p_n": "117",
            "p_symbol": "Au",
            "half_life": "6.1669",
            "unit_hl": "d",
            "half_life_sec": "532800",
        },
    ]
    x_rows = [
        {
            "energy": "65.122",
            "intensity": "21.9",
            "decay": "EC+B+",
            "p_energy": "0",
            "p_z": "79",
            "p_n": "117",
            "p_symbol": "Au",
        }
    ]
    (cache / "196au.json").write_text(json.dumps(g_rows), encoding="utf-8")
    (cache / "196au_x.json").write_text(json.dumps(x_rows), encoding="utf-8")
    client = GammaClient(cache_dir=cache)
    lines = client.top_gammas(parse_nuclide("196Au"), limit=10, exclude_xrays=True)
    energies = [g.energy_kev for g in lines if g.rank > 0]
    assert 355.73 in energies
    assert all(abs(e - 65.122) > 0.1 for e in energies)
    # Absolute Iγ kept (not 87*0.93)
    by_e = {g.energy_kev: g for g in lines if g.rank > 0}
    assert by_e[355.73].intensity == 87.0
    assert by_e[355.73].normalization_status == "absolute"


def test_xray_only_status_when_parent_matched(tmp_path: Path) -> None:
    """Parent matched but every g-line is an atomic X-ray → status=xray_only."""
    cache = tmp_path / "cache"
    cache.mkdir()
    g_rows = [
        {
            "energy": "6.93",
            "intensity": "30",
            "unc_i": "",
            "decay": "EC",
            "decay_%": "100",
            "p_energy": "0",
            "p_z": "28",
            "p_n": "31",
            "p_symbol": "Ni",
            "half_life": "76000",
            "unit_hl": "y",
            "half_life_sec": "2.4e12",
        }
    ]
    x_rows = [
        {
            "energy": "6.93",
            "intensity": "30",
            "decay": "EC",
            "p_energy": "0",
            "p_z": "28",
            "p_n": "31",
            "p_symbol": "Ni",
        }
    ]
    (cache / "59ni.json").write_text(json.dumps(g_rows), encoding="utf-8")
    (cache / "59ni_x.json").write_text(json.dumps(x_rows), encoding="utf-8")
    client = GammaClient(cache_dir=cache)
    result = client.query_gammas(parse_nuclide("59Ni"), limit=10, exclude_xrays=True)
    assert result.status == "xray_only"
    assert "Parent matched" in result.detail
    assert "no_matched_ensdf_parent" not in result.status

    lines = client.top_gammas(parse_nuclide("59Ni"), limit=10, exclude_xrays=True)
    assert len(lines) == 1
    assert lines[0].status == "xray_only"
    assert lines[0].rank is None

    from photonuclear_catalog.report import _gamma_status_message

    msg = _gamma_status_message(lines)
    assert "Üst seviye eşleşti" in msg or "üst seviye eşleşti" in msg.lower()
    assert "eski etiket" not in msg.lower()
    assert "decay_rads" not in msg.lower()
    assert "yalnızca atomik X-ışınları" in msg or "yalnız atomik" in msg


def test_branch_relative_only_when_detected(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    rows = [
        {
            "energy": "316.5",
            "intensity": "92.24",
            "unc_i": "1",
            "decay": "B-",
            "decay_%": "0.0175",
            "p_energy": "56.72",
            "p_z": "77",
            "p_n": "115",
            "p_symbol": "Ir",
            "half_life": "1.45",
            "unit_hl": "m",
            "half_life_sec": "87",
        },
        {
            "energy": "56.71",
            "intensity": "0.0351",
            "unc_i": "",
            "decay": "IT",
            "decay_%": "99.9825",
            "p_energy": "56.72",
            "p_z": "77",
            "p_n": "115",
            "p_symbol": "Ir",
            "half_life": "1.45",
            "unit_hl": "m",
            "half_life_sec": "87",
        },
    ]
    (cache / "192ir.json").write_text(json.dumps(rows), encoding="utf-8")
    (cache / "192ir_x.json").write_text(json.dumps([]), encoding="utf-8")
    client = GammaClient(cache_dir=cache)
    lines = client.top_gammas(
        parse_nuclide("192Irm"),
        limit=10,
        excitation_kev=56.72,
        exclude_xrays=True,
    )
    by_e = {round(g.energy_kev, 1): g for g in lines if g.rank > 0}
    assert by_e[316.5].normalization_status == "branch_relative"
    assert round(by_e[316.5].intensity, 5) == round(92.24 * 0.0175 / 100.0, 5)
    assert by_e[56.7].normalization_status == "absolute"
    assert by_e[56.7].intensity == 0.0351
