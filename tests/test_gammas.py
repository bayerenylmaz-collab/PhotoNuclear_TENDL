"""Tests for decay-gamma parent-energy selection."""

from __future__ import annotations

import json
from pathlib import Path

from photonuclear_catalog.gammas import GammaClient, _parent_energy_matches
from photonuclear_catalog.nuclide import parse_nuclide


def test_parent_energy_match_gs_and_isomer() -> None:
    assert _parent_energy_matches(0.0, want_isomer=False, excitation_kev=None)
    assert _parent_energy_matches(None, want_isomer=False, excitation_kev=None)
    assert not _parent_energy_matches(58.59, want_isomer=False, excitation_kev=None)

    assert _parent_energy_matches(
        1633.356, want_isomer=True, excitation_kev=1633.356
    )
    assert not _parent_energy_matches(0.0, want_isomer=True, excitation_kev=1633.356)
    assert _parent_energy_matches(100.0, want_isomer=True, excitation_kev=None)


def test_top_gammas_uses_gs_payload_for_isomer(tmp_path: Path) -> None:
    # Simulate LiveChart: isomer id empty; gs id carries both parents.
    cache = tmp_path / "cache"
    cache.mkdir()
    rows = [
        {
            "energy": "1332.492",
            "intensity": "99.9826",
            "unc_i": "",
            "decay": "B-",
            "decay_%": "100",
            "p_energy": "0",
            "half_life": "5.27",
            "unit_hl": "y",
            "half_life_sec": "1.66e8",
        },
        {
            "energy": "1173.228",
            "intensity": "99.85",
            "unc_i": "",
            "decay": "B-",
            "decay_%": "100",
            "p_energy": "0",
            "half_life": "5.27",
            "unit_hl": "y",
            "half_life_sec": "1.66e8",
        },
        {
            "energy": "58.59",
            "intensity": "2.0",
            "unc_i": "",
            "decay": "IT",
            "decay_%": "100",
            "p_energy": "58.59",
            "half_life": "10.5",
            "unit_hl": "m",
            "half_life_sec": "630",
        },
        {
            "energy": "1332.5",
            "intensity": "0.2",
            "unc_i": "",
            "decay": "IT",
            "decay_%": "100",
            "p_energy": "58.59",
            "half_life": "10.5",
            "unit_hl": "m",
            "half_life_sec": "630",
        },
        {
            "energy": "997.1",
            "intensity": "",
            "unc_i": "",
            "decay": "IT",
            "decay_%": "100",
            "p_energy": "1348.18",
            "half_life": "1.33",
            "unit_hl": "s",
            "half_life_sec": "1.33",
        },
    ]
    (cache / "60co.json").write_text(json.dumps(rows), encoding="utf-8")
    (cache / "207tl.json").write_text(
        json.dumps(
            [
                {
                    "energy": "897.77",
                    "intensity": "0.263",
                    "unc_i": "",
                    "decay": "B-",
                    "decay_%": "100",
                    "p_energy": "0",
                    "half_life": "4.77",
                    "unit_hl": "m",
                    "half_life_sec": "286",
                },
                {
                    "energy": "997.1",
                    "intensity": "",
                    "unc_i": "",
                    "decay": "IT",
                    "decay_%": "100",
                    "p_energy": "1348.18",
                    "half_life": "1.33",
                    "unit_hl": "s",
                    "half_life_sec": "1.33",
                },
                {
                    "energy": "351.07",
                    "intensity": "",
                    "unc_i": "",
                    "decay": "IT",
                    "decay_%": "100",
                    "p_energy": "1348.18",
                    "half_life": "1.33",
                    "unit_hl": "s",
                    "half_life_sec": "1.33",
                },
            ]
        ),
        encoding="utf-8",
    )
    client = GammaClient(cache_dir=cache)

    gs = client.top_gammas(
        parse_nuclide("60Co"), limit=10, excitation_kev=0.0, min_energy_kev=0.0
    )
    assert [g.energy_kev for g in gs[:2]] == [1332.492, 1173.228]
    assert all((g.parent_energy_kev or 0) < 1 for g in gs)

    iso = client.top_gammas(
        parse_nuclide("60Com"), limit=10, excitation_kev=58.59, min_energy_kev=0.0
    )
    assert iso
    assert all(abs((g.parent_energy_kev or 0) - 58.59) < 1 for g in iso)
    assert iso[0].energy_kev == 58.59

    tlm = client.top_gammas(
        parse_nuclide("207Tlm"),
        limit=10,
        excitation_kev=1348.18,
        min_energy_kev=0.0,
    )
    assert len(tlm) == 2
    assert {round(g.energy_kev, 2) for g in tlm} == {997.1, 351.07}
    assert all(g.intensity is None for g in tlm)
    assert all(g.rank is None for g in tlm)
    assert all(g.uncertainty_status == "intensity_missing" for g in tlm)


def test_isomer_adopted_levels_fallback_when_decay_rads_empty(tmp_path: Path) -> None:
    """Stable residual isomer: empty decay_rads → Adopted Levels relative γ."""
    from photonuclear_catalog.nuclide import Nuclide

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "54fe.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "54fe_x.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "54fe_levels.json").write_text(
        json.dumps(
            [
                {
                    "energy": "6527.1",
                    "jp": "10+",
                    "half_life": "364",
                    "unit_hl": "ns",
                    "half_life_sec": "3.64e-7",
                }
            ]
        ),
        encoding="utf-8",
    )
    (cache / "54fe_gammas.json").write_text(
        json.dumps(
            [
                {
                    "start_level_energy": "6527.1",
                    "end_level_energy": "6380.9",
                    "energy": "146.2",
                    "relative_intensity": "100",
                    "unc_ri": "0.3",
                    "multipolarity": "E2",
                },
                {
                    "start_level_energy": "6527.1",
                    "end_level_energy": "2949.2",
                    "energy": "3577.6",
                    "relative_intensity": "2",
                    "unc_ri": "0.2",
                    "multipolarity": "E4",
                },
                {
                    "start_level_energy": "0",
                    "energy": "999.0",
                    "relative_intensity": "50",
                    "unc_ri": "",
                    "multipolarity": "",
                },
            ]
        ),
        encoding="utf-8",
    )
    client = GammaClient(cache_dir=cache)
    nuc = Nuclide(z=26, a=54, isomer="m")
    result = client.query_gammas(
        nuc, limit=10, excitation_kev=6527.1, expected_half_life="364 ns"
    )
    assert result.status == "relative_adopted_level"
    assert len(result.lines) == 2
    assert [g.energy_kev for g in result.lines] == [146.2, 3577.6]
    assert result.lines[0].intensity == 100.0
    assert result.lines[0].normalization_status == "relative_adopted_level"
    assert result.lines[0].status == "relative_adopted_level"
    assert "no_matched_ensdf_parent" not in result.status

    lines = client.top_gammas(
        nuc, limit=10, excitation_kev=6527.1, expected_half_life="364 ns"
    )
    assert len(lines) == 2
    assert lines[0].rank == 1


def test_higher_isomer_adopted_fallback_when_only_other_parents_in_decay_rads(
    tmp_path: Path,
) -> None:
    """⁵⁸ᵐ²Co-like: decay_rads has gs/m1 only; m2 comes from Adopted Levels."""
    from photonuclear_catalog.nuclide import Nuclide

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "58co.json").write_text(
        json.dumps(
            [
                {
                    "energy": "810.76",
                    "intensity": "99.45",
                    "decay": "EC+B+",
                    "decay_%": "100",
                    "p_energy": "0",
                    "half_life": "70.86",
                    "unit_hl": "d",
                    "half_life_sec": "6.12e6",
                },
                {
                    "energy": "24.889",
                    "intensity": "0.0397",
                    "decay": "IT",
                    "decay_%": "100",
                    "p_energy": "24.95",
                    "half_life": "9.10",
                    "unit_hl": "h",
                    "half_life_sec": "32760",
                },
            ]
        ),
        encoding="utf-8",
    )
    (cache / "58co_x.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "58co_levels.json").write_text(
        json.dumps(
            [
                {
                    "energy": "53.15",
                    "jp": "5+",
                    "half_life": "9.10",
                    "unit_hl": "h",
                    "half_life_sec": "32760",
                }
            ]
        ),
        encoding="utf-8",
    )
    (cache / "58co_gammas.json").write_text(
        json.dumps(
            [
                {
                    "start_level_energy": "53.15",
                    "energy": "52.96",
                    "relative_intensity": "100",
                    "unc_ri": "",
                    "multipolarity": "[E2]",
                },
                {
                    "start_level_energy": "53.15",
                    "energy": "28.3",
                    "relative_intensity": "43",
                    "unc_ri": "5",
                    "multipolarity": "E2+M1",
                },
            ]
        ),
        encoding="utf-8",
    )
    client = GammaClient(cache_dir=cache)
    nuc = Nuclide(z=27, a=58, isomer="m2")
    result = client.query_gammas(
        nuc, limit=10, excitation_kev=53.15, expected_half_life="9.10 h"
    )
    assert result.status == "relative_adopted_level"
    assert [round(g.energy_kev, 2) for g in result.lines] == [52.96, 28.3]
    assert result.lines[0].normalization_status == "relative_adopted_level"


def test_adopted_level_without_gammas_is_not_no_matched_parent(tmp_path: Path) -> None:
    from photonuclear_catalog.nuclide import Nuclide

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "54fe.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "54fe_x.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "54fe_levels.json").write_text(
        json.dumps(
            [
                {
                    "energy": "6527.1",
                    "jp": "10+",
                    "half_life": "364",
                    "unit_hl": "ns",
                    "half_life_sec": "3.64e-7",
                }
            ]
        ),
        encoding="utf-8",
    )
    (cache / "54fe_gammas.json").write_text(json.dumps([]), encoding="utf-8")
    client = GammaClient(cache_dir=cache)
    result = client.query_gammas(
        Nuclide(z=26, a=54, isomer="m"),
        limit=10,
        excitation_kev=6527.1,
        expected_half_life="364 ns",
    )
    assert result.status == "no_absolute_decay_gamma_dataset"
    assert "no_matched_ensdf_parent" not in result.status


def test_select_closest_parent_among_decay_rads_tolerance_window(tmp_path: Path) -> None:
    """±2 keV may contain neighbors; keep only the closest parent (¹⁵²ᵐ³Eu)."""
    from photonuclear_catalog.nuclide import Nuclide

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "152eu.json").write_text(
        json.dumps(
            [
                {
                    "energy": "77.2583",
                    "intensity": "100",
                    "decay": "IT",
                    "decay_%": "100",
                    "p_energy": "77.2593",
                    "half_life": "38",
                    "unit_hl": "ns",
                    "half_life_sec": "3.8e-8",
                },
                {
                    "energy": "32.6341",
                    "intensity": "100",
                    "decay": "IT",
                    "decay_%": "100",
                    "p_energy": "78.2331",
                    "half_life": "165",
                    "unit_hl": "ns",
                    "half_life_sec": "1.65e-7",
                },
                {
                    "energy": "12.965",
                    "intensity": "50",
                    "decay": "IT",
                    "decay_%": "100",
                    "p_energy": "78.2331",
                    "half_life": "165",
                    "unit_hl": "ns",
                    "half_life_sec": "1.65e-7",
                },
            ]
        ),
        encoding="utf-8",
    )
    (cache / "152eu_x.json").write_text(json.dumps([]), encoding="utf-8")
    client = GammaClient(cache_dir=cache)
    result = client.query_gammas(
        Nuclide(z=63, a=152, isomer="m3"),
        limit=10,
        excitation_kev=78.2331,
        expected_half_life="165 ns",
        min_energy_kev=0.0,
    )
    energies = [round(g.energy_kev, 4) for g in result.lines]
    assert 77.2583 not in energies
    assert 32.6341 in energies
    assert 12.965 in energies
    assert all(abs((g.parent_energy_kev or 0) - 78.2331) < 0.05 for g in result.lines)
    assert all(g.parent_match_conflict for g in result.lines)


def test_adopted_levels_picks_closest_not_all_within_tolerance(tmp_path: Path) -> None:
    """¹⁵¹ᵐ¹Eu: do not merge neighbor level 196.54 into Exc=196.245."""
    from photonuclear_catalog.nuclide import Nuclide

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "151eu.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "151eu_x.json").write_text(json.dumps([]), encoding="utf-8")
    (cache / "151eu_levels.json").write_text(
        json.dumps(
            [
                {
                    "energy": "196.245",
                    "jp": "11/2-",
                    "half_life": "58.9",
                    "unit_hl": "us",
                    "half_life_sec": "5.89e-5",
                },
                {
                    "energy": "196.54",
                    "jp": "(7/2+)",
                    "half_life": "",
                    "unit_hl": "",
                    "half_life_sec": "",
                },
            ]
        ),
        encoding="utf-8",
    )
    (cache / "151eu_gammas.json").write_text(
        json.dumps(
            [
                {
                    "start_level_energy": "196.245",
                    "energy": "174.70",
                    "relative_intensity": "100",
                    "unc_ri": "",
                    "multipolarity": "M2",
                },
                {
                    "start_level_energy": "196.245",
                    "energy": "196.0",
                    "relative_intensity": "1.49",
                    "unc_ri": "",
                    "multipolarity": "E3",
                },
                {
                    "start_level_energy": "196.54",
                    "energy": "196.54",
                    "relative_intensity": "100",
                    "unc_ri": "",
                    "multipolarity": "",
                },
                {
                    "start_level_energy": "196.54",
                    "energy": "175.0",
                    "relative_intensity": "10",
                    "unc_ri": "",
                    "multipolarity": "",
                },
            ]
        ),
        encoding="utf-8",
    )
    client = GammaClient(cache_dir=cache)
    result = client.query_gammas(
        Nuclide(z=63, a=151, isomer="m"),
        limit=10,
        excitation_kev=196.245,
        expected_half_life="58.9 us",
        min_energy_kev=0.0,
    )
    assert result.status == "relative_adopted_level"
    energies = [round(g.energy_kev, 2) for g in result.lines]
    assert energies == [174.7, 196.0]
    assert 196.54 not in energies
    assert 175.0 not in energies
    assert all(abs((g.parent_energy_kev or 0) - 196.245) < 0.05 for g in result.lines)
    assert all(g.parent_match_conflict for g in result.lines)
    assert all(g.normalization_status == "relative_adopted_level" for g in result.lines)


def test_select_best_energy_helper_prefers_closest() -> None:
    from photonuclear_catalog.gammas import _select_best_energy

    selected, conflict = _select_best_energy(
        [(77.2593, 3.8e-8), (78.2331, 1.65e-7)],
        excitation_kev=78.2331,
        expected_sec=1.65e-7,
    )
    assert selected == 78.2331
    assert conflict is True

    selected2, conflict2 = _select_best_energy(
        [(89.6128, None), (89.8496, 3.84e-7)],
        excitation_kev=89.8496,
        expected_sec=3.84e-7,
    )
    assert selected2 == 89.8496
    assert conflict2 is True
