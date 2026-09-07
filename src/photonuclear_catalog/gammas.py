"""Fetch top decay gammas via IAEA LiveChart (ENSDF; NuDat-compatible source)."""

from __future__ import annotations

import csv
import io
import json
import re
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import requests

from photonuclear_catalog.nuclide import Nuclide

LIVECHART_URL = "https://nds.iaea.org/relnsd/v1/data"

# LiveChart isomer endpoints like "207mpb" / "60mco" usually return empty.
# Decay radiation for isomers is listed under the ground-state nuclide id,
# distinguished by the parent level energy (p_energy).
_PARENT_ENERGY_TOL_KEV = 2.0
_HL_RATIO_TOL = 5.0  # accept HL within factor of 5 when Exc matches
_ENERGY_KEY_TOL_KEV = 0.05
# After choosing one parent/level, keep only lines from that exact energy.
_SELECTED_LEVEL_TOL_KEV = 0.05


@dataclass(frozen=True)
class GammaLine:
    product: str
    energy_kev: float
    intensity: float | None  # absolute-or-best Iγ; None if unknown (not 0.0)
    intensity_unc: float | None  # numeric uncertainty only
    decay_mode: str
    parent_energy_kev: float | None
    half_life: str
    rank: int | None  # None when intensity unknown (not a true strength rank)
    intensity_ensdf: float | None = None
    intensity_relative: float | None = None
    branch_percent: float | None = None
    normalization_status: str = "absolute"  # absolute | branch_relative | uncertain | relative_adopted_level
    uncertainty_status: str = ""  # e.g. branch_unc_not_propagated | intensity_missing
    status: str = "ok"
    note: str = ""
    parent_match_conflict: bool = False

    def to_dict(self) -> dict:
        d = asdict(self)
        # CSV-friendly: unknown intensities/ranks/energies as empty, not 0.0 / "nan".
        for key in (
            "intensity",
            "intensity_unc",
            "intensity_ensdf",
            "intensity_relative",
            "rank",
            "energy_kev",
            "parent_energy_kev",
        ):
            val = d.get(key)
            if val is None:
                d[key] = ""
            elif isinstance(val, float) and val != val:  # NaN
                d[key] = ""
        return d


@dataclass
class GammaQueryResult:
    lines: list[GammaLine]
    status: str
    detail: str = ""


class GammaClient:
    def __init__(
        self,
        cache_dir: Path | str | None = None,
        timeout_s: float = 30.0,
        pause_s: float = 0.15,
        max_retries: int = 4,
        session: requests.Session | None = None,
    ) -> None:
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout_s = timeout_s
        self.pause_s = pause_s
        self.max_retries = max_retries
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": "photonuclear-catalog/0.2"})

    def top_gammas(
        self,
        nuclide: Nuclide,
        limit: int = 10,
        excitation_kev: float | None = None,
        min_energy_kev: float = 0.0,
        expected_half_life: str | None = None,
        exclude_xrays: bool = True,
    ) -> list[GammaLine]:
        """Return strongest decay gammas for a ground state or isomer."""
        result = self.query_gammas(
            nuclide,
            limit=limit,
            excitation_kev=excitation_kev,
            min_energy_kev=min_energy_kev,
            expected_half_life=expected_half_life,
            exclude_xrays=exclude_xrays,
        )
        # relative_adopted_level: real γ lines exist (relative Iγ from levels).
        if result.status in {"ok", "relative_adopted_level"}:
            return result.lines
        return [
            GammaLine(
                product=nuclide.symbol,
                energy_kev=float("nan"),
                intensity=None,
                intensity_unc=None,
                decay_mode=result.status.upper(),
                parent_energy_kev=None,
                half_life="",
                rank=None,
                status=result.status,
                note=result.detail,
                uncertainty_status="status_row",
                normalization_status="not_applicable",
                parent_match_conflict="conflict" in (result.detail or "").lower(),
            )
        ]

    def query_gammas(
        self,
        nuclide: Nuclide,
        limit: int = 10,
        excitation_kev: float | None = None,
        min_energy_kev: float = 0.0,
        expected_half_life: str | None = None,
        exclude_xrays: bool = True,
    ) -> GammaQueryResult:
        """Query with explicit status codes for missing / filtered data."""
        try:
            rows = self._fetch_rows(nuclide, rad_types="g")
            xrows = self._fetch_rows(nuclide, rad_types="x") if exclude_xrays else []
        except Exception as exc:  # noqa: BLE001
            return GammaQueryResult([], "network_error", str(exc))

        xray_keys = {_line_key(r) for r in xrows if _line_key(r) is not None}

        want_isomer = bool(nuclide.isomer)
        expected_sec = (
            _half_life_to_seconds(expected_half_life) if expected_half_life else None
        )

        # Collect decay_rads candidates, then keep only the single best parent.
        loose_rows: list[dict] = []
        for row in rows:
            parent_e = _parse_parent_energy(row.get("p_energy"))
            row_hl_sec = _parse_float(row.get("half_life_sec"))
            if not _parent_matches(
                parent_e,
                row_hl_sec=row_hl_sec,
                want_isomer=want_isomer,
                excitation_kev=excitation_kev,
                expected_sec=expected_sec,
            ):
                continue
            energy = _parse_float(row.get("energy"))
            if energy is None or abs(energy - 511.0) < 0.3:
                continue
            loose_rows.append(row)

        selected_parent_e: float | None = None
        parent_conflict = False
        matched_rows: list[dict] = []

        if loose_rows and want_isomer and excitation_kev is not None and excitation_kev > 0:
            selected_parent_e, parent_conflict = _select_best_parent_energy(
                loose_rows,
                excitation_kev=excitation_kev,
                expected_sec=expected_sec,
            )
            source_rows = [
                row
                for row in loose_rows
                if (pe := _parse_parent_energy(row.get("p_energy"))) is not None
                and selected_parent_e is not None
                and abs(pe - selected_parent_e) <= _SELECTED_LEVEL_TOL_KEV
            ]
        else:
            source_rows = list(loose_rows)

        for row in source_rows:
            key = _line_key(row)
            if exclude_xrays and key is not None and _key_in_xray_set(key, xray_keys):
                continue
            energy = _parse_float(row.get("energy"))
            assert energy is not None
            if energy < min_energy_kev:
                continue
            matched_rows.append(row)

        if not matched_rows:
            # Isomers with no usable decay_rads → Adopted Levels relative Iγ.
            if want_isomer:
                adopted = self._query_adopted_level_gammas(
                    nuclide,
                    limit=limit,
                    excitation_kev=excitation_kev,
                    min_energy_kev=min_energy_kev,
                    expected_half_life=expected_half_life,
                )
                if adopted is not None:
                    return adopted

            any_nuclear = False
            any_cut = False
            for row in source_rows or loose_rows:
                energy = _parse_float(row.get("energy"))
                if energy is None or abs(energy - 511.0) < 0.3:
                    continue
                key = _line_key(row)
                if exclude_xrays and key is not None and _key_in_xray_set(key, xray_keys):
                    continue
                any_nuclear = True
                if energy < min_energy_kev:
                    any_cut = True
            if any_cut and not any(
                (_parse_float(r.get("energy")) or -1) >= min_energy_kev
                for r in (source_rows or loose_rows)
                if _line_key(r) is None
                or not exclude_xrays
                or not _key_in_xray_set(_line_key(r), xray_keys)  # type: ignore[arg-type]
            ):
                return GammaQueryResult(
                    [],
                    "cut_excluded",
                    f"Nuclear lines exist but were excluded by Eγ ≥ {min_energy_kev:g} keV cut.",
                )
            if loose_rows and not any_nuclear:
                return GammaQueryResult(
                    [],
                    "xray_only",
                    "Parent matched; no nuclear γ remained after atomic X-ray filter "
                    "(LiveChart rad_types=x).",
                )
            return GammaQueryResult(
                [],
                "no_matched_ensdf_parent",
                "No ENSDF parent level matched (excitation / half-life).",
            )

        # Per-(parent, decay) normalization: LiveChart Iγ is often already absolute.
        # Only multiply by BR/100 when the mode is clearly branch-relative
        # (any line intensity exceeds the branching percent).
        norm_by_mode = _detect_normalization(matched_rows)

        candidates: list[tuple[float, float, float | None, float | None, str, bool, float, dict, str]] = []
        for row in matched_rows:
            raw_intensity = str(row.get("intensity") or "").strip()
            if raw_intensity == "":
                intensity_ensdf = 0.0
                intensity_known = False
            else:
                intensity_ensdf = _parse_float(raw_intensity)
                if intensity_ensdf is None or intensity_ensdf <= 0:
                    continue
                intensity_known = True

            branch = _parse_float(row.get("decay_%"))
            mode_key = _mode_key(row)
            norm = norm_by_mode.get(mode_key, "absolute")
            note_parts: list[str] = []

            if norm == "branch_relative" and branch is not None and branch > 0:
                intensity_abs = intensity_ensdf * (branch / 100.0)
                note_parts.append(f"branch_relative→abs×BR/100; BR={branch:g}%")
                unc_num, unc_status = _scale_uncertainty_fields(
                    str(row.get("unc_i") or ""), branch / 100.0
                )
            elif norm == "uncertain":
                intensity_abs = intensity_ensdf
                note_parts.append("normalization uncertain; ENSDF Iγ kept as-is")
                unc_num, unc_status = _parse_unc_fields(str(row.get("unc_i") or ""))
            else:
                intensity_abs = intensity_ensdf
                unc_num, unc_status = _parse_unc_fields(str(row.get("unc_i") or ""))
                if branch is not None and 0 < branch < 99.5:
                    note_parts.append(f"absolute ENSDF Iγ (BR={branch:g}% not reapplied)")

            if parent_conflict:
                pe_txt = (
                    f"{selected_parent_e:g}"
                    if selected_parent_e is not None
                    else "?"
                )
                note_parts.append(
                    "parent_match_conflict: multiple parents within Exc tolerance; "
                    f"selected closest (p_energy={pe_txt} keV"
                    + (
                        f", NUBASE HL={expected_half_life}"
                        if expected_half_life
                        else ""
                    )
                    + ")"
                )

            if not intensity_known:
                intensity_abs = None
                unc_num = None
                unc_status = "intensity_missing"

            sort_key = intensity_abs if intensity_known and intensity_abs is not None else -1.0
            energy = _parse_float(row.get("energy"))
            assert energy is not None
            candidates.append(
                (
                    sort_key,
                    intensity_abs,
                    intensity_ensdf if intensity_known else None,
                    branch,
                    norm,
                    intensity_known,
                    energy,
                    row,
                    "; ".join(note_parts),
                    unc_num,
                    unc_status,
                )
            )

        if not candidates:
            return GammaQueryResult(
                [],
                "intensity_missing",
                "Parent matched but Iγ not given in ENSDF.",
            )

        # Rank only lines with known intensity (strongest first).
        known = [c for c in candidates if c[5]]
        unknown = [c for c in candidates if not c[5]]
        known.sort(key=lambda t: (-t[0], t[6]))
        unknown.sort(key=lambda t: t[6])

        out: list[GammaLine] = []
        used_energies: list[float] = []

        def _append_line(item: tuple, *, ranked: bool) -> None:
            (
                _sort,
                intensity_abs,
                intensity_ensdf,
                branch,
                norm,
                intensity_known,
                energy,
                row,
                note,
                unc_num,
                unc_status,
            ) = item
            if any(abs(energy - e) < 0.05 for e in used_energies):
                return
            used_energies.append(energy)
            status = "ok"
            if parent_conflict:
                status = "parent_match_conflict"
            elif not intensity_known:
                status = "intensity_missing"
            rank: int | None
            if ranked and intensity_known:
                rank = len([g for g in out if g.rank is not None]) + 1
            else:
                rank = None
            out.append(
                GammaLine(
                    product=nuclide.symbol,
                    energy_kev=energy,
                    intensity=intensity_abs,
                    intensity_unc=unc_num,
                    decay_mode=str(row.get("decay") or ""),
                    parent_energy_kev=_parse_parent_energy(row.get("p_energy")),
                    half_life=_row_half_life(row),
                    rank=rank,
                    intensity_ensdf=intensity_ensdf,
                    intensity_relative=(
                        intensity_ensdf if norm == "branch_relative" else None
                    ),
                    branch_percent=branch,
                    normalization_status=norm,
                    uncertainty_status=unc_status,
                    status=status,
                    note=note,
                    parent_match_conflict=parent_conflict,
                )
            )

        for item in known:
            if len([g for g in out if g.rank is not None]) >= limit:
                break
            _append_line(item, ranked=True)

        # Keep unknown-Iγ lines (no fake rank) only when few/no ranked lines exist.
        if len([g for g in out if g.rank is not None]) < limit:
            for item in unknown:
                if len(out) >= limit:
                    break
                _append_line(item, ranked=False)

        if not out:
            return GammaQueryResult([], "no_matched_ensdf_parent", "No gamma lines selected.")
        return GammaQueryResult(out, "ok" if not parent_conflict else "ok")

    def _query_adopted_level_gammas(
        self,
        nuclide: Nuclide,
        *,
        limit: int,
        excitation_kev: float | None,
        min_energy_kev: float,
        expected_half_life: str | None,
    ) -> GammaQueryResult | None:
        """Fallback for isomers missing from decay_rads (stable GS / higher isomer).

        LiveChart ``fields=gammas`` carries Adopted Levels transitions with
        *relative* Iγ. This is not an absolute decay-radiation dataset.
        """
        if not nuclide.isomer:
            return None
        if excitation_kev is None or excitation_kev <= 0:
            return None

        try:
            gamma_rows = self._fetch_field(nuclide, fields="gammas")
            level_rows = self._fetch_field(nuclide, fields="levels")
        except Exception as exc:  # noqa: BLE001
            return GammaQueryResult([], "network_error", str(exc))

        expected_sec = (
            _half_life_to_seconds(expected_half_life) if expected_half_life else None
        )
        level, level_conflict = _match_adopted_level(
            level_rows,
            excitation_kev=excitation_kev,
            expected_half_life=expected_half_life,
        )

        selected_level_e: float | None = None
        parent_conflict = level_conflict
        if level is not None:
            selected_level_e = _parse_float(level.get("energy"))

        if selected_level_e is None:
            # Levels table sparse: choose one start_level among γ rows near Exc.
            start_rows = [
                r
                for r in gamma_rows
                if (se := _parse_float(r.get("start_level_energy"))) is not None
                and abs(se - excitation_kev) <= _PARENT_ENERGY_TOL_KEV
            ]
            if not start_rows:
                return None
            selected_level_e, start_conflict = _select_best_energy(
                [
                    (_parse_float(r.get("start_level_energy")), None)
                    for r in start_rows
                ],
                excitation_kev=excitation_kev,
                expected_sec=expected_sec,
            )
            parent_conflict = parent_conflict or start_conflict

        if selected_level_e is None:
            return None

        matched: list[dict] = []
        for row in gamma_rows:
            start_e = _parse_float(row.get("start_level_energy"))
            if start_e is None:
                continue
            # Only gammas from the single selected level — not every start in ±2 keV.
            if abs(start_e - selected_level_e) > _SELECTED_LEVEL_TOL_KEV:
                continue
            energy = _parse_float(row.get("energy"))
            if energy is None or abs(energy - 511.0) < 0.3:
                continue
            if energy < min_energy_kev:
                continue
            matched.append(row)

        if not matched:
            return GammaQueryResult(
                [],
                "no_absolute_decay_gamma_dataset",
                "Adopted Levels parent matched, but no absolute decay_rads dataset "
                "and no level-scheme γ transitions listed.",
            )

        # Rank by relative intensity (not absolute %).
        scored: list[tuple[float, float, float | None, dict]] = []
        for row in matched:
            rel = _parse_float(row.get("relative_intensity"))
            energy = _parse_float(row.get("energy"))
            assert energy is not None
            sort_key = rel if rel is not None and rel > 0 else -1.0
            scored.append((sort_key, energy, rel, row))
        scored.sort(key=lambda t: (-t[0], t[1]))

        hl = ""
        parent_e = selected_level_e
        if level is not None:
            hl = _level_half_life(level)

        out: list[GammaLine] = []
        used: list[float] = []
        for sort_key, energy, rel, row in scored:
            if len(out) >= limit:
                break
            if any(abs(energy - e) < 0.05 for e in used):
                continue
            used.append(energy)
            known = rel is not None and rel > 0
            unc_num, unc_status = _parse_unc_fields(str(row.get("unc_ri") or ""))
            multipole = str(row.get("multipolarity") or "").strip()
            note_parts = [
                "relative_adopted_level; not absolute Iγ% (IT/level scheme; ICC may apply)"
            ]
            if multipole:
                note_parts.append(f"mult={multipole}")
            if parent_conflict:
                note_parts.append(
                    "parent_match_conflict: multiple Adopted Levels near Exc; "
                    f"selected closest ({selected_level_e:g} keV)"
                )
            status = (
                "parent_match_conflict"
                if parent_conflict
                else "relative_adopted_level"
            )
            out.append(
                GammaLine(
                    product=nuclide.symbol,
                    energy_kev=energy,
                    intensity=rel if known else None,
                    intensity_unc=unc_num if known else None,
                    decay_mode="IT",
                    parent_energy_kev=parent_e,
                    half_life=hl or (expected_half_life or ""),
                    rank=(len(out) + 1) if known else None,
                    intensity_ensdf=rel if known else None,
                    intensity_relative=rel if known else None,
                    branch_percent=100.0,
                    normalization_status="relative_adopted_level",
                    uncertainty_status=unc_status if known else "intensity_missing",
                    status=status,
                    note="; ".join(note_parts),
                    parent_match_conflict=parent_conflict,
                )
            )

        if not out:
            return GammaQueryResult(
                [],
                "no_absolute_decay_gamma_dataset",
                "Adopted Levels parent matched; no usable relative γ intensities.",
            )
        return GammaQueryResult(
            out,
            "relative_adopted_level",
            "decay_rads parent missing; relative Iγ from LiveChart fields=gammas "
            f"(start_level = {selected_level_e:g} keV).",
        )

    def _fetch_rows(self, nuclide: Nuclide, rad_types: str = "g") -> list[dict]:
        return self._fetch_field(nuclide, fields="decay_rads", rad_types=rad_types)

    def _fetch_field(
        self,
        nuclide: Nuclide,
        *,
        fields: str,
        rad_types: str | None = None,
    ) -> list[dict]:
        gs_id = Nuclide(z=nuclide.z, a=nuclide.a, isomer="").livechart_id
        cache_path = None
        if self.cache_dir:
            if fields == "decay_rads":
                suffix = "" if rad_types == "g" else f"_{rad_types}"
            else:
                suffix = f"_{fields}"
            cache_path = self.cache_dir / f"{gs_id}{suffix}.json"
            if cache_path.exists():
                return json.loads(cache_path.read_text(encoding="utf-8"))

        params: dict[str, str] = {
            "fields": fields,
            "nuclides": gs_id,
        }
        if fields == "decay_rads":
            params["rad_types"] = rad_types or "g"
        last_error: Exception | None = None
        rows: list[dict] = []
        for attempt in range(self.max_retries):
            time.sleep(self.pause_s * (attempt + 1))
            try:
                resp = self.session.get(
                    LIVECHART_URL, params=params, timeout=self.timeout_s
                )
                if resp.status_code in {404}:
                    rows = []
                    break
                if resp.status_code in {429, 502, 503, 504}:
                    last_error = requests.HTTPError(
                        f"{resp.status_code} for {resp.url}", response=resp
                    )
                    continue
                resp.raise_for_status()
                text = resp.text
                if (
                    not text.strip()
                    or text.strip() == "0"
                    or text.lower().startswith("<!doctype")
                    or "<html" in text.lower()
                ):
                    # LiveChart returns literal "0" for empty decay_rads.
                    if text.strip() == "0":
                        rows = []
                        last_error = None
                        break
                    last_error = RuntimeError(
                        f"Empty or HTML response for {gs_id}/{fields} "
                        f"(status={resp.status_code})"
                    )
                    continue
                reader = csv.DictReader(io.StringIO(text))
                rows = list(reader)
                last_error = None
                break
            except requests.RequestException as exc:
                last_error = exc
        if last_error is not None and not rows:
            raise last_error

        if cache_path is not None:
            cache_path.write_text(json.dumps(rows), encoding="utf-8")
        return rows


def _match_adopted_level(
    level_rows: list[dict],
    *,
    excitation_kev: float,
    expected_half_life: str | None,
) -> tuple[dict | None, bool]:
    """Pick the single closest Adopted Levels row near Exc (optionally HL-aided).

    Returns ``(level_row_or_None, conflict)``. ``conflict`` is True when more
    than one distinct level energy lies inside the tolerance window.
    """
    expected_sec = (
        _half_life_to_seconds(expected_half_life) if expected_half_life else None
    )
    candidates: list[dict] = []
    for row in level_rows:
        energy = _parse_float(row.get("energy"))
        if energy is None:
            continue
        if abs(energy - excitation_kev) <= _PARENT_ENERGY_TOL_KEV:
            candidates.append(row)
    if not candidates:
        return None, False

    unique_energies = {
        round(e, 4)
        for e in (_parse_float(r.get("energy")) for r in candidates)
        if e is not None
    }
    conflict = len(unique_energies) > 1

    pool = candidates
    if expected_sec is not None:
        hl_ok = [
            row
            for row in candidates
            if (hl := _parse_float(row.get("half_life_sec")))
            and _hl_compatible(expected_sec, hl)
        ]
        if hl_ok:
            pool = hl_ok

    def _exc_delta(row: dict) -> float:
        energy = _parse_float(row.get("energy"))
        assert energy is not None
        return abs(energy - excitation_kev)

    best = min(pool, key=_exc_delta)
    return best, conflict


def _select_best_parent_energy(
    rows: list[dict],
    *,
    excitation_kev: float,
    expected_sec: float | None,
) -> tuple[float | None, bool]:
    """Among decay_rads parents in the candidate window, pick one closest Exc."""
    pairs: list[tuple[float | None, float | None]] = []
    for row in rows:
        pairs.append(
            (
                _parse_parent_energy(row.get("p_energy")),
                _parse_float(row.get("half_life_sec")),
            )
        )
    return _select_best_energy(
        pairs, excitation_kev=excitation_kev, expected_sec=expected_sec
    )


def _select_best_energy(
    pairs: list[tuple[float | None, float | None]],
    *,
    excitation_kev: float,
    expected_sec: float | None,
) -> tuple[float | None, bool]:
    """Select one energy from ``(energy, half_life_sec)`` pairs.

    Energy tolerance is assumed to have already filtered the input. When several
    distinct energies remain, prefer HL-compatible ones (if available), then the
    smallest |E − Exc|. ``conflict`` is True iff more than one distinct energy
    was present before selection.
    """
    by_energy: dict[float, float | None] = {}
    for energy, hl_sec in pairs:
        if energy is None:
            continue
        key = round(energy, 4)
        if key not in by_energy:
            by_energy[key] = hl_sec

    if not by_energy:
        return None, False

    conflict = len(by_energy) > 1
    candidates = list(by_energy.items())
    if expected_sec is not None:
        hl_ok = [
            (e, hl)
            for e, hl in candidates
            if hl is not None and _hl_compatible(expected_sec, hl)
        ]
        if hl_ok:
            candidates = hl_ok

    best_e = min(candidates, key=lambda t: abs(t[0] - excitation_kev))[0]
    return best_e, conflict


def _level_half_life(row: dict) -> str:
    hl = str(row.get("half_life") or "").strip()
    unit = str(row.get("unit_hl") or "").strip()
    if hl and unit:
        return f"{hl} {unit}"
    return hl


def _line_key(row: dict) -> tuple[str, str, str, float] | None:
    """Unique key (parent ZA-ish, p_energy, decay, E) for g/x set difference."""
    energy = _parse_float(row.get("energy"))
    if energy is None:
        return None
    parent = f"{row.get('p_z') or ''}|{row.get('p_n') or ''}|{row.get('p_symbol') or ''}"
    pe = str(row.get("p_energy") or "").strip()
    decay = str(row.get("decay") or "").strip()
    return (parent, pe, decay, round(energy, 3))


def _key_in_xray_set(
    key: tuple[str, str, str, float],
    xray_keys: set[tuple[str, str, str, float]],
) -> bool:
    if key in xray_keys:
        return True
    # Tolerance match on energy with same parent/decay.
    parent, pe, decay, energy = key
    for xp, xpe, xdec, xe in xray_keys:
        if xp == parent and xpe == pe and xdec == decay and abs(energy - xe) <= _ENERGY_KEY_TOL_KEV:
            return True
    return False


def _mode_key(row: dict) -> tuple[str, str]:
    return (str(row.get("p_energy") or "").strip(), str(row.get("decay") or "").strip())


def _detect_normalization(rows: list[dict]) -> dict[tuple[str, str], str]:
    """Classify each (p_energy, decay) group as absolute / branch_relative / uncertain.

    LiveChart often already stores absolute Iγ (per 100 parent decays). A clear
    signature of branch-relative intensities is any line with Iγ > BR%.
    """
    by_mode: dict[tuple[str, str], list[tuple[float | None, float | None]]] = defaultdict(list)
    for row in rows:
        i = _parse_float(row.get("intensity"))
        br = _parse_float(row.get("decay_%"))
        by_mode[_mode_key(row)].append((i, br))

    out: dict[tuple[str, str], str] = {}
    for mode, pairs in by_mode.items():
        intensities = [i for i, _br in pairs if i is not None and i > 0]
        branches = [br for _i, br in pairs if br is not None and br > 0]
        if not intensities:
            out[mode] = "uncertain"
            continue
        br = branches[0] if branches else None
        if br is None or br >= 99.5:
            out[mode] = "absolute"
            continue
        mx = max(intensities)
        if mx > br + 1e-9:
            out[mode] = "branch_relative"
        else:
            out[mode] = "absolute"
    return out


def _scale_uncertainty(unc: str, factor: float) -> str:
    """Legacy helper kept for tests; prefer _scale_uncertainty_fields."""
    num, status = _scale_uncertainty_fields(unc, factor)
    if num is None:
        return status
    return f"{num:g} ({status})" if status else f"{num:g}"


def _parse_unc_fields(unc: str) -> tuple[float | None, str]:
    text = (unc or "").strip()
    if not text:
        return None, ""
    try:
        return float(text), ""
    except ValueError:
        return None, text


def _scale_uncertainty_fields(
    unc: str, factor: float
) -> tuple[float | None, str]:
    """Return (numeric_unc, uncertainty_status) when intensity is scaled."""
    text = (unc or "").strip()
    if not text:
        return None, "branch_unc_not_propagated"
    try:
        val = float(text)
    except ValueError:
        return None, f"branch_unc_not_propagated; raw={text}"
    return val * factor, "branch_unc_not_propagated"


def _row_half_life(row: dict) -> str:
    return " ".join(
        str(row.get(k) or "").strip()
        for k in ("half_life", "unit_hl")
        if str(row.get(k) or "").strip()
    )


def _parse_parent_energy(value: object) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _parse_float(value: object) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        out = float(text)
    except ValueError:
        return None
    if out != out:  # NaN
        return None
    return out


def _parent_energy_matches(
    parent_e: float | None,
    *,
    want_isomer: bool,
    excitation_kev: float | None,
) -> bool:
    """Backward-compatible Exc-only matcher (tests)."""
    return _parent_matches(
        parent_e,
        row_hl_sec=None,
        want_isomer=want_isomer,
        excitation_kev=excitation_kev,
        expected_sec=None,
    )


def _parent_matches(
    parent_e: float | None,
    *,
    row_hl_sec: float | None,
    want_isomer: bool,
    excitation_kev: float | None,
    expected_sec: float | None,
) -> bool:
    if not want_isomer:
        return parent_e is None or parent_e < 1.0

    if excitation_kev is not None and excitation_kev > 0:
        if parent_e is None:
            return False
        if abs(parent_e - excitation_kev) <= _PARENT_ENERGY_TOL_KEV:
            return True
        if abs(parent_e - excitation_kev) <= 20.0 and expected_sec and row_hl_sec:
            return _hl_compatible(expected_sec, row_hl_sec)
        return False

    if parent_e is None or parent_e < 1.0:
        return False
    if expected_sec and row_hl_sec:
        return _hl_compatible(expected_sec, row_hl_sec)
    return True


def _hl_compatible(a: float, b: float) -> bool:
    if a <= 0 or b <= 0:
        return False
    return max(a, b) / min(a, b) <= _HL_RATIO_TOL


_HL_UNIT_SEC = {
    "ys": 1e-24,
    "zs": 1e-21,
    "as": 1e-18,
    "fs": 1e-15,
    "ps": 1e-12,
    "ns": 1e-9,
    "us": 1e-6,
    "µs": 1e-6,
    "ms": 1e-3,
    "s": 1.0,
    "m": 60.0,
    "h": 3600.0,
    "d": 86400.0,
    "y": 365.25 * 86400.0,
}


def _half_life_to_seconds(text: str | None) -> float | None:
    """Parse NUBASE-style half-life strings, including limits like '<500 ms'."""
    if not text:
        return None
    raw = text.strip().lower()
    if not raw or raw in {"stable", "p-unstable", "?"}:
        return None
    raw = raw.replace("#", "")
    raw = re.sub(r"^[<>~≈]+", "", raw).strip()
    m = re.match(
        r"^([0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)\s*([a-zµμ]+)?",
        raw,
    )
    if not m:
        return None
    value = float(m.group(1))
    unit = (m.group(2) or "s").replace("μ", "u").replace("µ", "u")
    factor = _HL_UNIT_SEC.get(unit)
    if factor is None:
        return None
    return value * factor
