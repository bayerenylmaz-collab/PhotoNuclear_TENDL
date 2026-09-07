"""TENDL-2023 photonuclear residual cross-section helper."""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import requests

from photonuclear_catalog.nuclide import z_to_element

TENDL_BASE = "https://tendl.imperial.ac.uk/tendl_2023/gamma_file"

# Prefer s60: explicit MF3/MF10 channels tabulated to 60 MeV (needed for Emax~44).
# Fallback: endfmt (mt200) then default s30 ENDF.
_TENDL_VARIANTS: tuple[tuple[str, str, str], ...] = (
    ("s60", "endfs60", "{tag}-s60"),  # lib/endfs60/g-Au197-s60.tendl
    ("mt200", "endfmt", "{tag}-mt"),  # lib/endfmt/g-Au197-mt.tendl
    ("s30", "endf", "{tag}"),  # lib/endf/g-Au197.tendl
)

# Free-nucleon (n+p only) ENDF MT map — TEFAL/TALYS Table 4.2 / ENDF-6 Table 7.
# Clustered channels (d, t, h, α, …) are intentionally absent so they cannot
# silently replace (γ,xn yp) when MF10 is keyed only by IZAP+LFS.
_FREE_NUCLEON_MT: dict[tuple[int, int], int] = {
    (1, 0): 4,  # (z,n)
    (2, 0): 16,  # (z,2n)
    (3, 0): 17,  # (z,3n)
    (4, 0): 37,  # (z,4n)
    (5, 0): 152,  # (z,5n)
    (6, 0): 153,  # (z,6n)
    (7, 0): 160,  # (z,7n)
    (8, 0): 161,  # (z,8n)
    (0, 1): 103,  # (z,p)
    (0, 2): 111,  # (z,2p)
    (1, 1): 28,  # (z,np)
    (2, 1): 41,  # (z,2np)
    (3, 1): 42,  # (z,3np)
    (4, 1): 156,  # (z,4np)
    (5, 1): 162,  # (z,5np)
    (6, 1): 163,  # (z,6np)
    (7, 1): 164,  # (z,7np)
    (1, 2): 44,  # (z,n2p)
    (2, 2): 190,  # (z,2n2p)
    (3, 2): 179,  # (z,3n2p)
}

_EXCITATION_MATCH_TOL_KEV = 5.0


def free_nucleon_mt(x_neutrons: int, y_protons: int) -> int | None:
    """Return ENDF MT for free-nucleon (γ, xn yp), or None if unmapped/clustered."""
    if x_neutrons < 0 or y_protons < 0:
        return None
    if x_neutrons == 0 and y_protons == 0:
        return None
    return _FREE_NUCLEON_MT.get((x_neutrons, y_protons))


@dataclass(frozen=True)
class SigmaPoint:
    energy_mev: float
    sigma_b: float


@dataclass
class ResidualSigma:
    z: int
    a: int
    level_index: int  # ENDF LFS
    points: list[SigmaPoint]
    qm_ev: float | None = None
    qi_ev: float | None = None
    excitation_kev: float | None = None  # (QM − QI) / 1000

    def max_in_window(
        self, e_lo_mev: float, e_hi_mev: float
    ) -> tuple[float | None, float | None]:
        """Return (sigma_mb, e_at_sigma_max_mev) on [e_lo, e_hi].

        Piecewise-linear max is attained at tabulated vertices or the window
        endpoints — evaluate endpoints via ENDF linear interpolation.
        """
        return _max_interp_mb(self.points, e_lo_mev, e_hi_mev)

    def max_sigma_mb(self, e_lo_mev: float, e_hi_mev: float) -> float | None:
        mb, _e = self.max_in_window(e_lo_mev, e_hi_mev)
        return mb


@dataclass(frozen=True)
class SigmaLookup:
    sigma_max_mb: float | None
    e_at_sigma_max_mev: float | None
    status: str
    scope: str  # state_specific | ground_state | channel_total | missing
    lfs: int | None = None
    mt: int | None = None
    sigma_ground_state_mb: float | None = None
    sigma_channel_total_mb: float | None = None  # MF3(MT)
    e_at_channel_total_max_mev: float | None = None
    sigma_sum_mf10_lfs_mb: float | None = None  # Σ_LFS MF10 at same MT+IZAP
    mf10_coverage_ratio: float | None = None  # ΣMF10 / MF3
    e_at_total_max_mev: float | None = None
    lfs_excitation_kev: float | None = None
    sigma_qi_mev: float | None = None
    source_variant: str | None = None
    source_url: str | None = None
    source_sha256: str | None = None


class TendlPhotonuclear:
    """Lazy TENDL-2023 ENDF reader: MF3 channel totals + MF10 (MT,IZAP,LFS)."""

    def __init__(
        self,
        cache_dir: Path | str | None = None,
        *,
        prefer_variant: str = "s60",
    ) -> None:
        self.cache_dir = Path(cache_dir) if cache_dir else Path("data/tendl_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.prefer_variant = prefer_variant
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "photonuclear-catalog/0.2"})
        # target -> {(mt, z, a): {lfs: ResidualSigma}}
        self._mf10: dict[
            tuple[int, int], dict[tuple[int, int, int], dict[int, ResidualSigma]]
        ] = {}
        # target -> {mt: ResidualSigma}  (MF3 channel total)
        self._mf3: dict[tuple[int, int], dict[int, ResidualSigma]] = {}
        # target -> provenance
        self._source: dict[tuple[int, int], dict[str, str]] = {}

    def source_info(self, z: int, a: int) -> dict[str, str]:
        self.ensure_target(z, a)
        return dict(self._source.get((z, a), {}))

    def ensure_target(self, z: int, a: int) -> bool:
        key = (z, a)
        if key in self._mf10:
            return bool(self._mf10[key]) or bool(self._mf3.get(key))
        downloaded = self._download(z, a)
        if downloaded is None:
            self._mf10[key] = {}
            self._mf3[key] = {}
            self._source[key] = {}
            return False
        path, variant, url, sha = downloaded
        mf10, mf3 = _parse_endf_mt_file(path)
        self._mf10[key] = mf10
        self._mf3[key] = mf3
        self._source[key] = {
            "variant": variant,
            "url": url,
            "sha256": sha,
            "path": str(path),
        }
        return bool(mf10) or bool(mf3)

    def _provenance(self, target_z: int, target_a: int) -> dict[str, str | None]:
        src = self._source.get((target_z, target_a), {})
        return {
            "source_variant": src.get("variant"),
            "source_url": src.get("url"),
            "source_sha256": src.get("sha256"),
        }

    def residual_lookup(
        self,
        target_z: int,
        target_a: int,
        residual_z: int,
        residual_a: int,
        e_lo_mev: float,
        e_hi_mev: float,
        *,
        mt: int | None = None,
        x_neutrons: int | None = None,
        y_protons: int | None = None,
        is_metastable: bool = False,
        isomer_ordinal: int = 0,
        excitation_kev: float | None = None,
        channel_e_lo_mev: float | None = None,
    ) -> SigmaLookup:
        """Look up residual σ keyed by free-nucleon MT + IZAP + LFS.

        Prefers the TENDL **s60** evaluation (explicit channels to 60 MeV).
        LFS=0 is ground-state production for that MT+IZAP.
        Channel totals (MF3, ΣMF10) use ``channel_e_lo_mev`` (default: same as
        ``e_lo_mev``); isomer Eth must *not* cut channel references.
        Metastable LFS is matched by QI excitation vs NUBASE Exc.
        """
        if mt is None:
            if x_neutrons is None or y_protons is None:
                return SigmaLookup(None, None, "mt_unspecified", "missing")
            mt = free_nucleon_mt(x_neutrons, y_protons)
        if mt is None:
            return SigmaLookup(None, None, "mt_unmapped", "missing")

        ok = self.ensure_target(target_z, target_a)
        prov = self._provenance(target_z, target_a)
        if not ok:
            return SigmaLookup(
                None, None, "tendl_missing", "missing", mt=mt, **prov  # type: ignore[arg-type]
            )

        chan_lo = e_lo_mev if channel_e_lo_mev is None else channel_e_lo_mev
        mf3_curve = self._mf3.get((target_z, target_a), {}).get(mt)
        channel_total_mb, chan_e = (
            (None, None)
            if mf3_curve is None
            else mf3_curve.max_in_window(chan_lo, e_hi_mev)
        )

        by_level = self._mf10.get((target_z, target_a), {}).get(
            (mt, residual_z, residual_a)
        )
        if not by_level:
            # MF3 may still exist: isotope-level production without state split.
            if mf3_curve is not None and channel_total_mb is not None:
                if is_metastable:
                    return SigmaLookup(
                        None,
                        None,
                        "state_specific_missing",
                        "missing",
                        mt=mt,
                        sigma_channel_total_mb=channel_total_mb,
                        e_at_channel_total_max_mev=chan_e,
                        **prov,  # type: ignore[arg-type]
                    )
                return SigmaLookup(
                    None,
                    None,
                    "channel_total_only",
                    "channel_total",
                    mt=mt,
                    sigma_channel_total_mb=channel_total_mb,
                    e_at_channel_total_max_mev=chan_e,
                    **prov,  # type: ignore[arg-type]
                )
            return SigmaLookup(
                None,
                None,
                "mt_not_in_tendl",
                "missing",
                mt=mt,
                **prov,  # type: ignore[arg-type]
            )

        gs_mb, gs_e = (None, None)
        if 0 in by_level:
            gs_mb, gs_e = by_level[0].max_in_window(chan_lo, e_hi_mev)
        sum_mb, sum_e = _sum_lfs_max(by_level, chan_lo, e_hi_mev)
        coverage = _coverage_ratio(sum_mb, channel_total_mb)

        def _qi_mev(lfs: int | None) -> float | None:
            if lfs is None or lfs not in by_level:
                return None
            qi = by_level[lfs].qi_ev
            return None if qi is None else qi / 1.0e6

        if not is_metastable:
            if 0 not in by_level:
                return SigmaLookup(
                    None,
                    None,
                    "state_specific_missing",
                    "missing",
                    mt=mt,
                    sigma_ground_state_mb=None,
                    sigma_channel_total_mb=channel_total_mb,
                    e_at_channel_total_max_mev=chan_e,
                    sigma_sum_mf10_lfs_mb=sum_mb,
                    mf10_coverage_ratio=coverage,
                    e_at_total_max_mev=sum_e,
                    **prov,  # type: ignore[arg-type]
                )
            if gs_mb is None:
                return SigmaLookup(
                    None,
                    None,
                    "no_sigma_in_window",
                    "ground_state",
                    lfs=0,
                    mt=mt,
                    sigma_ground_state_mb=None,
                    sigma_channel_total_mb=channel_total_mb,
                    e_at_channel_total_max_mev=chan_e,
                    sigma_sum_mf10_lfs_mb=sum_mb,
                    mf10_coverage_ratio=coverage,
                    e_at_total_max_mev=sum_e,
                    lfs_excitation_kev=0.0,
                    sigma_qi_mev=_qi_mev(0),
                    **prov,  # type: ignore[arg-type]
                )
            return SigmaLookup(
                gs_mb,
                gs_e,
                "ok",
                "ground_state",
                lfs=0,
                mt=mt,
                sigma_ground_state_mb=gs_mb,
                sigma_channel_total_mb=channel_total_mb,
                e_at_channel_total_max_mev=chan_e,
                sigma_sum_mf10_lfs_mb=sum_mb,
                mf10_coverage_ratio=coverage,
                e_at_total_max_mev=sum_e,
                lfs_excitation_kev=0.0,
                sigma_qi_mev=_qi_mev(0),
                **prov,  # type: ignore[arg-type]
            )

        lfs = _match_excited_lfs(by_level, excitation_kev, isomer_ordinal)
        if lfs is None:
            return SigmaLookup(
                None,
                None,
                "state_specific_missing",
                "missing",
                mt=mt,
                sigma_ground_state_mb=gs_mb,
                sigma_channel_total_mb=channel_total_mb,
                e_at_channel_total_max_mev=chan_e,
                sigma_sum_mf10_lfs_mb=sum_mb,
                mf10_coverage_ratio=coverage,
                e_at_total_max_mev=sum_e,
                **prov,  # type: ignore[arg-type]
            )

        mb, e_at = by_level[lfs].max_in_window(e_lo_mev, e_hi_mev)
        ex = by_level[lfs].excitation_kev
        if mb is None:
            return SigmaLookup(
                None,
                None,
                "no_sigma_in_window",
                "state_specific",
                lfs=lfs,
                mt=mt,
                sigma_ground_state_mb=gs_mb,
                sigma_channel_total_mb=channel_total_mb,
                e_at_channel_total_max_mev=chan_e,
                sigma_sum_mf10_lfs_mb=sum_mb,
                mf10_coverage_ratio=coverage,
                e_at_total_max_mev=sum_e,
                lfs_excitation_kev=ex,
                sigma_qi_mev=_qi_mev(lfs),
                **prov,  # type: ignore[arg-type]
            )
        return SigmaLookup(
            mb,
            e_at,
            "ok",
            "state_specific",
            lfs=lfs,
            mt=mt,
            sigma_ground_state_mb=gs_mb,
            sigma_channel_total_mb=channel_total_mb,
            e_at_channel_total_max_mev=chan_e,
            sigma_sum_mf10_lfs_mb=sum_mb,
            mf10_coverage_ratio=coverage,
            e_at_total_max_mev=sum_e,
            lfs_excitation_kev=ex,
            sigma_qi_mev=_qi_mev(lfs),
            **prov,  # type: ignore[arg-type]
        )

    def residual_max_mb(
        self,
        target_z: int,
        target_a: int,
        residual_z: int,
        residual_a: int,
        e_lo_mev: float,
        e_hi_mev: float,
        prefer_level: int | None = None,
        *,
        mt: int | None = None,
        x_neutrons: int | None = None,
        y_protons: int | None = None,
    ) -> tuple[float | None, str]:
        """Backward-compatible wrapper."""
        lookup = self.residual_lookup(
            target_z,
            target_a,
            residual_z,
            residual_a,
            e_lo_mev,
            e_hi_mev,
            mt=mt,
            x_neutrons=x_neutrons,
            y_protons=y_protons,
            is_metastable=prefer_level is not None and prefer_level > 0,
            isomer_ordinal=prefer_level or 0,
        )
        return lookup.sigma_max_mb, lookup.status

    def _download(self, z: int, a: int) -> tuple[Path, str, str, str] | None:
        """Download preferred TENDL variant; return (path, variant, url, sha256)."""
        el = z_to_element(z)
        tag = f"{el}{a:03d}"
        ordered = sorted(
            _TENDL_VARIANTS,
            key=lambda v: 0 if v[0] == self.prefer_variant else 1,
        )
        for variant, lib, suffix_template in ordered:
            # suffix_template uses {tag}=Au197 → file g-Au197-s60.tendl
            filename = f"g-{suffix_template.format(tag=tag)}.tendl"
            dest = self.cache_dir / filename
            url = f"{TENDL_BASE}/{el}/{tag}/lib/{lib}/{filename}"
            if dest.exists() and dest.stat().st_size > 1000:
                sha = hashlib.sha256(dest.read_bytes()).hexdigest()
                return dest, variant, url, sha
            try:
                resp = self.session.get(url, timeout=180)
                if resp.status_code != 200 or len(resp.content) < 1000:
                    continue
                dest.write_bytes(resp.content)
                sha = hashlib.sha256(resp.content).hexdigest()
                return dest, variant, url, sha
            except requests.RequestException:
                continue
        return None


def _match_excited_lfs(
    by_level: dict[int, ResidualSigma],
    excitation_kev: float | None,
    isomer_ordinal: int,
) -> int | None:
    """Match metastable to LFS>0 via QI excitation; ordinal only as last resort."""
    excited = sorted(lfs for lfs in by_level if lfs != 0)
    if not excited:
        return None
    if excitation_kev is not None and excitation_kev > 0:
        best: int | None = None
        best_d = None
        for lfs in excited:
            ex = by_level[lfs].excitation_kev
            if ex is None:
                continue
            d = abs(ex - excitation_kev)
            if d <= _EXCITATION_MATCH_TOL_KEV and (best_d is None or d < best_d):
                best, best_d = lfs, d
        if best is not None:
            return best
        # No excitation match → do not invent an LFS.
        return None
    idx = max(1, isomer_ordinal) - 1
    if idx < len(excited):
        return excited[idx]
    return None


def _coverage_ratio(
    sum_mf10_mb: float | None, channel_total_mb: float | None
) -> float | None:
    if sum_mf10_mb is None or channel_total_mb is None:
        return None
    if channel_total_mb <= 0:
        return None
    return round(sum_mf10_mb / channel_total_mb, 6)


def _parse_endf_mt_file(
    path: Path,
) -> tuple[
    dict[tuple[int, int, int], dict[int, ResidualSigma]],
    dict[int, ResidualSigma],
]:
    """Parse MF10 residuals and MF3 channel totals from a TENDL mt file."""
    mf10: dict[tuple[int, int, int], dict[int, ResidualSigma]] = {}
    mf3: dict[int, ResidualSigma] = {}
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if len(line) < 75:
            i += 1
            continue
        try:
            mf = int(line[70:72])
            mt = int(line[72:75])
        except ValueError:
            i += 1
            continue
        if mf == 3 and mt != 0:
            i, curve = _parse_mf3_section(lines, i, mt)
            if curve is not None:
                mf3[mt] = curve
            continue
        if mf != 10 or mt == 0:
            i += 1
            continue
        i += 1
        while i < len(lines):
            ln = lines[i]
            if len(ln) < 75:
                i += 1
                continue
            try:
                mf2 = int(ln[70:72])
                mt2 = int(ln[72:75])
            except ValueError:
                i += 1
                continue
            if mf2 != 10 or mt2 != mt:
                break
            vals = _endf_floats(ln)
            if not _is_izap_header(vals):
                i += 1
                continue
            izap = int(round(vals[2]))  # type: ignore[arg-type]
            lfs = int(round(vals[3] or 0))
            qm, qi = vals[0], vals[1]
            i += 1
            points: list[SigmaPoint] = []
            while i < len(lines):
                ln2 = lines[i]
                if len(ln2) < 75:
                    i += 1
                    continue
                try:
                    mf3x = int(ln2[70:72])
                    mt3 = int(ln2[72:75])
                except ValueError:
                    i += 1
                    continue
                if mf3x != 10 or mt3 != mt:
                    break
                vals2 = _endf_floats(ln2)
                if _is_izap_header(vals2):
                    break
                clean = [v for v in vals2 if v is not None]
                j = 0
                while j + 1 < len(clean):
                    e, s = clean[j], clean[j + 1]
                    if 1e5 <= e <= 5e8 and 0 <= s < 1e3:
                        points.append(SigmaPoint(energy_mev=e / 1e6, sigma_b=s))
                        j += 2
                    else:
                        j += 1
                i += 1
            if not points:
                continue
            z = izap // 1000
            a = izap % 1000
            ex = None
            if qm is not None and qi is not None:
                ex = (qm - qi) / 1000.0  # eV → keV
            key = (mt, z, a)
            bucket = mf10.setdefault(key, {})
            # Same MT+IZAP+LFS should be unique; keep denser table if duplicated.
            prev = bucket.get(lfs)
            if prev is None or len(points) >= len(prev.points):
                bucket[lfs] = ResidualSigma(
                    z=z,
                    a=a,
                    level_index=lfs,
                    points=points,
                    qm_ev=qm,
                    qi_ev=qi,
                    excitation_kev=ex,
                )
    return mf10, mf3


def _parse_mf3_section(
    lines: list[str], start: int, mt: int
) -> tuple[int, ResidualSigma | None]:
    """Parse one MF3/MT section into an energy–σ curve; return (next_index, curve)."""
    i = start + 1
    points: list[SigmaPoint] = []
    while i < len(lines):
        ln = lines[i]
        if len(ln) < 75:
            i += 1
            continue
        try:
            mf2 = int(ln[70:72])
            mt2 = int(ln[72:75])
        except ValueError:
            i += 1
            continue
        if mf2 != 3 or mt2 != mt:
            break
        vals = _endf_floats(ln)
        clean = [v for v in vals if v is not None]
        j = 0
        while j + 1 < len(clean):
            e, s = clean[j], clean[j + 1]
            # MF3 σ in barns; photonuclear channels are typically ≪ 10 b.
            if 1e5 <= e <= 5e8 and 0 <= s < 10:
                points.append(SigmaPoint(energy_mev=e / 1e6, sigma_b=s))
                j += 2
            else:
                j += 1
        i += 1
    if not points:
        return i, None
    by_e = {p.energy_mev: p.sigma_b for p in points}
    ordered = [SigmaPoint(e, by_e[e]) for e in sorted(by_e)]
    return i, ResidualSigma(z=0, a=0, level_index=0, points=ordered)


def _is_izap_header(vals: list[float | None]) -> bool:
    if len(vals) < 4:
        return False
    qm, _qi, izap, _lfs = vals[0], vals[1], vals[2], vals[3]
    if izap is None or qm is None:
        return False
    if not (1000 <= izap <= 120000 and abs(izap - round(izap)) < 1e-6):
        return False
    return qm < 0


def _interp_sigma_b_points(points: list[SigmaPoint], energy_mev: float) -> float:
    """Linear interpolate σ(b); 0 outside the tabulated range."""
    if not points:
        return 0.0
    if energy_mev < points[0].energy_mev - 1e-12 or energy_mev > points[-1].energy_mev + 1e-12:
        return 0.0
    if abs(energy_mev - points[0].energy_mev) < 1e-12:
        return points[0].sigma_b
    for i in range(1, len(points)):
        e0, e1 = points[i - 1].energy_mev, points[i].energy_mev
        if e0 <= energy_mev <= e1:
            if e1 <= e0:
                return points[i].sigma_b
            t = (energy_mev - e0) / (e1 - e0)
            return points[i - 1].sigma_b * (1 - t) + points[i].sigma_b * t
    return points[-1].sigma_b


def _interp_sigma_b(res: ResidualSigma, energy_mev: float) -> float:
    return _interp_sigma_b_points(res.points, energy_mev)


def _max_interp_mb(
    points: list[SigmaPoint], e_lo_mev: float, e_hi_mev: float
) -> tuple[float | None, float | None]:
    if not points or e_hi_mev < e_lo_mev:
        return None, None
    candidates: list[float] = []
    for p in points:
        if e_lo_mev - 1e-9 <= p.energy_mev <= e_hi_mev + 1e-9:
            candidates.append(p.energy_mev)
    # Endpoints (Emax interpolation) when inside the tabulated span.
    for edge in (e_lo_mev, e_hi_mev):
        if points[0].energy_mev - 1e-12 <= edge <= points[-1].energy_mev + 1e-12:
            candidates.append(edge)
    if not candidates:
        return None, None
    best_mb: float | None = None
    best_e: float | None = None
    for e in sorted(set(candidates)):
        s = _interp_sigma_b_points(points, e)
        if s <= 0:
            continue
        mb = s * 1e3
        if best_mb is None or mb > best_mb:
            best_mb = mb
            best_e = e
    return best_mb, best_e


def _sum_lfs_max(
    by_level: dict[int, ResidualSigma],
    e_lo_mev: float,
    e_hi_mev: float,
) -> tuple[float | None, float | None]:
    """σ_sum(E)=Σ_LFS σ_LFS(E); return (max_mb, E_at_max) in [e_lo, e_hi]."""
    if not by_level:
        return None, None
    energies: set[float] = set()
    spans: list[tuple[float, float]] = []
    for res in by_level.values():
        if not res.points:
            continue
        spans.append((res.points[0].energy_mev, res.points[-1].energy_mev))
        for p in res.points:
            if e_lo_mev - 1e-9 <= p.energy_mev <= e_hi_mev + 1e-9:
                energies.add(p.energy_mev)
    if spans:
        lo_tab = min(s[0] for s in spans)
        hi_tab = max(s[1] for s in spans)
        for edge in (e_lo_mev, e_hi_mev):
            if lo_tab - 1e-12 <= edge <= hi_tab + 1e-12:
                energies.add(edge)
    if not energies:
        return None, None
    best_mb: float | None = None
    best_e: float | None = None
    for e in sorted(energies):
        total_b = sum(_interp_sigma_b(res, e) for res in by_level.values())
        if total_b <= 0:
            continue
        mb = total_b * 1e3
        if best_mb is None or mb > best_mb:
            best_mb = mb
            best_e = e
    return best_mb, best_e


def _endf_floats(line: str) -> list[float | None]:
    out: list[float | None] = []
    body = line[:66]
    for j in range(0, 66, 11):
        field = body[j : j + 11]
        if not field.strip():
            out.append(None)
            continue
        out.append(_endf_float(field))
    return out


def _endf_float(field: str) -> float | None:
    s = field.strip().replace(" ", "")
    if not s:
        return None
    m = re.match(r"^([+-]?\d*\.\d+|[+-]?\d+)([+-]\d+)$", s)
    if m:
        return float(f"{m.group(1)}e{m.group(2)}")
    try:
        return float(s)
    except ValueError:
        return None
