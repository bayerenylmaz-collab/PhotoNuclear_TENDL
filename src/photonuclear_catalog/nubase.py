"""NUBASE2020 parser for half-lives and isomers."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from photonuclear_catalog.nuclide import Nuclide, z_to_element


@dataclass(frozen=True)
class NubaseRecord:
    z: int
    a: int
    isomer_index: int
    isomer_label: str  # "", "m", "n", ...
    half_life_raw: str
    half_life_unit: str
    is_stable: bool
    is_particle_unstable: bool
    excitation_kev: float | None
    abundance: float | None
    jp: str
    decay: str
    half_life_estimated: bool = False
    excitation_estimated: bool = False

    @property
    def nuclide(self) -> Nuclide:
        return Nuclide(z=self.z, a=self.a, isomer=self.isomer_label)

    @property
    def is_metastable(self) -> bool:
        """True for nuclear isomers (not IAS/resonance placeholders)."""
        if self.isomer_index in {5, 8, 9}:
            return False
        if self.isomer_label in {"r", "i", "j"}:
            return False
        if self.isomer_index > 0 and self.isomer_label in {"m", "n", "p", "q", "x"}:
            return True
        return self.isomer_index in {1, 2} and bool(self.isomer_label)

    @property
    def is_radioactive(self) -> bool:
        return (not self.is_stable) and (not self.is_particle_unstable)

    @property
    def include_in_catalog(self) -> bool:
        """Ground-state radioactives and true metastable isomers only."""
        if self.is_particle_unstable:
            return False
        if self.isomer_index == 0 and not self.isomer_label:
            return self.is_radioactive
        return self.is_metastable and (self.is_radioactive or not self.is_stable)

    @property
    def display_isomer(self) -> str:
        """Canonical isomer tag for reports: m1, m2, m3, ..."""
        return isomer_display_label(self.isomer_label, self.isomer_index)


class NubaseTable:
    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self.records: list[NubaseRecord] = []
        self._by_nuclide: dict[tuple[int, int, str], NubaseRecord] = {}
        self._isomers_by_za: dict[tuple[int, int], list[NubaseRecord]] = {}
        self._load()

    def _load(self) -> None:
        for line in self.path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line or line.startswith("#"):
                continue
            if len(line) < 80:
                continue
            try:
                a = int(line[0:3])
                zzzi = line[4:8]
                z = int(zzzi[0:3])
                isomer_index = int(zzzi[3])
                label_field = line[11:16]
                isomer_flag = line[16:17].strip()
                # Excitation energy columns 43:54
                exc_raw = line[42:54].strip()
                hl_raw = line[69:78].strip()
                hl_unit = line[78:80].strip()
                jp = line[88:102].strip()
                decay = line[119:].strip() if len(line) > 119 else ""
            except ValueError:
                continue

            isomer_label = _isomer_label(isomer_index, isomer_flag, label_field)
            is_stable = hl_raw.lower().startswith("stbl")
            is_p_unst = "p-unst" in hl_raw.lower()
            hl_estimated = "#" in hl_raw
            exc_estimated = "#" in exc_raw
            excitation = _parse_optional_float(exc_raw) if isomer_index > 0 else 0.0
            if isomer_index == 0:
                excitation = 0.0
                exc_estimated = False
            abundance = _parse_abundance(decay)

            rec = NubaseRecord(
                z=z,
                a=a,
                isomer_index=isomer_index,
                isomer_label=isomer_label,
                half_life_raw=hl_raw,
                half_life_unit=hl_unit,
                is_stable=is_stable,
                is_particle_unstable=is_p_unst,
                excitation_kev=excitation,
                abundance=abundance,
                jp=jp,
                decay=decay,
                half_life_estimated=hl_estimated,
                excitation_estimated=exc_estimated,
            )
            self.records.append(rec)
            self._by_nuclide[(z, a, isomer_label)] = rec
            self._isomers_by_za.setdefault((z, a), []).append(rec)

        if not self.records:
            raise RuntimeError(f"No NUBASE records parsed from {self.path}")

    def get(self, z: int, a: int, isomer: str = "") -> NubaseRecord | None:
        return self._by_nuclide.get((z, a, isomer))

    def states_for(self, z: int, a: int) -> list[NubaseRecord]:
        return list(self._isomers_by_za.get((z, a), []))

    def is_stable_ground(self, z: int, a: int) -> bool:
        gs = self.get(z, a, "")
        return bool(gs and gs.is_stable)

    def stable_with_abundance(self) -> list[NubaseRecord]:
        out: list[NubaseRecord] = []
        for rec in self.records:
            if rec.isomer_index == 0 and rec.is_stable:
                out.append(rec)
        return out


def _isomer_label(isomer_index: int, flag: str, label_field: str) -> str:
    if isomer_index == 0 and not flag:
        return ""
    flag = flag.lower()
    if flag in {"m", "n", "p", "q", "r", "i", "j", "x"}:
        # NUBASE: first isomer often "m", second "n"
        return flag
    if isomer_index == 1:
        return "m"
    if isomer_index == 2:
        return "n"
    if isomer_index > 0:
        return f"m{isomer_index}"
    # Sometimes label embeds m, e.g. "99mTc" style in A El field
    compact = label_field.strip()
    for ch in ("m", "n"):
        if ch in compact.lower():
            return ch
    return ""


def _parse_optional_float(text: str) -> float | None:
    if not text or text == "*":
        return None
    try:
        return float(text.replace("#", ""))
    except ValueError:
        return None


def _parse_abundance(decay: str) -> float | None:
    if "IS=" not in decay:
        return None
    try:
        part = decay.split("IS=", 1)[1].strip().split()[0]
        return float(part)
    except (IndexError, ValueError):
        return None


def default_nubase_path() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "nubase_4.mas20.txt"


def format_half_life(rec: NubaseRecord) -> str:
    if rec.is_stable:
        return "stable"
    if rec.is_particle_unstable:
        return "p-unstable"
    if not rec.half_life_raw:
        return ""
    text = f"{rec.half_life_raw} {rec.half_life_unit}".strip()
    if rec.half_life_estimated and "#" not in text:
        text = f"{text} #"
    return text


def isomer_display_label(isomer_label: str, isomer_index: int = 0) -> str:
    """Map NUBASE m/n/p/q (or m2) onto m1/m2/m3 display tags."""
    lab = (isomer_label or "").lower().strip()
    if not lab:
        return ""
    mapping = {"m": "m1", "n": "m2", "p": "m3", "q": "m4", "x": "mx"}
    if lab in mapping:
        return mapping[lab]
    if re.fullmatch(r"m\d+", lab):
        return lab if lab != "m1" else "m1"
    if isomer_index > 0:
        return f"m{isomer_index}"
    return lab


def pretty_name(
    z: int,
    a: int,
    isomer: str = "",
    *,
    isomer_index: int = 0,
    display: bool = True,
) -> str:
    try:
        el = z_to_element(z)
    except ValueError:
        el = f"Z{z}"
    tag = isomer_display_label(isomer, isomer_index) if display else (isomer or "")
    return f"{a}{el}{tag}"


def pretty_name_nuclear(
    z: int,
    a: int,
    isomer: str = "",
    *,
    isomer_index: int = 0,
) -> str:
    """Human-report form like ¹⁹⁶ᵐ¹Au (ASCII fallback: ^{196m1}Au)."""
    try:
        el = z_to_element(z)
    except ValueError:
        el = f"Z{z}"
    tag = isomer_display_label(isomer, isomer_index)
    # Prefer compact ASCII nuclear notation requested in reviews.
    if tag:
        return f"^{{{a}{tag}}}{el}"
    return f"^{{{a}}}{el}"


def pretty_name_html(
    z: int,
    a: int,
    isomer: str = "",
    *,
    isomer_index: int = 0,
) -> str:
    """HTML superscript form: <sup>196m1</sup>Au."""
    try:
        el = z_to_element(z)
    except ValueError:
        el = f"Z{z}"
    tag = isomer_display_label(isomer, isomer_index)
    return f"<sup>{a}{tag}</sup>{el}"
