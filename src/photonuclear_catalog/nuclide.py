"""Nuclide parsing and element symbol helpers."""

from __future__ import annotations

import re
from dataclasses import dataclass

_ELEMENT_TO_Z: dict[str, int] = {
    "n": 0,
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Nh": 113,
    "Fl": 114,
    "Mc": 115,
    "Lv": 116,
    "Ts": 117,
    "Og": 118,
}

_Z_TO_ELEMENT = {z: el for el, z in _ELEMENT_TO_Z.items() if el != "n"}

# Supports: 208Pb, Pb-208, 99mTc, 207Pbm, Tc99m, 180mTa, 192Irm1
_NUCLIDE_RE = re.compile(
    r"""
    ^\s*
    (?:
        (?P<a1>\d{1,3})\s*[-_]?\s*(?P<body1>[A-Za-z]{1,5})
      | (?P<body2>[A-Za-z]{1,3})\s*[-_]?\s*(?P<a2>\d{1,3})\s*(?P<meta2>m\d*|n|p|q)?
    )
    \s*$
    """,
    re.VERBOSE,
)


@dataclass(frozen=True, order=True)
class Nuclide:
    z: int
    a: int
    isomer: str = ""  # "", "m", "n", "m2", ...

    @property
    def n(self) -> int:
        return self.a - self.z

    @property
    def symbol(self) -> str:
        if self.z == 0:
            el = "n"
        else:
            el = _Z_TO_ELEMENT.get(self.z, f"Z{self.z}")
        return f"{self.a}{el}{self.isomer}"

    @property
    def livechart_id(self) -> str:
        """IAEA LiveChart nuclide id, e.g. 60co, 99mtc."""
        if self.z == 0:
            return "1n"
        el = _Z_TO_ELEMENT[self.z].lower()
        meta = self.isomer.lower()
        return f"{self.a}{meta}{el}"


def element_to_z(symbol: str) -> int:
    key = symbol.strip()
    # Only bare lowercase "n" is the neutron; "N" is nitrogen.
    if key == "n":
        return 0
    if key in _ELEMENT_TO_Z:
        return _ELEMENT_TO_Z[key]
    titled = key[:1].upper() + key[1:].lower()
    if titled in _ELEMENT_TO_Z:
        return _ELEMENT_TO_Z[titled]
    raise ValueError(f"Unknown element symbol: {symbol!r}")


def z_to_element(z: int) -> str:
    if z == 0:
        return "n"
    if z not in _Z_TO_ELEMENT:
        raise ValueError(f"Unknown Z: {z}")
    return _Z_TO_ELEMENT[z]


def _split_element_and_isomer(body: str) -> tuple[str, str]:
    """Split 'Pbm', 'mTc', 'Pb', 'Tc', 'Pbm1' into (element, isomer).

    Prefer an exact element match first so two-letter symbols that start with an
    isomer letter (Ni, Nb, …) are not misread as prefix-isomer + element
    (n+I, n+B, …).
    """
    body = body.strip()
    # Bare element first (Ni, Na, N, Pb, …).
    try:
        element_to_z(body)
        return body[:1].upper() + body[1:].lower() if body != "n" else "n", ""
    except ValueError:
        pass
    # Prefix isomer style: mTc, nAg, m1Tc
    for prefix in ("m1", "m2", "m3", "m4", "m", "n"):
        if body.lower().startswith(prefix) and len(body) > len(prefix):
            rest = body[len(prefix) :]
            try:
                element_to_z(rest)
                el = rest[:1].upper() + rest[1:].lower()
                return el, _normalize_isomer(prefix)
            except ValueError:
                pass
    # Suffix isomer style: Pbm, Pbn, Pbm1, Tlq
    for suffix in ("m1", "m2", "m3", "m4", "m", "n", "p", "q"):
        if body.lower().endswith(suffix) and len(body) > len(suffix):
            el = body[: -len(suffix)]
            try:
                element_to_z(el)
                titled = el[:1].upper() + el[1:].lower()
                return titled, _normalize_isomer(suffix)
            except ValueError:
                continue
    raise ValueError(f"Cannot split element/isomer from {body!r}")


def _normalize_isomer(tag: str) -> str:
    """Store m1 as m (NUBASE first isomer), keep m2+ and n/p/q."""
    t = tag.lower()
    if t == "m1":
        return "m"
    if t == "m2":
        return "n"
    if t == "m3":
        return "p"
    if t == "m4":
        return "q"
    return t


def parse_nuclide(text: str) -> Nuclide:
    """Parse forms like 208Pb, Pb-208, 99mTc, 207Pbm, Tc99m, 192Irm1."""
    raw = text.strip()
    m = _NUCLIDE_RE.match(raw)
    if not m:
        raise ValueError(
            f"Cannot parse nuclide {text!r}. Use forms like 208Pb, Pb-208, 99mTc."
        )
    if m.group("a1"):
        a = int(m.group("a1"))
        el, isomer = _split_element_and_isomer(m.group("body1"))
    else:
        a = int(m.group("a2"))
        el = m.group("body2")
        isomer = _normalize_isomer(m.group("meta2") or "")
        element_to_z(el)
    z = element_to_z(el)
    isomer = isomer.lower()
    if isomer == "m1":
        isomer = "m"
    return Nuclide(z=z, a=a, isomer=isomer)
