"""AME2020 atomic mass table parser."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Neutrons and hydrogen mass excesses (keV) from AME2020.
ME_NEUTRON_KEV = 8071.31806
ME_HYDROGEN_KEV = 7288.971064


@dataclass(frozen=True)
class MassRecord:
    z: int
    a: int
    n: int
    element: str
    mass_excess_kev: float
    estimated: bool


class MassTable:
    """Ground-state atomic mass excesses from AME2020 mass_1.mas20."""

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self._by_za: dict[tuple[int, int], MassRecord] = {}
        self._load()

    def _load(self) -> None:
        lines = self.path.read_text(encoding="utf-8", errors="replace").splitlines()
        for line in lines:
            # Fixed-width AME line with Fortran carriage control in column 1.
            if len(line) < 42:
                continue
            try:
                n = int(line[4:9])
                z = int(line[9:14])
                a = int(line[14:19])
            except ValueError:
                continue
            el = line[20:23].strip()
            if not el:
                continue
            me_field = line[28:42].strip()
            if not me_field or me_field == "*":
                continue
            estimated = "#" in me_field
            try:
                me = float(me_field.replace("#", "."))
            except ValueError:
                continue
            if a < 1 or z < 0:
                continue
            self._by_za[(z, a)] = MassRecord(
                z=z,
                a=a,
                n=n,
                element=el,
                mass_excess_kev=me,
                estimated=estimated,
            )

        if not self._by_za:
            raise RuntimeError(f"No mass records parsed from {self.path}")

    def get(self, z: int, a: int) -> MassRecord | None:
        return self._by_za.get((z, a))

    def has(self, z: int, a: int) -> bool:
        return (z, a) in self._by_za

    def __len__(self) -> int:
        return len(self._by_za)


def default_mass_path() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "mass_1.mas20.txt"
