"""CSV / JSON exporters."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable, Sequence


def write_csv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        fields = list(fieldnames or [])
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
        return
    fields = list(fieldnames or rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def rows_from_objects(objects: Iterable[object]) -> list[dict]:
    out: list[dict] = []
    for obj in objects:
        if hasattr(obj, "to_dict"):
            out.append(obj.to_dict())  # type: ignore[no-untyped-call]
        elif isinstance(obj, dict):
            out.append(obj)
        else:
            raise TypeError(f"Cannot serialize {type(obj)!r}")
    return out
