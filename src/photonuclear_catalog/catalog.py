"""High-level catalog build orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from photonuclear_catalog.gammas import GammaClient, GammaLine
from photonuclear_catalog.io_export import rows_from_objects, write_csv, write_json
from photonuclear_catalog.masses import MassTable, default_mass_path
from photonuclear_catalog.nubase import NubaseTable, default_nubase_path
from photonuclear_catalog.nuclide import Nuclide, parse_nuclide
from photonuclear_catalog.reactions import (
    DEFAULT_SIGMA_NEGLIGIBLE_MB,
    Product,
    enumerate_products,
    filter_by_emax,
)
from photonuclear_catalog.report import write_reports
from photonuclear_catalog.tendl import TendlPhotonuclear


@dataclass
class CatalogResult:
    target: str
    emax_mev: float
    products: list[Product]
    gammas: list[GammaLine]
    min_gamma_kev: float = 0.0
    gamma_limit: int = 10
    use_tendl: bool = True
    meta: dict = field(default_factory=dict)


def load_tables(
    mass_path: Path | None = None,
    nubase_path: Path | None = None,
) -> tuple[MassTable, NubaseTable]:
    masses = MassTable(mass_path or default_mass_path())
    nubase = NubaseTable(nubase_path or default_nubase_path())
    return masses, nubase


def build_catalog(
    target_text: str,
    emax_mev: float = 45.0,
    gamma_limit: int = 10,
    fetch_gammas: bool = True,
    cache_dir: Path | None = None,
    masses: MassTable | None = None,
    nubase: NubaseTable | None = None,
    gamma_client: GammaClient | None = None,
    min_gamma_kev: float = 0.0,
    use_tendl: bool = True,
    tendl_cache_dir: Path | None = None,
    exclude_xrays: bool = True,
    include_stable_products: bool = True,
    sigma_negligible_mb: float = DEFAULT_SIGMA_NEGLIGIBLE_MB,
) -> CatalogResult:
    if emax_mev <= 0:
        raise ValueError("emax must be positive")
    if emax_mev > 45.0:
        raise ValueError("emax must be <= 45 MeV for this task")

    target = parse_nuclide(target_text)
    if masses is None or nubase is None:
        masses, nubase = load_tables()

    gs = nubase.get(target.z, target.a, "")
    if gs is None:
        raise ValueError(f"Target {target.symbol} not found in NUBASE2020.")

    tendl = None
    if use_tendl:
        tendl = TendlPhotonuclear(cache_dir=tendl_cache_dir or Path("data/tendl_cache"))

    products = enumerate_products(
        target,
        masses,
        nubase,
        emax_mev=emax_mev,
        include_stable_products=include_stable_products,
        tendl=tendl,
        sigma_negligible_mb=sigma_negligible_mb,
    )
    gammas: list[GammaLine] = []
    if fetch_gammas and products:
        client = gamma_client or GammaClient(cache_dir=cache_dir)
        unique: dict[tuple[int, int, str], Product] = {}
        for p in products:
            if p.is_stable:
                continue  # stable residuals: no decay gammas
            unique.setdefault((p.z, p.a, p.isomer), p)
        for p in unique.values():
            nuc = Nuclide(z=p.z, a=p.a, isomer=p.isomer)
            try:
                lines = client.top_gammas(
                    nuc,
                    limit=gamma_limit,
                    excitation_kev=p.excitation_kev,
                    min_energy_kev=min_gamma_kev,
                    expected_half_life=p.half_life,
                    exclude_xrays=exclude_xrays,
                )
                for g in lines:
                    gammas.append(
                        GammaLine(
                            product=p.product,
                            energy_kev=g.energy_kev,
                            intensity=g.intensity,
                            intensity_unc=g.intensity_unc,
                            decay_mode=g.decay_mode,
                            parent_energy_kev=g.parent_energy_kev,
                            half_life=g.half_life,
                            rank=g.rank,
                            intensity_ensdf=g.intensity_ensdf,
                            intensity_relative=g.intensity_relative,
                            branch_percent=g.branch_percent,
                            normalization_status=g.normalization_status,
                            uncertainty_status=g.uncertainty_status,
                            status=g.status,
                            note=g.note,
                            parent_match_conflict=g.parent_match_conflict,
                        )
                    )
            except Exception as exc:  # noqa: BLE001
                gammas.append(
                    GammaLine(
                        product=p.product,
                        energy_kev=float("nan"),
                        intensity=None,
                        intensity_unc=None,
                        decay_mode="NETWORK_ERROR",
                        parent_energy_kev=None,
                        half_life="",
                        rank=None,
                        status="network_error",
                        uncertainty_status="status_row",
                        note=str(exc),
                    )
                )

    return CatalogResult(
        target=target.symbol,
        emax_mev=emax_mev,
        products=products,
        gammas=gammas,
        min_gamma_kev=min_gamma_kev,
        gamma_limit=gamma_limit,
        use_tendl=use_tendl,
        meta={
            "exclude_xrays": exclude_xrays,
            "include_stable_products": include_stable_products,
            "sigma_negligible_mb": sigma_negligible_mb,
            "intensity_policy": (
                "Keep ENSDF/LiveChart Iγ as absolute unless a decay mode is "
                "detected as branch_relative (any Iγ > BR%); only then I_abs=I×BR/100. "
                "If decay_rads has no parent for an isomer, fall back to LiveChart "
                "fields=gammas (Adopted Levels) with relative_adopted_level (not absolute %)."
            ),
            "xray_policy": "Subtract LiveChart rad_types=x keys from rad_types=g.",
            "tendl_policy": (
                "LFS=0 is ground-state production; residual total = Σ_LFS σ(E). "
                "Isomers without LFS match → state_specific_missing (no LFS=0 copy)."
            ),
            "sources": {
                "reference_ui": "https://www.nndc.bnl.gov/nudat3/",
                "masses": "AME2020",
                "half_lives_isomers": "NUBASE2020",
                "decay_gammas": "IAEA LiveChart / ENSDF (NuDat-compatible)",
                "cross_sections": "TENDL-2023 MF10 residual production (LFS-aware)",
            },
        },
    )


def filter_catalog(result: CatalogResult, emax_mev: float) -> CatalogResult:
    products = filter_by_emax(result.products, emax_mev)
    allowed = {(p.z, p.a, p.isomer) for p in products}
    gammas = [g for g in result.gammas if _gamma_key(g.product) in allowed]
    return CatalogResult(
        target=result.target,
        emax_mev=emax_mev,
        products=products,
        gammas=gammas,
        min_gamma_kev=result.min_gamma_kev,
        gamma_limit=result.gamma_limit,
        use_tendl=result.use_tendl,
        meta=result.meta,
    )


def _gamma_key(product_symbol: str) -> tuple[int, int, str]:
    n = parse_nuclide(product_symbol)
    return n.z, n.a, n.isomer


def export_result(result: CatalogResult, out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    product_rows = rows_from_objects(result.products)
    gamma_rows = rows_from_objects(result.gammas)
    paths = {
        "products": out_dir / "products.csv",
        "gammas": out_dir / "gammas.csv",
        "catalog": out_dir / "catalog.json",
    }
    write_csv(paths["products"], product_rows)
    write_csv(paths["gammas"], gamma_rows)
    write_json(
        paths["catalog"],
        {
            "target": result.target,
            "emax_mev": result.emax_mev,
            "n_products": len(result.products),
            "n_gamma_rows": len(result.gammas),
            "min_gamma_kev": result.min_gamma_kev,
            "gamma_limit": result.gamma_limit,
            "use_tendl": result.use_tendl,
            "products": product_rows,
            "gammas": gamma_rows,
            "meta": result.meta,
            "sources": result.meta.get("sources", {}),
            "notes": [
                "Eth is the kinematic Q threshold.",
                "Barrier indicator = Eth + Vc (approximate; not a sharp second threshold).",
                "When available, TENDL residual cross sections take precedence over the Coulomb estimate.",
                "TENDL matching uses MT + IZAP + LFS (free-nucleon channels; no clustered fallback).",
                "LFS=0 is ground-state production; MF3 is the channel total; Σ_LFS is the MF10 level sum.",
                "Default TENDL source is s60, with URL and sha256 provenance.",
                "Iγ values are absolute ENSDF intensities unless a branch-relative dataset is detected.",
                "For isomers missing from decay radiation tables, Adopted Levels relative intensities may be used.",
                "Atomic X-rays are excluded from the nuclear γ list.",
                "Channel labels are free-nucleon balances; clustered channels such as (γ,α) are out of scope.",
            ],
        },
    )
    md_path, html_path = write_reports(
        result,
        out_dir,
        min_gamma_kev=result.min_gamma_kev,
        gamma_limit=result.gamma_limit,
    )
    paths["report_md"] = md_path
    paths["report_html"] = html_path
    return paths
