"""Command-line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from photonuclear_catalog import __version__
from photonuclear_catalog.catalog import build_catalog, export_result, load_tables
from photonuclear_catalog.gammas import GammaClient
from photonuclear_catalog.nubase import pretty_name
from photonuclear_catalog.nuclide import parse_nuclide


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="photonuclear-catalog",
        description=(
            "Kararlı bir hedeften (γ,xn yp) kinematik radyoaktif adayları, "
            "Coulomb/TENDL notları ve en güçlü mutlak bozunma gamalarını listeler."
        ),
    )
    p.add_argument("--target", help="Hedef çekirdek, örn. 208Pb veya Pb-208")
    p.add_argument(
        "--emax",
        type=float,
        default=None,
        help="Maksimum foton / eşik enerjisi (MeV), varsayılan 45, üst sınır 45",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Çıktı klasörü (varsayılan: out/<hedef>_E<max>)",
    )
    p.add_argument(
        "--gamma-limit",
        type=int,
        default=10,
        help="Ürün başına en fazla gama sayısı (varsayılan: 10)",
    )
    p.add_argument(
        "--min-gamma-kev",
        type=float,
        default=0.0,
        help=(
            "İsteğe bağlı nükleer γ alt enerji kesimi keV (varsayılan: 0). "
            "X-ışınları enerji kesimiyle değil karakteristik X-ray kimliğiyle elenir."
        ),
    )
    p.add_argument(
        "--keep-xrays",
        action="store_true",
        help="Atomik X-ışını filtresini kapat (debug)",
    )
    p.add_argument(
        "--no-tendl",
        action="store_true",
        help="TENDL-2023 artıksal σ indirmesini/atamasını atla",
    )
    p.add_argument(
        "--tendl-cache",
        type=Path,
        default=Path("data/tendl_cache"),
        help="TENDL ENDF önbellek klasörü",
    )
    p.add_argument(
        "--no-gammas",
        action="store_true",
        help="Gama çekimini atla; yalnız ürün listesi yaz",
    )
    p.add_argument(
        "--interactive",
        "-i",
        action="store_true",
        help="Hedef ve enerjiyi sorarak çalış (e-posta ile paylaşım için kolay yol)",
    )
    p.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(".gamma_cache"),
        help="Gama API önbellek klasörü",
    )
    p.add_argument("--mass-file", type=Path, default=None)
    p.add_argument("--nubase-file", type=Path, default=None)
    p.add_argument(
        "--list-stable",
        action="store_true",
        help="NUBASE kararlı temel durumları yazdır ve çık",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def _prompt(text: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{text}{suffix}: ").strip()
        if raw:
            return raw
        if default is not None:
            return default
        print("  Lütfen bir değer girin.")


def run_interactive(args: argparse.Namespace) -> argparse.Namespace:
    print()
    print("=== Photonuclear Katalog (etkileşimli) ===")
    print(
        "Kararlı hedef + enerji limiti verin; kinematik adaylar, "
        "Coulomb/TENDL notları ve mutlak Iγ üretilir."
    )
    print()
    target = _prompt("Hedef çekirdek (örn. 208Pb)", "208Pb")
    emax_s = _prompt("Enerji limiti E_max MeV (max 45)", "45")
    try:
        emax = float(emax_s.replace(",", "."))
    except ValueError:
        print("Geçersiz enerji; 45 MeV kullanılacak.")
        emax = 45.0
    tendl = _prompt("TENDL-2023 σ kullanılsın mı? (e/h)", "e").lower()
    args.target = target
    args.emax = emax
    args.no_tendl = not (tendl.startswith("e") or tendl.startswith("y"))
    if args.out is None:
        safe = target.replace("-", "").replace(" ", "")
        args.out = Path("out") / f"{safe}_E{emax:g}"
    return args


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    # No args => interactive (mail-friendly double-click / run.bat flow)
    if not argv:
        argv = ["--interactive"]
    args = parser.parse_args(argv)

    if args.list_stable:
        _, nubase = load_tables(args.mass_file, args.nubase_file)
        for rec in nubase.stable_with_abundance():
            abund = f" IS={rec.abundance}" if rec.abundance is not None else ""
            print(f"{pretty_name(rec.z, rec.a)}{abund}")
        return 0

    if args.interactive or not args.target:
        try:
            args = run_interactive(args)
        except EOFError:
            print("error: etkileşimli giriş iptal", file=sys.stderr)
            return 2

    if args.emax is None:
        args.emax = 45.0
    if args.emax > 45.0:
        print("error: --emax en fazla 45 MeV olabilir", file=sys.stderr)
        return 2
    if args.emax <= 0:
        print("error: --emax pozitif olmalı", file=sys.stderr)
        return 2

    try:
        target = parse_nuclide(args.target)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.out is None:
        args.out = Path("out") / f"{target.symbol}_E{args.emax:g}"

    masses, nubase = load_tables(args.mass_file, args.nubase_file)
    gs = nubase.get(target.z, target.a, "")
    if gs is None:
        print(f"error: hedef {target.symbol} NUBASE2020'de yok", file=sys.stderr)
        return 2
    if not gs.is_stable:
        print(
            f"uyarı: {target.symbol} NUBASE'de kararlı değil; yine de devam",
            file=sys.stderr,
        )

    print()
    print(f"Hedef : {target.symbol}")
    print(f"E_max : {args.emax:g} MeV")
    print(
        f"Gama  : en fazla {args.gamma_limit} çizgi; ENSDF Iγ "
        f"(branch_relative ise ×BR); X-ray=LiveChart rad_types=x"
        + (
            f"; Eγ ≥ {args.min_gamma_kev:g} keV"
            if args.min_gamma_kev > 0
            else ""
        )
    )
    print(f"TENDL : {'kapalı' if args.no_tendl else 'açık (MF10 LFS residual σ)'}")
    print("Çalışıyor...")

    client = None if args.no_gammas else GammaClient(cache_dir=args.cache_dir)
    result = build_catalog(
        target_text=args.target,
        emax_mev=args.emax,
        gamma_limit=args.gamma_limit,
        fetch_gammas=not args.no_gammas,
        cache_dir=args.cache_dir,
        masses=masses,
        nubase=nubase,
        gamma_client=client,
        min_gamma_kev=args.min_gamma_kev,
        use_tendl=not args.no_tendl,
        tendl_cache_dir=args.tendl_cache,
        exclude_xrays=not args.keep_xrays,
    )
    paths = export_result(result, args.out)

    open_n = sum(1 for p in result.products if p.viability == "kinematic_open")
    stable_n = sum(1 for p in result.products if p.is_stable)
    print()
    print("Bitti.")
    print(
        f"  aday ürün  : {len(result.products)} "
        f"(open: {open_n}, stable: {stable_n})"
    )
    print(f"  gama satırı : {len(result.gammas)}")
    print(f"  klasör      : {args.out.resolve()}")
    print(f"  özet (HTML) : {paths['report_html'].resolve()}")
    print(f"  özet (MD)   : {paths['report_md'].resolve()}")
    print(f"  products.csv: {paths['products'].resolve()}")
    print(f"  gammas.csv  : {paths['gammas'].resolve()}")
    print()
    print("İpucu: report.html dosyasını çift tıklayıp tarayıcıda açabilirsiniz.")
    print("Not: Eth kinematiktir; TENDL σ, Coulomb heuristic’ten önceliklidir.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
