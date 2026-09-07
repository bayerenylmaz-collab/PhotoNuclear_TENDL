"""Human-readable Markdown / HTML summary reports."""

from __future__ import annotations

import html
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

from photonuclear_catalog.nubase import pretty_name_html
from photonuclear_catalog.reactions import DEFAULT_SIGMA_NEGLIGIBLE_MB

if TYPE_CHECKING:
    from photonuclear_catalog.catalog import CatalogResult
    from photonuclear_catalog.gammas import GammaLine
    from photonuclear_catalog.reactions import Product


def write_reports(
    result: "CatalogResult",
    out_dir: Path,
    *,
    min_gamma_kev: float,
    gamma_limit: int,
) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / "report.md"
    html_path = out_dir / "report.html"
    md_path.write_text(
        render_markdown(result, min_gamma_kev=min_gamma_kev, gamma_limit=gamma_limit),
        encoding="utf-8",
    )
    html_path.write_text(
        render_html(result, min_gamma_kev=min_gamma_kev, gamma_limit=gamma_limit),
        encoding="utf-8",
    )
    return md_path, html_path


def render_markdown(
    result: "CatalogResult",
    *,
    min_gamma_kev: float,
    gamma_limit: int,
) -> str:
    by_gamma = _group_gammas(result.gammas)
    lines: list[str] = []
    lines.append(f"# Photonuclear kinematik aday kataloğu — {result.target}")
    lines.append("")
    lines.append(f"- **Hedef çekirdek:** `{result.target}`")
    lines.append(f"- **Enerji limiti (E_max):** {result.emax_mev:g} MeV")
    lines.append(f"- **Aday ürün sayısı:** {len(result.products)}")
    negl = float(
        (result.meta or {}).get("sigma_negligible_mb", DEFAULT_SIGMA_NEGLIGIBLE_MB)
    )
    lines.append(
        f"- **σ_negligible eşiği:** {negl:g} mb "
        f"(TENDL σ_max bu eşiğin altındaysa ürün `sigma_negligible` olarak işaretlenir)"
    )
    lines.append(
        f"- **Gama seçimi:** en fazla {gamma_limit} çizgi; "
        f"öncelik mutlak Iγ; X-ışınları elenir"
        + (
            f"; ek kesim Eγ ≥ {min_gamma_kev:g} keV"
            if min_gamma_kev > 0
            else ""
        )
    )
    lines.append(
        f"- **Oluşturma:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}"
    )
    lines.append("")
    lines.append("## Okuma notu")
    lines.append("")
    lines.append(
        "- Eth: kinematik (Q) eşiği.  \n"
        "- Vc: Coulomb bariyer ölçeği; Eth+Vc yalnızca kaba bariyer göstergesidir "
        "(keskin ikinci eşik değildir; tünelleme olabilir).  \n"
        "- TENDL artıksal kesit varsa karar için önceliklidir.  \n"
        "- TENDL eşlemesi MT + IZAP + LFS ile yapılır; serbest nükleon kanalları "
        "küme kanallarıyla `(γ,d)/(γ,α)` yedeklenmez.  \n"
        "- LFS=0 temel durum üretimini; MF3 kanal toplamını; Σ_LFS ise MF10 seviye "
        "toplamını temsil eder.  \n"
        "- İzomer için seviye-özgü LFS yoksa yalnızca kinematik bilgi verilir "
        "(temel durum kesiti kopyalanmaz).  \n"
        "- Kanal etiketi serbest nükleon dengesidir; aynı (Z,A) başka kanallarla da "
        "oluşabilir."
    )
    lines.append("")
    lines.append("## Veri kaynakları")
    lines.append("")
    lines.append(
        "- Referans arayüz: [NuDat 3](https://www.nndc.bnl.gov/nudat3/)  \n"
        "- Eşik enerjileri: AME2020  \n"
        "- Yarı ömür / izomer: NUBASE2020  \n"
        "- Bozunma gamaları: ENSDF (IAEA LiveChart)  \n"
        "- Artıksal σ: TENDL-2023 MF3 (kanal) + MF10 (MT+IZAP+LFS)"
    )
    lines.append("")
    lines.append("## Ürünler ve en güçlü bozunma gamaları")
    lines.append("")

    if not result.products:
        lines.append("_Bu enerji limitinde kinematik aday bulunamadı._")
        lines.append("")
        return "\n".join(lines)

    for idx, product in enumerate(result.products, start=1):
        gammas = by_gamma.get(product.product, [])
        title = product.product_pretty or product.product
        lines.append(f"### {idx}. `{title}` (`{product.product}`)")
        lines.append("")
        lines.append(
            f"- Kanal: `{product.channel}` — {product.channel_note}  \n"
            f"- Eth (kinematik): **{product.eth_mev:.3f} MeV**  \n"
            f"- Vc (Coulomb ölçeği): {product.coulomb_barrier_mev:.2f} MeV  \n"
            f"- Bariyer göstergesi (Eth+Vc): "
            f"**{product.heuristic_barrier_indicator_mev:.3f} MeV**  \n"
            f"- σ_max / E(σ_max): {_fmt_sigma(product)}  \n"
            f"- σ kapsamı: `{product.sigma_scope}`  \n"
            f"- TENDL MT: `{product.sigma_mt if product.sigma_mt is not None else '—'}`  \n"
            f"- σ_gs (LFS=0, mb): "
            f"{_fmt_optional_num(product.sigma_ground_state_reference_mb)}  \n"
            f"- σ_kanal_toplam = MF3(MT) (mb): "
            f"{_fmt_optional_num(product.sigma_channel_total_mb)}"
            f"{'' if product.e_at_channel_total_max_mev is None else f' @ {product.e_at_channel_total_max_mev:g} MeV'}  \n"
            f"- σ_ΣLFS = max Σ_LFS (mb): "
            f"{_fmt_optional_num(product.sigma_sum_mf10_lfs_mb)}"
            f"{'' if product.e_at_total_max_mev is None else f' @ {product.e_at_total_max_mev:g} MeV'}  \n"
            f"- MF10/MF3 oranı: "
            f"{_fmt_optional_num(product.mf10_coverage_ratio)}  \n"
            f"- TENDL kaynağı: `{product.sigma_source_variant or '—'}`  \n"
            f"- Uygunluk: `{product.viability}`  \n"
            f"- Yarı ömür: {product.half_life or '—'}  \n"
            f"- Durum: {_state_label(product)}  \n"
            f"- Jπ: {product.jp or '—'}  \n"
            f"- Bozunma: {product.decay or '—'}  \n"
            f"- Tahmin bayrakları: {_estimate_flags(product)}"
        )
        lines.append("")

        if product.is_stable:
            lines.append("> Kararlı artıksal çekirdek — bozunma gamasi yok.")
            lines.append("")
            continue

        real_gammas = _real_gamma_lines(gammas)
        if not real_gammas:
            lines.append(f"> **{_gamma_status_message(gammas)}**")
            lines.append("")
            continue

        conflict = any(g.parent_match_conflict for g in real_gammas)
        if conflict:
            lines.append(
                "> **Uyarı:** NUBASE yarıömrü ile ENSDF üst seviye yarıömrü uyuşmuyor; "
                "eşleşme uyarma enerjisine göre korunmuştur."
            )
            lines.append("")
        if any(g.normalization_status == "relative_adopted_level" for g in real_gammas):
            lines.append(
                "> **Not:** Bu çizgiler Adopted Levels göreli şiddetleridir "
                "(mutlak Iγ% değildir; iç dönüşüm düzeltmesi gerekebilir)."
            )
            lines.append("")

        lines.append(
            "| Sıra | Eγ (keV) | Iγ | ENSDF Iγ | Normlama | Belirsizlik | Bozunma | Not |"
        )
        lines.append("| ---: | ---: | ---: | ---: | --- | --- | --- | --- |")
        for g in real_gammas:
            lines.append(
                f"| {_fmt_rank(g.rank)} | {_fmt_energy(g.energy_kev)} | "
                f"{_fmt_intensity(g.intensity, g.uncertainty_status)} | "
                f"{_fmt_optional(g.intensity_ensdf)} | "
                f"{_fmt_norm(g.normalization_status)} | "
                f"{_fmt_unc_status(g.uncertainty_status)} | "
                f"{g.decay_mode or '—'} | {_gamma_note(g)} |"
            )
        if len(real_gammas) < gamma_limit:
            lines.append("")
            lines.append(
                f"_Not: Bu ürün için {gamma_limit} yerine {len(real_gammas)} çizgi listelendi._"
            )
        lines.append("")

    lines.append("## Dosyalar")
    lines.append("")
    lines.append(
        "- `products.csv` — ürün eşikleri, bariyer göstergesi, σ özeti  \n"
        "- `gammas.csv` — seçilen γ çizgileri ve normlama bilgisi  \n"
        "- `catalog.json` — makine-okur tam çıktı  \n"
        "- `report.html` — tarayıcı özeti"
    )
    lines.append("")
    return "\n".join(lines)


def render_html(
    result: "CatalogResult",
    *,
    min_gamma_kev: float,
    gamma_limit: int,
) -> str:
    by_gamma = _group_gammas(result.gammas)
    esc = html.escape
    negl = float(
        (result.meta or {}).get("sigma_negligible_mb", DEFAULT_SIGMA_NEGLIGIBLE_MB)
    )
    blocks: list[str] = []
    for idx, product in enumerate(result.products, start=1):
        gammas = by_gamma.get(product.product, [])
        real_gammas = _real_gamma_lines(gammas)
        title_html = pretty_name_html(
            product.z, product.a, product.isomer
        )
        head = (
            f"<section class='card'>"
            f"<h3>{idx}. {title_html} <code>{esc(product.product)}</code></h3>"
            f"<ul>"
            f"<li>Kanal: <code>{esc(product.channel)}</code> — {esc(product.channel_note)}</li>"
            f"<li>Eth (kinematik): <strong>{product.eth_mev:.3f} MeV</strong></li>"
            f"<li>Vc (Coulomb ölçeği): {product.coulomb_barrier_mev:.2f} MeV</li>"
            f"<li>Bariyer göstergesi (Eth+Vc): "
            f"<strong>{product.heuristic_barrier_indicator_mev:.3f} MeV</strong></li>"
            f"<li>σ_max / E(σ_max): {esc(_fmt_sigma(product))}</li>"
            f"<li>σ kapsamı: <code>{esc(product.sigma_scope)}</code></li>"
            f"<li>TENDL MT: <code>{esc(str(product.sigma_mt) if product.sigma_mt is not None else '—')}</code></li>"
            f"<li>σ_gs (LFS=0): "
            f"{esc(_fmt_optional_num(product.sigma_ground_state_reference_mb))}</li>"
            f"<li>σ_kanal_toplam (MF3): "
            f"{esc(_fmt_optional_num(product.sigma_channel_total_mb))}"
            f"{'' if product.e_at_channel_total_max_mev is None else esc(f' @ {product.e_at_channel_total_max_mev:g} MeV')}</li>"
            f"<li>σ_ΣLFS: "
            f"{esc(_fmt_optional_num(product.sigma_sum_mf10_lfs_mb))}"
            f"{'' if product.e_at_total_max_mev is None else esc(f' @ {product.e_at_total_max_mev:g} MeV')}</li>"
            f"<li>MF10/MF3 oranı: "
            f"{esc(_fmt_optional_num(product.mf10_coverage_ratio))}</li>"
            f"<li>TENDL kaynağı: <code>{esc(product.sigma_source_variant or '—')}</code></li>"
            f"<li>Uygunluk: <code>{esc(product.viability)}</code></li>"
            f"<li>Yarı ömür: {esc(product.half_life or '—')}</li>"
            f"<li>Durum: {esc(_state_label(product))}</li>"
            f"<li>Jπ: {esc(product.jp or '—')}</li>"
            f"<li>Bozunma: {esc(product.decay or '—')}</li>"
            f"<li>Tahmin: {esc(_estimate_flags(product))}</li>"
            f"</ul>"
        )
        if product.is_stable:
            blocks.append(
                head
                + "<p class='muted'>Kararlı artıksal çekirdek — bozunma gamasi yok.</p>"
                "</section>"
            )
            continue
        if not real_gammas:
            blocks.append(
                head
                + f"<p class='warn'><strong>{esc(_gamma_status_message(gammas))}</strong></p>"
                "</section>"
            )
            continue
        conflict = any(g.parent_match_conflict for g in real_gammas)
        conflict_html = ""
        if conflict:
            conflict_html = (
                "<p class='warn'><strong>Uyarı:</strong> "
                "NUBASE yarıömrü ile ENSDF üst seviye yarıömrü uyuşmuyor; "
                "eşleşme uyarma enerjisine göre korunmuştur.</p>"
            )
        if any(g.normalization_status == "relative_adopted_level" for g in real_gammas):
            conflict_html += (
                "<p class='muted'><strong>Not:</strong> Bu çizgiler Adopted Levels "
                "göreli şiddetleridir (mutlak Iγ% değildir; iç dönüşüm düzeltmesi "
                "gerekebilir).</p>"
            )
        rows = []
        for g in real_gammas:
            rows.append(
                "<tr>"
                f"<td>{esc(_fmt_rank(g.rank))}</td>"
                f"<td>{esc(_fmt_energy(g.energy_kev))}</td>"
                f"<td>{esc(_fmt_intensity(g.intensity, g.uncertainty_status))}</td>"
                f"<td>{esc(_fmt_optional(g.intensity_ensdf))}</td>"
                f"<td>{esc(_fmt_norm(g.normalization_status))}</td>"
                f"<td>{esc(_fmt_unc_status(g.uncertainty_status))}</td>"
                f"<td>{esc(g.decay_mode or '—')}</td>"
                f"<td>{esc(_gamma_note(g))}</td>"
                "</tr>"
            )
        note = ""
        if len(real_gammas) < gamma_limit:
            note = (
                f"<p class='muted'>Bu ürün için {gamma_limit} yerine "
                f"{len(real_gammas)} çizgi listelendi.</p>"
            )
        blocks.append(
            head
            + conflict_html
            + "<table><thead><tr>"
            "<th>Sıra</th><th>Eγ (keV)</th><th>Iγ</th><th>ENSDF Iγ</th>"
            "<th>Normlama</th><th>Belirsizlik</th><th>Bozunma</th><th>Not</th>"
            "</tr></thead><tbody>"
            + "".join(rows)
            + "</tbody></table>"
            + note
            + "</section>"
        )

    body = "\n".join(blocks) if blocks else "<p>Bu enerji limitinde aday yok.</p>"
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    gamma_cut = (
        f"ek kesim Eγ ≥ {min_gamma_kev:g} keV"
        if min_gamma_kev > 0
        else "atomik X-ışınları elendi"
    )
    return f"""<!DOCTYPE html>
<html lang="tr">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Photonuclear katalog — {esc(result.target)}</title>
<style>
:root {{
  --bg: #f6f1e8;
  --ink: #1c1917;
  --card: #fffdf8;
  --line: #d6cfc4;
  --accent: #0f766e;
  --muted: #78716c;
  --warn: #9a3412;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  color: var(--ink);
  background:
    radial-gradient(circle at top left, #e7f0ea, transparent 40%),
    linear-gradient(180deg, #f8f4ec, var(--bg));
  line-height: 1.45;
}}
main {{ max-width: 960px; margin: 0 auto; padding: 2rem 1.25rem 3rem; }}
h1 {{ font-size: 1.8rem; margin: 0 0 0.4rem; letter-spacing: -0.02em; }}
h3 {{ margin: 0 0 0.6rem; color: var(--accent); }}
.meta, .sources, .note {{ color: var(--muted); margin-bottom: 1.25rem; }}
.meta div, .sources li, .note li {{ margin: 0.15rem 0; }}
.card {{
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 14px;
  padding: 1rem 1.1rem 1.1rem;
  margin: 1rem 0;
  box-shadow: 0 8px 24px rgba(28, 25, 23, 0.04);
}}
ul {{ padding-left: 1.1rem; margin: 0.2rem 0 0.8rem; }}
table {{ width: 100%; border-collapse: collapse; font-size: 0.92rem; table-layout: fixed; }}
th, td {{
  border-bottom: 1px solid var(--line);
  padding: 0.45rem 0.55rem;
  text-align: left;
  vertical-align: top;
}}
th {{ font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); }}
/* Keep digit widths consistent; do NOT right-align — that pulls values away from headers */
td:nth-child(1), td:nth-child(2), td:nth-child(3), td:nth-child(4) {{
  font-variant-numeric: tabular-nums;
}}
th:nth-child(1), td:nth-child(1) {{ width: 3.5rem; }}
th:nth-child(2), td:nth-child(2) {{ width: 7.5rem; }}
th:nth-child(3), td:nth-child(3) {{ width: 5.5rem; }}
th:nth-child(4), td:nth-child(4) {{ width: 7rem; }}
.warn {{ color: var(--warn); }}
.muted {{ color: var(--muted); font-size: 0.9rem; }}
code {{ background: #efeae2; padding: 0.1rem 0.35rem; border-radius: 6px; }}
sup {{ font-size: 0.75em; }}
</style>
</head>
<body>
<main>
  <h1>Photonuclear kinematik aday kataloğu — {esc(result.target)}</h1>
  <div class="meta">
    <div><strong>Enerji limiti:</strong> {result.emax_mev:g} MeV</div>
    <div><strong>Aday ürün sayısı:</strong> {len(result.products)}</div>
    <div><strong>σ_negligible eşiği:</strong> {negl:g} mb</div>
    <div><strong>Gama:</strong> mutlak/tespitli Iγ; {esc(gamma_cut)}</div>
    <div><strong>Oluşturma:</strong> {esc(stamp)}</div>
  </div>
  <div class="note">
    <strong>Okuma notu:</strong>
    <ul>
      <li>Eth kinematik eşiktir; Eth+Vc yalnızca kaba bariyer göstergesidir.</li>
      <li>TENDL artıksal kesit varsa Coulomb göstergesinden önceliklidir.</li>
      <li>TENDL eşlemesi MT + IZAP + LFS ile yapılır (varsayılan kaynak s60).</li>
      <li>LFS=0 temel durum; MF3 kanal toplamı; Σ_LFS seviye toplamıdır.</li>
      <li>Kanal toplamı olup seviye ayrımı yoksa channel_total_only;
          çok küçük kesitler sigma_negligible olarak işaretlenir.</li>
      <li>Iγ öncelikle mutlak ENSDF değeridir; göreli dal şiddeti gerektiğinde
          dal oranı ile ölçeklenir.</li>
      <li>Atomik X-ışınları nükleer γ listesinden ayrılır.</li>
    </ul>
  </div>
  <div class="sources">
    <strong>Kaynaklar:</strong>
    <ul>
      <li>Referans arayüz: NuDat 3</li>
      <li>Eşikler: AME2020 · Yarıömür/izomer: NUBASE2020</li>
      <li>Gamalar: ENSDF (IAEA LiveChart)</li>
      <li>σ: TENDL-2023 MF3 + MF10 (MT+IZAP+LFS)</li>
    </ul>
  </div>
  {body}
  <p class="muted">Aynı klasörde <code>products.csv</code>, <code>gammas.csv</code>,
  <code>catalog.json</code> ve <code>report.md</code> dosyaları da vardır.</p>
</main>
</body>
</html>
"""


def _real_gamma_lines(gammas: list["GammaLine"]) -> list["GammaLine"]:
    out: list[GammaLine] = []
    for g in gammas:
        if g.energy_kev != g.energy_kev:  # NaN placeholder
            continue
        if g.status in {
            "network_error",
            "no_matched_ensdf_parent",
            "no_absolute_decay_gamma_dataset",
            "no_gamma",
            "xray_only",
            "cut_excluded",
        } and g.rank is None:
            continue
        out.append(g)
    out.sort(
        key=lambda g: (
            0 if g.rank is not None else 1,
            g.rank if g.rank is not None else 999,
            g.energy_kev,
        )
    )
    return out


def _group_gammas(gammas: list["GammaLine"]) -> dict[str, list["GammaLine"]]:
    by: dict[str, list[GammaLine]] = defaultdict(list)
    for g in gammas:
        by[g.product].append(g)
    for key in by:
        by[key].sort(
            key=lambda g: (
                0 if g.rank is not None else 1,
                g.rank if g.rank is not None else 999,
                g.energy_kev if g.energy_kev == g.energy_kev else 1e99,
            )
        )
    return by


def _fmt_rank(value: int | None) -> str:
    return "—" if value is None else str(value)


def _fmt_energy(value: float) -> str:
    if value != value:
        return "—"
    if value >= 100:
        return f"{value:.1f}"
    if value >= 10:
        return f"{value:.2f}"
    return f"{value:.3f}"


def _fmt_intensity(value: float | None, unc_or_status: object) -> str:
    if value is None:
        return "yok (ENSDF)"
    if isinstance(value, float) and value != value:
        return "—"
    if isinstance(unc_or_status, str) and unc_or_status == "intensity_missing":
        return "yok (ENSDF)"
    if value >= 10:
        return f"{value:.2f}"
    if value >= 1:
        return f"{value:.3f}"
    if value >= 0.01:
        return f"{value:.4f}"
    return f"{value:.3g}"


def _fmt_optional(value: float | None) -> str:
    if value is None:
        return "—"
    return _fmt_intensity(value, "")


def _fmt_optional_num(value: float | None) -> str:
    if value is None:
        return "—"
    return f"{value:g}"


def _fmt_sigma(product: "Product") -> str:
    if product.sigma_max_mb is None:
        return f"— ({product.sigma_status})"
    epart = (
        "—"
        if product.e_at_sigma_max_mev is None
        else f"{product.e_at_sigma_max_mev:g} MeV"
    )
    return f"{product.sigma_max_mb:g} mb @ {epart} ({product.sigma_status})"


def _estimate_flags(product: "Product") -> str:
    flags = []
    if product.mass_estimated:
        flags.append("kütle #")
    if product.half_life_estimated:
        flags.append("T½ #")
    if product.excitation_estimated:
        flags.append("Exc #")
    return ", ".join(flags) if flags else "yok"


def _state_label(product: "Product") -> str:
    if product.is_stable:
        return "kararlı gs"
    if product.is_metastable:
        return f"metastabil ({product.isomer_display or 'm'})"
    return "radyoaktif gs"


def _fmt_norm(value: str | None) -> str:
    mapping = {
        "absolute": "mutlak",
        "branch_relative": "dal göreli → mutlak",
        "relative_adopted_level": "göreli (Adopted Levels)",
        "uncertain": "belirsiz",
        "not_applicable": "—",
    }
    if not value:
        return "—"
    return mapping.get(value, value)


def _fmt_unc_status(value: str | None) -> str:
    mapping = {
        "intensity_missing": "Iγ yok",
        "branch_unc_not_propagated": "dal belirsizliği taşınmadı",
        "status_row": "—",
        "": "—",
    }
    if value is None:
        return "—"
    return mapping.get(value, value)


def _gamma_status_message(gammas: list["GammaLine"]) -> str:
    if not gammas:
        return "Uygun nükleer γ kaydı bulunamadı."
    g = gammas[0]
    status = g.status or ""
    detail = g.note or g.decay_mode or ""
    mapping = {
        "no_matched_ensdf_parent": "Eşleşen ENSDF üst seviyesi bulunamadı.",
        "no_gamma": "Eşleşen ENSDF üst seviyesi bulunamadı.",
        "no_absolute_decay_gamma_dataset": (
            "Üst seviye Adopted Levels’ta var; mutlak bozunma γ veri seti ve "
            "listelenebilir göreli geçiş bulunamadı."
        ),
        "relative_adopted_level": (
            "Çizgiler Adopted Levels göreli şiddetleridir (mutlak Iγ% değil)."
        ),
        "xray_only": (
            "Üst seviye eşleşti; nükleer γ bulunmadı, yalnızca atomik X-ışınları vardı."
        ),
        "intensity_missing": "γ çizgisi var; Iγ ENSDF’te verilmemiş.",
        "cut_excluded": "Çizgiler enerji kesimi nedeniyle listelenmedi.",
        "network_error": "Gama verisi alınamadı (ağ/veri kaynağı).",
        "parent_match_conflict": "Üst seviye eşleşti; yarıömür değerleri çelişiyor.",
        "state_specific_missing": "TENDL’de izomere özgü seviye kaydı yok.",
    }
    base = mapping.get(status, "Gama verisi bulunamadı.")
    # Prefer the clean Turkish message; skip raw English implementation notes.
    if status in {
        "xray_only",
        "relative_adopted_level",
        "no_absolute_decay_gamma_dataset",
        "no_matched_ensdf_parent",
        "no_gamma",
        "network_error",
        "cut_excluded",
        "intensity_missing",
        "parent_match_conflict",
        "state_specific_missing",
    }:
        return base
    if detail and detail not in base and "ERROR" not in detail.upper():
        # Avoid dumping snake_case status tokens into the report.
        if "_" in detail and detail.replace("_", "").isalpha():
            return base
        return f"{base} {detail}"
    return base


def _gamma_note(g: "GammaLine") -> str:
    if g.status == "network_error" or "ERROR" in (g.decay_mode or ""):
        return "veri alınamadı"
    if g.status == "xray_only":
        return "yalnızca atomik X-ışınları"
    if g.normalization_status == "relative_adopted_level" or g.status == "relative_adopted_level":
        return "göreli Adopted Levels"
    if g.uncertainty_status == "intensity_missing" or g.status == "intensity_missing":
        return "Iγ ENSDF’te yok"
    if g.parent_match_conflict:
        return "yarıömür çelişkisi"
    # Hide raw implementation notes from the human report.
    note = (g.note or "").strip()
    if not note or "relative_adopted_level" in note or "parent_match_conflict" in note:
        return "—"
    if "branch_relative" in note:
        return "dal oranı ile ölçeklendi"
    if "absolute ENSDF" in note:
        return "—"
    if "normalization uncertain" in note:
        return "normlama belirsiz"
    return note if len(note) < 80 else "—"
