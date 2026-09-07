# TASK1 — Photonuclear Product Catalog

Kararlı bir hedef çekirdekten, fotonnükleer `(γ, xn yp)` tepkimeleriyle **45 MeV** (veya daha düşük `E_max`) altında **kinematik olarak açık** radyoaktif izotopları ve metastabil durumları listeler; Coulomb bariyeri + TENDL-2023 artıksal σ ile sınıflar; her biri için en güçlü en fazla 10 **mutlak** bozunma gamasını çıkarır.

Kaynak referansı: [NuDat 3](https://www.nndc.bnl.gov/nudat3/).  
Plan: [`PLAN.md`](PLAN.md) · Kolay kullanım: [`KULLANIM.md`](KULLANIM.md)

## En kolay kullanım (Windows)

1. Bu klasörü zip olarak gönderin / indirin  
2. `run.bat` dosyasına çift tıklayın  
3. Hedef (`197Au`) ve enerji (`44`) girin  
4. `out\...\report.html` dosyasını tarayıcıda açın  

Linux/macOS/WSL: `./run.sh`

## Ne üretir?

Her çalıştırmada çıktı klasöründe:

| Dosya | Açıklama |
|-------|----------|
| **`report.html`** | İnsan okur özet (önerilen) |
| `report.md` | Aynı özet Markdown |
| `products.csv` | Ürün listesi (`Eth`, `E_eff`, `σ_max`, `viability`) |
| `gammas.csv` | Mutlak Iγ tablosu |
| `catalog.json` | Tam JSON |

### Fizik düzeltmeleri (v0.2.4)

- **Eth kinematiktir**; ölçülebilir verim iddiası değildir.
- **TENDL kaynağı:** varsayılan **s60** (explicit kanallar ≤60 MeV); `sigma_source_variant/url/sha256`.
- **TENDL anahtar:** `MT + IZAP + LFS` (serbest nükleon MT; küme kanalı yedeklenmez).
- **MF3/MF10:** `sigma_channel_total_mb` + `E_at`; `sigma_sum_mf10_lfs_mb`; `mf10_coverage_ratio`.
- MF3 var / MF10 yok → `channel_total_only` / `channel_supported_state_split_missing`.
- Küçük ama mevcut σ → `sigma_negligible` (yanlış `mt_not_in_tendl` değil).
- CSV: `no_matched_ensdf_parent`; `normalization_status=not_applicable` status satırlarında.

## Komut satırı

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e .

python -m photonuclear_catalog --interactive
python -m photonuclear_catalog --target 197Au --emax 44
python -m photonuclear_catalog --target 208Pb --emax 20 --no-tendl
```

İlk TENDL indirmesi internet ister; dosyalar `data/tendl_cache/` altına yazılır.

## Kaynak notu

NuDat 3’ün public API’si olmadığı için hesap motoru AME2020 + NUBASE2020 + ENSDF (LiveChart) + TENDL-2023 kullanır. Kritik γ çizgileri NuDat decay çıktısıyla uyumludur; Iγ mutlak normalizasyonu NuDat’taki “absolute intensity” ile aynı mantıktadır.

## Test

```bash
pip install -e ".[dev]"
pytest -q
```
