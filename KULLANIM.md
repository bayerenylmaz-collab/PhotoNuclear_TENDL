# Kullanım Kılavuzu (e-posta ile paylaşım)

Bu klasörü zipleyip gönderebilirsiniz. Karşı taraf **Python 3.10+** kurulu bir bilgisayarda çok az adımla çalıştırır.

## Windows (en kolay)

1. Zip’i açın.
2. **Python 3.10+ kurulu olmalı**  
   - https://www.python.org/downloads/  
   - Kurulumda **“Add python.exe to PATH”** işaretleyin  
   - Microsoft Store’daki boş “python” kısayolu yetmez
3. `run.bat` dosyasına **çift tıklayın**.
4. Sorulan sorulara cevap verin:
   - Hedef çekirdek: örn. `197Au` veya `208Pb`
   - Enerji limiti: örn. `44` veya `45`
   - TENDL-2023 σ: genelde `e` (evet; internet gerekir)
5. Bitince `out\...` klasöründe sonuçlar oluşur.
6. **`report.html`** dosyasını çift tıklayıp tarayıcıda okuyun.

İlk çalışmada paket kurulumu 1–2 dakika sürebilir. İnternet gerekir (bozunma gamaları + isteğe bağlı TENDL).

`Python was not found` görürseniz: gerçek Python kurulu değildir. Yukarıdaki 2. adımı yapıp `run.bat`’i yeniden açın.

## Linux / WSL / macOS

```bash
chmod +x run.sh
./run.sh
```

veya:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python -m photonuclear_catalog --interactive
```

## Komut satırı (ileri)

```bash
python -m photonuclear_catalog --target 197Au --emax 44
python -m photonuclear_catalog --target 63Cu --emax 45 --out out/cu63_45
python -m photonuclear_catalog --target 208Pb --emax 20 --no-tendl
```

## Çıktı dosyaları

| Dosya | Ne işe yarar |
|-------|----------------|
| `report.html` | İnsan okur özet (önerilen) |
| `report.md` | Aynı özet, Markdown |
| `products.csv` | Ürün listesi: Eth, E_eff, σ_max, viability |
| `gammas.csv` | En güçlü **mutlak** Iγ çizgileri |
| `catalog.json` | Tam makine-okur çıktı |

Her ürün için en fazla **10** bozunma gamasi, **mutlak Iγ** (`Iγ × BR`) azalan sırada listelenir. Atomik X-ışınları karakteristik enerjiyle elenir; varsayılan sert `Eγ ≥ 100 keV` kesimi yoktur.

`Eth` yalnız kinematik adaylıktır. Ölçülebilirlik için `E_eff`, TENDL `σ_max` ve `viability` alanlarına bakın.

## Kaynak notu

- Referans arayüz: [NuDat 3](https://www.nndc.bnl.gov/nudat3/)
- Eşik hesabı: AME2020
- Yarı ömür / izomer: NUBASE2020
- Gamalar: ENSDF (IAEA LiveChart; mutlak Iγ)
- Artıksal σ: TENDL-2023 MF10

## Sorun giderme

- **Python yok:** https://www.python.org/downloads/ — kurulumda PATH seçeneğini işaretleyin.
- **Gama satırları az/boş:** interneti kontrol edin; bazı izotoplarda ENSDF’te γ yoktur (raporda ayrı durum mesajı çıkar).
- **TENDL indirilemedi:** `--no-tendl` ile yalnız kinematik + Coulomb çıktısı alınır.
- **WSL kullanıyorsanız:** proje Windows masaüstündeyse `/mnt/c/Users/...` yolundan çalıştırın.
