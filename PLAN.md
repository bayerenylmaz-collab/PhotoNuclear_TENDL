# TASK1 — Kilitlenen Plan

Bu belge, Cursor ile birlikte tartışıp üzerinde anlaştığımız kapsamı ve çalışma planını kayda geçirir. Bundan sonra uygulamayı bu plana göre yürütürüz.

## Benden istenen iş (özet)

Kararlı bir çekirdekten yola çıkarak, 45 MeV’nin altında parçacık ayırma / üretim eşiğine sahip olabilecek radyoaktif izotopları ve metastabil durumları listelemek. `(γ,n)`, `(γ,2n)`, …, `(γ,p)`, `(γ,2p)`, …, `(γ,pn)`, … gibi fotonnükleer tepkimeler dikkate alınacak. Liste için her izotopun en fazla ilk/en güçlü 10 bozunma gamasını çıkarmak. Ardından bana 45 MeV’den daha düşük bir enerji limiti verdiğimde o listeden daha dar bir alt küme seçebilecek bir kod yazmak. Kaynak: [NuDat 3](https://www.nndc.bnl.gov/nudat3/).

## Ne tartıştık, neye karar verdik?

### 1) Hangi kararlı çekirdekler?

İşi “birkaç çekirdeği elle seçip sabit liste üretmek” olarak değil, **kullanıcının verdiği kararlı çekirdek(ler) için çalışan genel bir araç** olarak okuyoruz.

- Kod, hedef çekirdeği parametre olarak alır.
- Böylece “neden bu çekirdek?” sorusu kodun sorumluluğu olmaz; kullanıcı seçer, biz fiziği doğru ve tekrarlanabilir uygularız.
- README’de birkaç örnek hedefle (ör. `208Pb`, `63Cu`) çalıştırma göstereceğiz; bu örnekler tercihten ziyade doğrulama / demo içindir.

### 2) Hangi tepkime kanalları?

Metne sadık kalıyoruz:

- Dahil: `(γ, xn)`, `(γ, yp)`, `(γ, xn yp)` — yani yalnızca nötron ve proton koparma kombinasyonları.
- Bilinçli olarak özel kanal işlemeyenler: `(γ,α)`, `(γ,d)`, `(γ,t)` vb. (metinde yok).
- Not: `(γ,2p2n)` kütlece alfa ile aynı parçacık sayısına denk gelebilir; biz bunu **n+p kombinasyonu** olarak hesaplarız, ayrı bir alfa kanalı açmayız.

### 3) Metastabil durum

Metastabil / izomer: aynı izotopun (aynı Z, N) uzun ömürlü uyarılmış hali (`99mTc`, `180mTa` gibi).

Listede:

- fotonnükleer yolla ulaşılabilen **radyoaktif ürünler**, ve
- NuDat / NUBASE’de bilinen **metastabil (izomerik) durumları**

yer alacak. Üretim eşiği pratikte çoğunlukla temel durumla aynı kabul edilir; izomerler ayrı satır olur çünkü bozunma gamaları farklı olabilir.

### 4) “En güçlü 10 gama” kriteri

Varsayılan sıralama: **γ şiddeti (`Iγ`) azalan**. En fazla 10 çizgi.

- Doz veya E×I gibi alternatifleri bilinçli kullanmıyoruz; görev spektroskopik “en güçlü çizgi” diyor, `Iγ` bunun karşılığı.
- Birden fazla bozunma modu varsa, ilgili izotopun bozunma radyasyonlarından γ’ları toplayıp `Iγ`’ya göre sıralarız.

### 5) Enerji filtresi

- Ana tarama tavanı: **45 MeV**.
- Sonra kullanıcı `E_max < 45` verdiğinde, aynı sonuç tablosundan `Eth ≤ E_max` olanlar seçilir.
- Yani önce 45 MeV’ye kadar üret → sonra yerel filtre. Bu, istenen 3. maddeyle birebir örtüşür.

### 6) Çıktı formatı

**Ana ürün: Python CLI + CSV/JSON.**

- Notebook zorunlu değil; istenirse sonradan ince bir örnek eklenebilir.
- Tipik kullanım:

```bash
python -m photonuclear_catalog --target 208Pb --emax 45
python -m photonuclear_catalog --target 208Pb --emax 30
```

- Çıktılar: ürün listesi + her ürün için top-10 gamalar (`CSV` / `JSON`).

### 7) Yazılım yığını (bilgisayarımda ne var?)

Bilgisayarımda WSL, Fortran, MATLAB, ROOT ve C++ var; VS Code üzerinden Python’u da kurabilirim.

Bu görev için kararımız:

- **Python** kullanıyoruz (veri işleme, CLI, CSV/JSON, NuDat/AME-NUBASE ile çalışma için en uygun yol).
- ROOT / Fortran / MATLAB / C++ bu işte zorunlu değil; spektrum analizi veya FLUKA simülasyonu değil, nükleer veri kataloğu + filtre aracı yazıyoruz.
- Ortam: WSL + Python 3 + `pip` yeterli.

## Fizik kuralı (uygulamada kullanacağımız)

Hedef çekirdek `T = (Z, N)` için kanal `(γ, xn yp)`:

\[
E_{\mathrm{th}} = \big[M(Z-y,\,N-x) + x\,M_n + y\,M_p - M(Z,N)\big]c^2
\]

`Eth ≤ E_max` ise ürün **kinematik adaydır**. Kararlı (kalıcı) ürünleri eleeriz; radyoaktif gs + metastabil durumları tutarız.

Ek fizik katmanları (v0.2):

- `E_eff ≈ Eth + Vc` — yüklü emisyon için Coulomb bariyeri tahmini.
- TENDL-2023 MF10 artıksal σ — `[Eth, Emax]` tepe değeri; `viability` sınıflaması.
- Mutlak Iγ = `Iγ_rel × BR`; X-ışınları karakteristik enerjiyle elenir (sert 100 keV kesimi yok).

## Veri kaynakları

- **NuDat 3** — istenen referans arayüz / bozunma ve yapı verisi.
- **AME2020** — kütleler ve eşik / ayırma enerjisi hesabı için.
- **NUBASE2020** — yarı ömür, kararlılık, izomer (`m1/m2/…`) bilgisi için.
- **TENDL-2023** — fotonnükleer artıksal üretim kesitleri (MF10).

NuDat’ın düzgün bir public REST API’si yok; bu yüzden NuDat ile uyumlu değerlendirilmiş tabloları (AME/NUBASE/ENSDF LiveChart) + TENDL kullanırız. Kaynaklar README’de açıkça yazılacak.

## Uygulama adımları

1. [x] Python proje iskeleti: paket, CLI, bağımlılıklar, CSV/JSON yazıcı.
2. [x] Kütle tablosu + `Eth(γ, xn yp)` motoru; verilen hedef ve `E_max` için ürün listesi.
3. [x] Radyoaktif + metastabil filtre.
4. [x] Her ürün için en güçlü ≤10 bozunma gamması (mutlak `Iγ`).
5. [x] `E_max` yerel filtresi; örnek çalıştırma; kısa testler; dokümantasyon.
6. [x] Coulomb + TENDL viability; X-ray / abs Iγ / m1–m3 gösterim düzeltmeleri.

## Bilinçli kapsam dışı

- `(γ,α)`, `(γ,d)`, `(γ,t)` özel kanalları.
- Tam TALYS yeniden koşumu / izomer besleme oranlarının deneysel kalibrasyonu (TENDL artıksal σ kullanılır; izomer LFS varsa tercih edilir).
- FLUKA / ROOT spektrum analizi (önceki işlerimiz; bu repo farklı).
- Sabit, tek seferlik “seçilmiş çekirdek kataloğu” (araç genel olacak).

---

*Bu plan Cursor ile birlikte kilitlenmiştir. Uygulama bu belgeye göre ilerler; sapma olursa önce plan güncellenir.*

## Bilinen teknik not (gama)

LiveChart’ta `207mpb` / `60mco` gibi **izomer id’leri boş** döner. Gamalar gs id üzerinden çekilip `p_energy` ≈ NUBASE uyarılma enerjisi ile ayrılır. Bazı IT çizgilerinde ENSDF `Iγ` vermez; bunlar listede tutulur ve not düşülür. ENSDF’te hiç γ kaydı olmayan nüklitler boş kalabilir (raporda ayrı durum mesajı).

Atomik X-ışınları karakteristik K/L enerjisiyle elenir; isteğe bağlı `--min-gamma-kev` ekstra kesim ekler.

Teslimde `report.html` / `report.md` özet raporları ve Windows `run.bat` kolay çalıştırıcı vardır (`KULLANIM.md`).
