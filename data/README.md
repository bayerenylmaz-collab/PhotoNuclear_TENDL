# Bundled evaluated nuclear data

These ASCII tables are included so the catalog tool works offline for
threshold and half-life/isomer logic.

| File | Evaluation | Use in this project |
|------|------------|---------------------|
| `mass_1.mas20.txt` | AME2020 | Atomic mass excess → Eth for (γ, xn yp) |
| `nubase_4.mas20.txt` | NUBASE2020 | Stability, half-life, metastable states |

## Citations

- W.J. Huang et al., Chin. Phys. C **45**, 030002 (2021); M. Wang et al., Chin. Phys. C **45**, 030003 (2021) — AME2020
- F.G. Kondev et al., Chin. Phys. C **45**, 030001 (2021) — NUBASE2020

Download source: https://www-nds.iaea.org/amdc/

## Relation to NuDat 3

The task reference UI is [NuDat 3](https://www.nndc.bnl.gov/nudat3/). NuDat is backed by ENSDF / Wallet Card / mass evaluations. For programmatic use we:

1. compute thresholds from AME2020,
2. classify radioactive / metastable states with NUBASE2020,
3. pull decay gammas from the IAEA LiveChart API (`fields=decay_rads`, ENSDF), which is NuDat-compatible evaluated decay radiation data.
