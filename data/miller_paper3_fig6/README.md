# Miller Paper III Fig 6 — Digitized Reference

Source: S.B.S. Miller, A. Ekström, C. Forssén,
"Posterior predictive distributions of neutron-deuteron cross sections",
**Phys. Rev. C 107, 014002 (2023)** — Fig. 6 (p. 7).

Also appears as Fig. 4.2 in Miller's PhD thesis (Chalmers 2022).

## Content

PPDs for elastic nd neutron vector analyzing power A_y(n) at three
laboratory energies, four chiral orders, plus EXFOR experimental markers.

- `Ay_n_N3LO_Elab10MeV.csv` — N³LO (red) curve at E_Lab = 10 MeV
- `Ay_n_N3LO_Elab35MeV.csv` — N³LO (red) curve at E_Lab = 35 MeV
- `Ay_n_N3LO_Elab67MeV.csv` — N³LO (red) curve at E_Lab = 67 MeV
- `Ay_n_N2LO_Elab10MeV.csv` — N²LO (blue) at 10 MeV (Ay puzzle context)
- `Ay_n_NLO_Elab10MeV.csv`  — NLO (green) at 10 MeV
- `Ay_n_expt_Elab10MeV.csv` — EXFOR markers at 10 MeV (Bunker 1968)
- `Ay_n_expt_Elab35MeV.csv` — markers at 35 MeV (Binder/LENPIC 2016)
- `Ay_n_expt_Elab67MeV.csv` — markers at ~66.6 MeV (Bunker 1968)

## Column format

```
theta_cm_deg, Ay_value[, Ay_error]
```

`Ay_error` is present only for experimental files.

## Accuracy status

**Draft values from visual reading of the PDF figure.**
Precision ≈ ±0.01 in A_y, ±2° in θ_cm.

**TODO**: refine using WebPlotDigitizer (https://apps.automeris.io/wpd/)
or similar tool once we have network access. For experimental markers, the
canonical source is the EXFOR database entry for Bunker 1968
(Nucl. Phys. A 113, 461).

## Key Ay puzzle observation (from Paper III §IV.D)

> "At E_Lab = 10 MeV we do not reproduce the experimental data at any
> chiral order. ... the N³LO prediction has significantly lower accuracy
> at lower energies."

- N³LO peak: A_y ≈ 0.12 at θ_cm ≈ 120°
- EXFOR peak: A_y ≈ 0.19 at θ_cm ≈ 120°
- Puzzle magnitude: ΔA_y ≈ 0.07 (absolute), ≈ 58% relative

This gap is the signal 3NF is supposed to close.
