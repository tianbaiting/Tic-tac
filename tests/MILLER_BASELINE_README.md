# Miller 2NF Baseline Guardrail

This directory holds the bit-for-bit guardrail that pins Tic-tac's
2NF-only code path (`three_nucleon_force=none`) before the upcoming
chiral N2LO 3NF rewrite.

## Files

| Path                                              | Role                                                                  |
|---------------------------------------------------|-----------------------------------------------------------------------|
| `tests/test_2nf_miller_baseline.py`               | Regression test — loads the hash JSON, asserts no drift.              |
| `tests/fixtures/freeze_2nf_baseline.py`           | One-shot fixture script — runs the solver and writes the hash JSON.   |
| `tests/fixtures/2nf_baseline_hashes.json`         | Frozen SHA-256 hashes of every `U_PW_elements_*.txt` (checked in).    |

The guardrail exists because the Faddeev solver runs for hours; we
cannot re-solve the problem in CI on every commit. Instead we freeze a
digest of the 2NF output once and compare each future run against it.

## Benchmark sources (cited in the hash file)

1. Miller, Ekström, Hebeler, **Phys. Rev. C 106, 024001 (2022)** —
   primary upstream WPCD paper. Target values:
   - 13 MeV nd elastic phase shifts (Nijmegen-I) vs. Glöckle–Witała
     *standard Faddeev*, within ~1%.
   - Ay(n) minimum ≈ −0.50 ± 0.02 at θ_cm ≈ 120° for E_lab = 35 MeV,
     Nijmegen-I.
2. Glöckle, Witała, Huber, Kamada, Golak, **Phys. Rep. 274, 107
   (1996)** — the "standard Faddeev" reference Miller benchmarks against.
3. Nogga, Kamada, Glöckle, Barrett, **Phys. Rev. C 65, 054003 (2002)**
   — Idaho-N3LO triton binding ≈ −7.855 MeV (NN-only), the community
   triton-E_B target.

## How to re-freeze (rarely — see below)

```bash
python3 tests/fixtures/freeze_2nf_baseline.py \
    --work-dir output/2nf_miller_baseline \
    --target-tlab-mev 190
```

This overwrites `tests/fixtures/2nf_baseline_hashes.json`.

## When to re-freeze

**Essentially never.** Re-freezing the baseline is equivalent to
removing the guardrail. Re-freeze only when:

1. The 2NF physics of the solver is *intentionally* changed (e.g. a
   corrected isospin factor in the 2N potential, a fixed quadrature
   bug), AND
2. The change has been independently reviewed against Miller PRC 106
   Fig. 1 (13 MeV Nijmegen-I phase shifts) and Fig. 2 (35 MeV Ay(n))
   and still agrees within the ~1% / ±0.02 envelopes respectively.

Re-freezing to silence a failing `test_2nf_miller_baseline.py` is
exactly the mistake this guardrail is designed to prevent. If the
guardrail fires during 3NF-rewrite work, the 3NF code is leaking into
the `three_nucleon_force=none` path.

## Known gaps vs. Miller (follow-up TODOs, not in scope for this guardrail)

* Miller Fig. 1 / Fig. 2 values are read off plotted curves rather
  than tabulated; digitization of `tools/check_3nf_normalization/
  miller_benchmark.md` into numeric references would convert the
  `test_nijmegen_i_phase_shift_reference` skip stub into a real
  assertion.
* `CPP/run` currently targets elastic nd scattering above threshold;
  reproducing the Nogga triton binding requires a separate bound-state
  driver, which is why `test_triton_binding_nogga_nijmegen` is
  decorated `@unittest.skip`.
