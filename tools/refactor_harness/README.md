# Bitwise regression harness

Guards **byte-for-byte stability** of deterministic solver artifacts across the
structural refactor. Two reduced-grid solves cover the full pipeline:

| Case | Input | Covers |
|------|-------|--------|
| 2NF-only | `input_golden_2nf.txt`  | P123, V, SWP, resolvent, Neumann, Padé, on-shell U |
| complete-3NF | `input_golden_3nf.txt` | all of the above **plus** the exact factorized W^(1) build + AGS kernel `K = P V + (1+P) W1` |

Both use the documented reduced grid `Np_WP=4, Nq_WP=3, J_2N_max=1,
two_J_3N_max=1, channel_idx=0` (single `JP=1/2+` block) with `potential_model=N2LOopt`.

## Usage

```bash
# Re-run both solves and compare every text artifact to the baseline
tools/refactor_harness/run_and_compare.sh

# (Re)write the baseline after an intentional, verified change
tools/refactor_harness/run_and_compare.sh --record
```

The script hashes every non-HDF5 output (U elements, convergence sidecar,
Neumann history, run parameters, q kinematics/boundaries) and `diff`s against
`baseline_hashes.txt`. HDF5 P123/W1 caches are excluded because their internal
byte layout is not a contract; their *content* parity is covered by the CTest
`cache_layer_test` and `chiral_n2lo_w1_cache` suites.

## Determinism note

Both solves are bitwise **self-deterministic** across fresh-cache reruns on this
machine (verified: two independent runs produce identical hashes). The stored
`tmp/factorized_conv/*` artifacts are **not** a valid baseline here: they
predate the cE-normalization / W1-schema-v8/v9 physics fixes and therefore
differ from current-code output by design.
