# 2NF-only Np=Nq=5 regression baseline

This baseline uses `potential_model=N2LOopt`, `three_nucleon_force=none`,
`J_2N_max=1`, `two_J_3N_max=1`, and `CPP/Input/lab_energies.txt`.

It was refreshed on 2026-08-11 after fixing the production Padé convergence
selector. The old selector could stop at `P[0/0]` and read order `-1`; simply
taking `P[14/14]` was also poor for this asymptotic sequence. The corrected
solver evaluates every order through 14, certifies convergence only from the
final three-order tail, and uses the closest finite consecutive pair only as a
labelled optimal truncation when the final tail is not stable.

Independent validation used the same kernel and grid with `solve_dense=true`.
Across both parities and all 40 stored complex amplitudes,

```text
max |U_Pade - U_dense| = 6.57420273e-09 MeV
mean|U_Pade - U_dense| = 4.82005695e-10 MeV
```

Every element is conservatively marked `Conv=2` because the final three Padé
orders are not stable, even though the selected optimal truncations agree with
the dense oracle. The sidecar reports convergence status, while these files
lock the independently validated numerical 2NF result.

Reproduction commands are listed in `docs/3nf_final_audit.md`.
