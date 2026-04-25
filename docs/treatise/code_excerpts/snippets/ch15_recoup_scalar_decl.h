// ===============================================================
// 抽取自仓库 [current]: src/utils/chiral_3nf_recoupling.h
// 行号区段：22..55
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
// Rank-0 (scalar pair operator) 3NF recoupling coefficient.
//
// Applies to operators of the form
//   O_scalar = (spatial scalar) * (spin scalar) * (isospin scalar)
// acting on the pair (2,3) with the spectator (1) unchanged. The three
// scalar factors include:
//   - spin:    sigma_2 . sigma_3      (eigenvalue 2 S_2N (S_2N+1) - 3)
//   - isospin: tau_2 . tau_3          (eigenvalue 2 T_2N (T_2N+1) - 3)
//   - spatial: any channel-diagonal scalar (pair orbital unchanged)
//
// The returned coefficient is the Kronecker product of selection rules
// times the pair spin and isospin eigenvalues. Since the operator is
// scalar in every subspace, the recoupling reduces to a product of
// delta-functions times the algebraic eigenvalues (all other recoupling
// 6j's collapse to 1/(2J+1) normalisations absorbed into the convention).
//
// This form matches A_c1_rank0 in formula_reference.md §3.4 and is also
// the limiting case of A_cE in §1.4 and A_c3_rank0 in §4 (up to the
// overall LEC- and momentum-kernel dependent prefactors handled by the
// caller).
//
// Arguments (row = bra "r", column = ket "c"):
//   L_2N_*, S_2N_*, J_2N_*, T_2N_* : pair partial-wave integers
//   L_1N_*, two_J_1N_*             : spectator orbital integer, 2*spectator total (half-integer)
//   two_J_3N, two_T_3N             : 2*3N conserved quantum numbers
//
// Returns: real-valued recoupling coefficient (no regulator, no LEC, no
// momentum kernel).
double recoupling_3nf_scalar(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N);
