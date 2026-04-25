// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：1..50
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
#ifndef CHIRAL_3NF_PW_KERNELS_H
#define CHIRAL_3NF_PW_KERNELS_H

#include <cmath>

// ============================================================================
// [EN] Scalar momentum-space kernels for the partial-wave decomposition of the
//      chiral N2LO three-nucleon force.
//
//      These helpers factor out the pure *momentum-space* pieces (pion
//      propagators, LEC prefactors, regulator) so that the outer partial-wave
//      machinery (recoupling coefficients, spin/isospin matrix elements,
//      x-integration) can stay clean. All functions are side-effect free,
//      header-only, and inline.
//
//      Formula sources (transcribed in
//      tools/check_3nf_normalization/formula_reference.md):
//        - [E2002] Epelbaum et al., PRC 66 (2002) 064001 — eqs. (2.10), (3.19),
//                  App. A eq. (A-2).
//        - [G2010] Golak et al., EPJA 43 (2010) 241 — eqs. (17), (18), (22).
//        - [H2015] Hebeler et al., PRC 91 (2015) 044001 — eqs. (13)–(15).
//
//      Conventions (see §5 of formula_reference.md):
//        - Jacobi momenta `p`, `q` in fm⁻¹ (spectator-1 convention).
//        - After fixing p̂ = ẑ and φ_q = 0, the only non-trivial angular
//          variable in the scalar partial-wave projection is
//          x ≡ cos(q̂, q̂'). Azimuthal averages are absorbed implicitly.
//        - Momentum transfers:
//            Δp = p' − p,    Δq = q' − q
//            q₁ = Δq                       (transferred to spectator 1)
//            q₂ = Δp + Δq/2                (to pair-particle 2)
//            q₃ = Δp − Δq/2                (to pair-particle 3)
//          With p̂ = ẑ and the azimuthal average over the p'-direction
//          orthogonal to the x-integration, cross-terms Δp·Δq drop out at
//          the scalar (rank-0) level, leaving:
//            |Δq|² = q² + q'² − 2 q q' x       ← [E2002] eq. (A-2)
//            ⟨|Δp|²⟩ = p² + p'²
//            ⟨q₂²⟩ = ⟨q₃²⟩ = ⟨|Δp|²⟩ + |Δq|²/4
//            ⟨q₂·q₃⟩ = ⟨|Δp|²⟩ − |Δq|²/4
//
//      Unit conventions (see §Units of formula_reference.md):
//        - momenta (p, q, p', q', m_π, Λ): fm⁻¹.
//        - fπ in fm (converted from MeV via fπ/ħc).
//        - Λ_χ in fm (converted from MeV via Λ_χ/ħc).
//        - c_E dimensionless; c_D dimensionless; c_1, c_3 in fm (converted
//          from GeV⁻¹ via ×ħc/1000 = 0.19733).
//
// [CN] 手征 N2LO 三体核力分波分解的标量动量空间核函数。
//      这些辅助函数将纯动量空间部分（π 传播子、LEC 前因子、正规化因子）
//      分离出来，使得外层分波机制（重耦合系数、自旋/同位旋矩阵元、x 积分）
