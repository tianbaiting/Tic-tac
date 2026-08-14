// W1WorkPlan implementation: deterministic evaluate/transpose classification.
//
// Mirrors the block enumeration and Hermitian-triangle rule already used inside
// W1_PW_cache::build, lifted into a pure, cache-state-independent, testable
// description. No physics or value change -- this only names the work.

#include "w1_work_plan.h"

#include <stdexcept>

namespace tictac::interactions {

pw_3N_statespace make_channel_view(const pw_3N_statespace& pw, int chn)
{
	// Same slice as solver_pipeline / w1_worker: non-owning view over the
	// channel's contiguous alpha range.  No allocation; the result aliases pw.
	pw_3N_statespace v = {};
	if (chn < 0 || chn >= pw.N_chn_3N) {
		throw std::runtime_error("make_channel_view: channel index out of range");
	}
	const int lo = pw.chn_3N_idx_array[chn];
	const int hi = pw.chn_3N_idx_array[chn + 1];
	v.Nalpha = hi - lo;
	v.J_2N_max = pw.J_2N_max;
	v.L_2N_array     = &pw.L_2N_array[lo];
	v.S_2N_array     = &pw.S_2N_array[lo];
	v.J_2N_array     = &pw.J_2N_array[lo];
	v.T_2N_array     = &pw.T_2N_array[lo];
	v.L_1N_array     = &pw.L_1N_array[lo];
	v.two_J_1N_array = &pw.two_J_1N_array[lo];
	v.two_T_3N_array = &pw.two_T_3N_array[lo];
	v.two_J_3N_array = &pw.two_J_3N_array[lo];
	v.P_3N_array     = &pw.P_3N_array[lo];
	return v;
}

W1WorkPlan::W1WorkPlan(const pw_3N_statespace& pw, bool hermitian)
{
	m_Nalpha = pw.Nalpha;
	m_block_index.assign(static_cast<std::size_t>(m_Nalpha) * m_Nalpha, -1);
	m_units.clear();
	m_units.reserve(static_cast<std::size_t>(m_Nalpha) * m_Nalpha);

	// Enumerate allowed (a_r,a_c) pairs in the SAME order as
	// W1_PW_cache::build: outer a_r, inner a_c, accept iff the conserved
	// quantum numbers (2J_3N, 2T_3N, parity) match. This fixes the
	// sector-local block id (== position in m_units) deterministically.
	for (int a_r = 0; a_r < m_Nalpha; ++a_r) {
		for (int a_c = 0; a_c < m_Nalpha; ++a_c) {
			if (pw.two_J_3N_array[a_r] != pw.two_J_3N_array[a_c]) continue;
			if (pw.two_T_3N_array[a_r] != pw.two_T_3N_array[a_c]) continue;
			if (pw.P_3N_array[a_r]     != pw.P_3N_array[a_c])     continue;

			W1WorkUnit u;
			u.id  = static_cast<int>(m_units.size());
			u.a_r = a_r;
			u.a_c = a_c;
			m_block_index[static_cast<std::size_t>(a_r) * m_Nalpha + a_c] = u.id;
			m_units.push_back(u);
		}
	}

	// Classify evaluate / transpose_fill using the absolute canonical rule:
	// the member of a Hermitian pair with the smaller sector-local block id is
	// `evaluate`; its reverse is `transpose_fill`.  Diagonals and pairs with
	// no reverse block are `evaluate`.  This is the cache-state-independent
	// form of W1_PW_cache::build's `blk < reverse` test; the two agree because
	// pre-building every canonical block makes every reverse block a
	// transpose-fill-at-solve (see header doc).
	if (hermitian) {
		for (auto& u : m_units) {
			const int rev = block_index(u.a_c, u.a_r);
			if (rev < 0 || rev == u.id) {
				u.role = W1UnitRole::evaluate;
			} else if (u.id < rev) {
				u.role = W1UnitRole::evaluate;            // canonical orientation
			} else {
				u.role = W1UnitRole::transpose_fill;
				u.conj_a_r = u.a_c;
				u.conj_a_c = u.a_r;
			}
		}
	} else {
		// Non-Hermitian model: every block must be integrated directly,
		// matching W1_PW_cache::build's fallback (evaluation_blocks = missing).
		for (auto& u : m_units) u.role = W1UnitRole::evaluate;
	}

	m_num_eval = 0;
	for (const auto& u : m_units)
		if (u.role == W1UnitRole::evaluate) ++m_num_eval;
}

std::vector<const W1WorkUnit*> W1WorkPlan::evaluate_units() const
{
	std::vector<const W1WorkUnit*> out;
	out.reserve(m_num_eval);
	for (const auto& u : m_units)
		if (u.role == W1UnitRole::evaluate) out.push_back(&u);
	return out;
}

} // namespace tictac::interactions
