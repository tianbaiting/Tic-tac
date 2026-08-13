#ifndef TICTAC_CORE_PARTIAL_WAVE_BASIS_VIEW_H
#define TICTAC_CORE_PARTIAL_WAVE_BASIS_VIEW_H

#include <cstddef>

#include "three_body_channel.h"
#include "type_defs.h"

namespace tictac::core {

// [EN] Read-only, non-owning semantic view over the legacy channel-by-channel
// SoA storage `pw_3N_statespace`. It exposes the partial-wave basis as the
// mathematical object it is -- a sequence of three-body channels alpha grouped
// into independently-solvable conserved (2J, 2T, parity) blocks -- instead of a
// bag of raw integer arrays.
//
// STORAGE IS UNCHANGED: the view only holds pointers into the existing
// `pw_3N_statespace` arrays, so it adds zero copy and zero numerical impact.
// It is the public semantic interface; the SoA layout remains the
// performance-critical storage representation.
// / [CN] 对旧式按通道 SoA 存储 pw_3N_statespace 的只读非拥有语义视图。
// 把分波基暴露为其数学本质（按守恒 (2J,2T,parity) 分块的通道序列），而非一袋裸数组。
// 存储不变：仅持有指向现有数组的指针，零拷贝、零数值影响。
class PartialWaveBasisView {
public:
	PartialWaveBasisView() = default;
	explicit PartialWaveBasisView(const pw_3N_statespace& pw) : pw_(&pw) {}

	// Total number of partial waves (the alpha dimension).
	std::size_t num_channels() const { return static_cast<std::size_t>(pw_->Nalpha); }

	// Number of independently-solvable conserved (2J,2T,parity) blocks.
	std::size_t num_blocks() const { return static_cast<std::size_t>(pw_->N_chn_3N); }

	// One three-body partial wave alpha.
	ThreeBodyChannel channel(std::size_t alpha) const
	{
		const int a = static_cast<int>(alpha);
		return ThreeBodyChannel{
			pw_->L_2N_array[a],
			pw_->S_2N_array[a],
			pw_->J_2N_array[a],
			pw_->T_2N_array[a],
			pw_->L_1N_array[a],
			pw_->two_J_1N_array[a],
			pw_->two_J_3N_array[a],
			pw_->two_T_3N_array[a],
			pw_->P_3N_array[a],
		};
	}

	// Half-open alpha range [begin, end) of conserved block `block`.
	struct BlockRange { std::size_t begin; std::size_t end; };
	BlockRange block_range(std::size_t block) const
	{
		return BlockRange{
			static_cast<std::size_t>(pw_->chn_3N_idx_array[block]),
			static_cast<std::size_t>(pw_->chn_3N_idx_array[block + 1]),
		};
	}
	std::size_t block_size(std::size_t block) const
	{
		const auto r = block_range(block);
		return r.end - r.begin;
	}

	// Conserved quantum numbers are constant across a block; read from its
	// first channel.
	int block_two_J(std::size_t block) const
	{ return pw_->two_J_3N_array[pw_->chn_3N_idx_array[block]]; }
	int block_two_T(std::size_t block) const
	{ return pw_->two_T_3N_array[pw_->chn_3N_idx_array[block]]; }
	int block_parity(std::size_t block) const
	{ return pw_->P_3N_array[pw_->chn_3N_idx_array[block]]; }

	int J_2N_max() const { return pw_->J_2N_max; }

	const pw_3N_statespace& storage() const { return *pw_; }

private:
	const pw_3N_statespace* pw_ = nullptr;
};

} // namespace tictac::core

#endif // TICTAC_CORE_PARTIAL_WAVE_BASIS_VIEW_H
