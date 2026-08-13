#ifndef TICTAC_CORE_PACKET_GRID_VIEW_H
#define TICTAC_CORE_PACKET_GRID_VIEW_H

#include <cstddef>

#include "type_defs.h"

namespace tictac::core {

// [EN] Read-only, non-owning view over one Jacobi-momentum wave-packet axis
// (p = pair-relative or q = spectator). A packet axis is a sequence of bins
// `[k_i, k_{i+1})` each carrying midpoint/average momentum quadrature nodes.
//
// The boundary array has N_WP+1 entries; the per-bin node arrays have
// N_WP * N_per_WP entries laid out bin-major. Storage is unchanged -- the view
// only aliases the existing fwp_statespace arrays.
// / [CN] 对单个 Jacobi 动量波包轴（p 或 q）的只读非拥有视图。
// 边界数组 N_WP+1 个；逐 bin 节点数组 N_WP*N_per_WP 个，按 bin 为主序。存储不变。
class PacketAxisView {
public:
	PacketAxisView() = default;
	PacketAxisView(int n_wp, int n_per_wp,
	               const double* boundaries,
	               const double* nodes, const double* weights,
	               const double* fvalues, const double* norms)
		: n_wp_(n_wp), n_per_wp_(n_per_wp),
		  boundaries_(boundaries), nodes_(nodes), weights_(weights),
		  fvalues_(fvalues), norms_(norms) {}

	int num_packets() const { return n_wp_; }
	int nodes_per_packet() const { return n_per_wp_; }

	// Lower/upper momentum boundary of packet i [fm^-1].
	double lower(std::size_t i) const { return boundaries_[i]; }
	double upper(std::size_t i) const { return boundaries_[i + 1]; }
	double width(std::size_t i) const { return upper(i) - lower(i); }

	// j-th quadrature node/weight inside packet i.
	double node(std::size_t i, std::size_t j) const { return nodes_[i * n_per_wp_ + j]; }
	double weight(std::size_t i, std::size_t j) const { return weights_[i * n_per_wp_ + j]; }
	double fvalue(std::size_t i, std::size_t j) const { return fvalues_[i * n_per_wp_ + j]; }
	double norm(std::size_t i, std::size_t j) const { return norms_[i * n_per_wp_ + j]; }

	const double* boundaries() const { return boundaries_; }
	const double* nodes() const { return nodes_; }

private:
	int n_wp_ = 0;
	int n_per_wp_ = 1;
	const double* boundaries_ = nullptr;
	const double* nodes_ = nullptr;
	const double* weights_ = nullptr;
	const double* fvalues_ = nullptr;
	const double* norms_ = nullptr;
};

// [EN] Read-only, non-owning view over the full two-axis (p, q) free wave-packet
// mesh `fwp_statespace`. The discretized operator lattice is the tensor product
// of the p and q axes; later operators are projected onto this p-q lattice.
// / [CN] 对完整 (p,q) 双轴自由波包网格 fwp_statespace 的只读非拥有视图。
class PacketGridView {
public:
	PacketGridView() = default;
	explicit PacketGridView(const fwp_statespace& fwp)
		: fwp_(&fwp),
		  p_(fwp.Np_WP, fwp.Np_per_WP, fwp.p_WP_array, fwp.p_array,
		     fwp.wp_array, fwp.fp_array, fwp.norm_p_array),
		  q_(fwp.Nq_WP, fwp.Nq_per_WP, fwp.q_WP_array, fwp.q_array,
		     fwp.wq_array, fwp.fq_array, fwp.norm_q_array) {}

	const PacketAxisView& p() const { return p_; }
	const PacketAxisView& q() const { return q_; }

	int Np_WP() const { return fwp_->Np_WP; }
	int Nq_WP() const { return fwp_->Nq_WP; }

	// Dense packet-space dimension = Nalpha * Nq * Np (set once the basis is known).
	std::size_t dense_dim(std::size_t nalpha) const
	{ return nalpha * static_cast<std::size_t>(Nq_WP()) * static_cast<std::size_t>(Np_WP()); }

	const fwp_statespace& storage() const { return *fwp_; }

private:
	const fwp_statespace* fwp_ = nullptr;
	PacketAxisView p_;
	PacketAxisView q_;
};

} // namespace tictac::core

#endif // TICTAC_CORE_PACKET_GRID_VIEW_H
