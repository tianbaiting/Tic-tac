#ifndef TICTAC_CORE_PERMUTATION_OPERATOR_VIEW_H
#define TICTAC_CORE_PERMUTATION_OPERATOR_VIEW_H

#include <cstddef>

namespace tictac::core {

// [EN] Read-only, non-owning view over the finite-dimensional cyclic
// permutation operator
//
//     P = P_123 + P_132           (sum of cyclic particle permutations;
//                                  P_123 = P_132 for antisymmetric pair states)
//
// stored as a CSC sparse matrix on the packet lattice of one conserved J^pi
// block (docs/three_nf_equation_contract.md §4.1 / §4.2). The CSC storage is:
//   values[k]      nonzero entry k
//   row_index[k]   row of nonzero k
//   col_ptr[j]     index into values/row_index of the first nonzero of column j
//                  (length dim+1; col_ptr[dim] == num_nonzeros)
// Storage is unchanged -- the view only aliases the existing arrays.
// / [CN] 对有限维循环置换算符 P = P_123 + P_132 的只读非拥有视图（按 J^pi 块的 CSC 稀疏存储）。
class PermutationOperatorView {
public:
	PermutationOperatorView() = default;
	PermutationOperatorView(const double* values,
	                        const int* row_index,
	                        const std::size_t* col_ptr,
	                        std::size_t dim)
		: values_(values), row_index_(row_index), col_ptr_(col_ptr), dim_(dim) {}

	std::size_t dimension() const { return dim_; }
	std::size_t num_nonzeros() const { return col_ptr_[dim_]; }

	// Half-open nonzero range [begin, end) of column j (for sparse matvec / the
	// PVC = P·(V·C) contraction in cpvc_kernel).
	struct ColumnNonzeros {
		std::size_t begin;
		std::size_t end;
		const double* values;
		const int* row_index;
	};
	ColumnNonzeros column(std::size_t j) const
	{
		return ColumnNonzeros{
			col_ptr_[j], col_ptr_[j + 1], values_, row_index_
		};
	}

	const double*     values()     const { return values_; }
	const int*        row_index()  const { return row_index_; }
	const std::size_t* col_ptr()   const { return col_ptr_; }

private:
	const double*      values_ = nullptr;
	const int*         row_index_ = nullptr;
	const std::size_t* col_ptr_ = nullptr;
	std::size_t        dim_ = 0;
};

} // namespace tictac::core

#endif // TICTAC_CORE_PERMUTATION_OPERATOR_VIEW_H
