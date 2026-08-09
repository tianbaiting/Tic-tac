#ifndef COUPLED_CHANNEL_TRANSFORM_H
#define COUPLED_CHANNEL_TRANSFORM_H

#include <cstddef>

// After a 2Np x 2Np coupled-channel matrix has been transposed, split it into
// four contiguous Np x Np blocks.  Each block remains row-major in its two
// packet indices.  The block ordering is
//
//     (lower,lower), (lower,upper), (upper,lower), (upper,upper)
//
// where the first label is the external row channel of the transposed matrix.
void restructure_coupled_matrix_blocks(double* matrix, std::size_t Np_WP);

inline constexpr std::size_t coupled_matrix_block_index(bool row_is_upper,
                                                        bool col_is_upper) {
    return 2U * static_cast<std::size_t>(row_is_upper)
         +      static_cast<std::size_t>(col_is_upper);
}

#endif
