#include "coupled_channel_transform.h"

#include <algorithm>
#include <vector>

void restructure_coupled_matrix_blocks(double* matrix, std::size_t Np_WP) {
    const std::size_t block_size = Np_WP * Np_WP;
    std::vector<double> blocks(4U * block_size);

    for (std::size_t col_block = 0; col_block < 2; ++col_block) {
        for (std::size_t row_block = 0; row_block < 2; ++row_block) {
            const std::size_t block = col_block * 2U + row_block;
            for (std::size_t col = 0; col < Np_WP; ++col) {
                for (std::size_t row = 0; row < Np_WP; ++row) {
                    const std::size_t source_col = col + Np_WP * col_block;
                    const std::size_t source_row = row + Np_WP * row_block;
                    blocks[block * block_size + col * Np_WP + row]
                        = matrix[source_col * (2U * Np_WP) + source_row];
                }
            }
        }
    }

    std::copy(blocks.begin(), blocks.end(), matrix);
}
