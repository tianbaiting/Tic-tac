#pragma once
#include "cache_keys.h"
#include <string>
#include <vector>

namespace tictac::cache {

struct W1Block {
    int                 Nq;        // q-WP count for this block
    int                 Np;        // p-WP count for this block
    int                 a_r;       // alpha row index (informational)
    int                 a_c;       // alpha col index (informational)
    // 3NF audit B6 (2026-06-21): storage changed from float to double.
    std::vector<double> data;      // size = Nq*Nq*Np*Np
};

bool read_w1_h5 (const std::string& path,
                 const W1Key&       expected_key,
                 W1Block*           out,
                 std::string*       miss_reason);

void write_w1_h5(const std::string& path,
                 const W1Key&       key,
                 const W1Block&     in,
                 const std::string& writer_version);

}  // namespace tictac::cache
