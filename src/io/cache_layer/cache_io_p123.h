#pragma once
#include "cache_keys.h"
#include <string>
#include <cstdint>

namespace tictac::cache {

struct P123Arrays {
    double*  val_array;   // owned by caller
    int*     row_array;
    int*     col_array;
    size_t   nnz;
};

// Reads the .h5 file at `path`. Validates kind/version/key match `key`.
// On success, fills `out` (allocates val/row/col with new[]) and returns true.
// On failure (mismatch, corrupt, missing), returns false; out is untouched.
// `miss_reason` is set to a short tag.
bool read_p123_h5(const std::string& path,
                  const P123Key&     expected_key,
                  P123Arrays*        out,
                  std::string*       miss_reason);

// Writes the file at `path` with sparse arrays + root attributes.
// Overwrites any existing file at path.
void write_p123_h5(const std::string& path,
                   const P123Key&     key,
                   const P123Arrays&  in,
                   const std::string& writer_version);

}  // namespace tictac::cache
