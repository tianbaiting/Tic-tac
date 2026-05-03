#pragma once
#include <string>
#include <cstdint>

namespace tictac::cache {

struct P123Key {
    int    schema_version;
    int    Np_WP, Nq_WP;
    int    J_2N_max;
    int    two_J_3N, P_3N;
    int    Nphi, Nx;
    bool   tensor_force;
    bool   isospin_breaking_1S0;
    double chebyshev_s, chebyshev_t;

    bool operator==(const P123Key& o) const;
};

struct W1Key {
    int    schema_version;
    std::string potential_model;
    std::string tnf_model;
    int    Np_WP, Nq_WP;
    int    J_2N_max, two_J_3N_max;
    int    two_J_3N, P_3N;
    int    a_r, a_c;          // alpha row/col indices for the block
    double c_D, c_E, Lambda_3NF;
    std::string regulator_kind;
    // Grid + channel parameters: change them, the WP bin midpoints change,
    // hence the cached W^(1) values change.
    double chebyshev_s, chebyshev_t;
    bool   tensor_force;
    bool   isospin_breaking_1S0;
    // Cell-quadrature order used when building the W^(1) bin matrix elements.
    // 1 = legacy midpoint (single-point Gauss), >=2 = cell-averaged via Gauss-Legendre
    // quadrature. Different orders produce different cached values, so they must
    // appear in the key.
    int    Np_per_WP_W1, Nq_per_WP_W1;

    bool operator==(const W1Key& o) const;
};

// Canonical JSON: keys sorted lexicographically; doubles formatted via
// snprintf("%.9f", floor(x * 1e9 + 0.5) / 1e9).
std::string canonical_json(const P123Key& k);
std::string canonical_json(const W1Key& k);

// Returns 64-char lowercase hex SHA-256 over canonical_json.
std::string hash_full(const P123Key& k);
std::string hash_full(const W1Key& k);

// Filename-safe 6-character prefix for use as filename suffix.
std::string hash_short(const std::string& full_hex);

// Human-readable filename prefix (lossless of int fields, lossy of doubles).
std::string filename_prefix(const P123Key& k);
std::string filename_prefix(const W1Key& k);

// Quantize a double to 1e-9 precision.
double quantize_1e9(double x);

}  // namespace tictac::cache
