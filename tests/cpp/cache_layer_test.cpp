#include "cache_keys.h"
#include "third_party/sha256.h"
#include <cassert>
#include <iostream>
#include <string>

using tictac::cache::P123Key;
using tictac::cache::W1Key;

static int failures = 0;

#define EXPECT(cond) do { if (!(cond)) { \
    std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << " " << #cond << "\n"; \
    ++failures; } } while(0)

#define EXPECT_EQ(a, b) do { auto _a = (a); auto _b = (b); if (!(_a == _b)) { \
    std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << " expected " \
              << #a << " == " << #b << " but got " << _a << " vs " << _b << "\n"; \
    ++failures; } } while(0)

static P123Key make_p123_key() {
    P123Key k{};
    k.schema_version = 1;
    k.Np_WP = 30; k.Nq_WP = 30;
    k.J_2N_max = 2; k.two_J_3N = 1; k.P_3N = 1;
    k.Nphi = 24; k.Nx = 24;
    k.tensor_force = true;
    k.isospin_breaking_1S0 = false;
    k.chebyshev_s = 1.5; k.chebyshev_t = 1.0;
    return k;
}

static W1Key make_w1_key() {
    W1Key k{};
    k.schema_version = 1;
    k.potential_model = "N2LOopt";
    k.tnf_model = "chiral_N2LO";
    k.Np_WP = 30; k.Nq_WP = 30;
    k.J_2N_max = 2; k.two_J_3N_max = 1;
    k.two_J_3N = 1; k.P_3N = 1;
    k.c_D = -0.20; k.c_E = -0.205;
    k.Lambda_3NF = 500.0;
    k.regulator_kind = "gaussian";
    return k;
}

void test_sha256_known_vectors() {
    auto h = tictac::cache::sha256::hex_digest("");
    EXPECT_EQ(h, std::string("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
    h = tictac::cache::sha256::hex_digest("abc");
    EXPECT_EQ(h, std::string("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"));
}

void test_p123_canonical_json_keys_sorted() {
    auto k = make_p123_key();
    auto j = tictac::cache::canonical_json(k);
    auto pos_J2 = j.find("\"J_2N_max\"");
    auto pos_Np = j.find("\"Np_WP\"");
    auto pos_Nq = j.find("\"Nq_WP\"");
    auto pos_P3 = j.find("\"P_3N\"");
    EXPECT(pos_J2 < pos_Np);
    EXPECT(pos_Np < pos_Nq);
    EXPECT(pos_Nq < pos_P3);
}

void test_p123_hash_stable() {
    auto k = make_p123_key();
    auto h1 = tictac::cache::hash_full(k);
    auto h2 = tictac::cache::hash_full(k);
    EXPECT_EQ(h1, h2);
    EXPECT_EQ(h1.size(), (size_t)64);
}

void test_p123_hash_changes_on_field_change() {
    auto k1 = make_p123_key();
    auto k2 = make_p123_key();
    k2.Np_WP = 50;
    EXPECT(tictac::cache::hash_full(k1) != tictac::cache::hash_full(k2));
}

void test_w1_double_quantization() {
    auto k1 = make_w1_key();
    auto k2 = make_w1_key();
    k2.c_D = -0.20000000001;  // below 1e-9 quantum
    EXPECT_EQ(tictac::cache::hash_full(k1), tictac::cache::hash_full(k2));
    EXPECT(k1 == k2);
}

void test_w1_double_above_quantum_changes_hash() {
    auto k1 = make_w1_key();
    auto k2 = make_w1_key();
    k2.c_D = -0.2001;  // well above 1e-9 quantum
    EXPECT(tictac::cache::hash_full(k1) != tictac::cache::hash_full(k2));
}

void test_filename_prefix_p123() {
    auto k = make_p123_key();
    auto p = tictac::cache::filename_prefix(k);
    EXPECT_EQ(p, std::string("Np30_Nq30_J2max2_JP1+1"));
}

int main() {
    test_sha256_known_vectors();
    test_p123_canonical_json_keys_sorted();
    test_p123_hash_stable();
    test_p123_hash_changes_on_field_change();
    test_w1_double_quantization();
    test_w1_double_above_quantum_changes_hash();
    test_filename_prefix_p123();

    if (failures > 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "all key tests passed\n";
    return 0;
}
