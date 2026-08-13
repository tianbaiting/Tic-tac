#include "cache_keys.h"
#include "cache_manifest.h"
#include "cache_io_p123.h"
#include "cache_io_w1.h"
#include "cache_layer.h"
#include "third_party/sha256.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

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
    k.schema_version = 2;
    k.Np_WP = 30; k.Nq_WP = 30;
    k.J_2N_max = 2; k.two_J_3N = 1; k.P_3N = 1;
    k.Nphi = 24; k.Nx = 24;
    k.tensor_force = true;
    k.isospin_breaking_1S0 = false;
    k.p_grid_hash = "p-grid-a"; k.q_grid_hash = "q-grid-a";
    return k;
}

static W1Key make_w1_key() {
    W1Key k{};
    k.schema_version = 6;
    k.potential_model = "N2LOopt";
    k.tnf_model = "chiral_N2LO";
    k.Np_WP = 30; k.Nq_WP = 30;
    k.J_2N_max = 2; k.two_J_3N_max = 1;
    k.two_J_3N = 1; k.P_3N = 1;
    k.c_D = -0.20; k.c_E = -0.205;
    k.Lambda_3NF = 500.0;
    // B5: NEW fields — c_1, c_3, c_4 must be initialised (defaults break
    // aggregator-style usage). N2LOopt values from constants.h.
    k.c_1 = -0.86; k.c_3 = -3.40; k.c_4 = 0.0;
    k.g_A = 1.289; k.f_pi_MeV = 92.2; k.m_pi_MeV = 138.039;
    k.Lambda_chi_MeV = 700.0; k.hbarc_MeV_fm = 197.327;
    k.regulator_kind = "gaussian";
    k.a_r = 0; k.a_c = 0;
    k.chebyshev_s = 1.5; k.chebyshev_t = 1.0;
    k.tensor_force = true;
    k.isospin_breaking_1S0 = false;
    // B5: NEW grid hashes — empty disables the check (back-compat).
    k.p_grid_hash = ""; k.q_grid_hash = "";
    k.Np_per_WP_W1 = 1; k.Nq_per_WP_W1 = 1;
    k.Nangle_3NF = 0;
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
    auto pos_p_grid = j.find("\"p_grid_hash\"");
    auto pos_q_grid = j.find("\"q_grid_hash\"");
    auto pos_schema = j.find("\"schema_version\"");
    EXPECT(pos_J2 < pos_Np);
    EXPECT(pos_Np < pos_Nq);
    EXPECT(pos_Nq < pos_P3);
    EXPECT(pos_P3 < pos_p_grid);
    EXPECT(pos_p_grid < pos_q_grid);
    EXPECT(pos_q_grid < pos_schema);
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

void test_p123_grid_hashes_change_identity() {
    auto base = make_p123_key();
    auto changed_p = base;
    changed_p.p_grid_hash = "p-grid-b";
    EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(changed_p));
    EXPECT(!(base == changed_p));

    auto changed_q = base;
    changed_q.q_grid_hash = "q-grid-b";
    EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(changed_q));
    EXPECT(!(base == changed_q));
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

void test_w1_a_indices_change_hash() {
    auto k1 = make_w1_key();
    auto k2 = make_w1_key();
    k2.a_r = 1;  // different alpha row
    EXPECT(tictac::cache::hash_full(k1) != tictac::cache::hash_full(k2));
}

void test_w1_grid_and_channel_fields_change_hash() {
    auto base = make_w1_key();
    {
        auto k = make_w1_key(); k.chebyshev_s = 2.0;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.chebyshev_t = 0.5;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.tensor_force = false;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.isospin_breaking_1S0 = true;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.Np_per_WP_W1 = 4;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.Nq_per_WP_W1 = 4;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.Nangle_3NF = 4;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
    {
        auto k = make_w1_key(); k.g_A = 1.276;
        EXPECT(tictac::cache::hash_full(base) != tictac::cache::hash_full(k));
    }
}

void test_filename_prefix_p123() {
    auto k = make_p123_key();
    auto p = tictac::cache::filename_prefix(k);
    EXPECT_EQ(p, std::string("Np30_Nq30_J2max2_JP1+1"));
}

static std::string make_tmpdir() {
    char tmpl[] = "/tmp/tictac_cache_test_XXXXXX";
    char* d = mkdtemp(tmpl);
    if (!d) { std::cerr << "mkdtemp failed\n"; std::exit(2); }
    return std::string(d);
}

static void rmrf(const std::string& path) {
    std::string cmd = "rm -rf '" + path + "'";
    int rc = std::system(cmd.c_str());
    (void)rc;
}

void test_manifest_roundtrip() {
    auto root = make_tmpdir();
    tictac::cache::Manifest m;
    tictac::cache::ManifestEntry e{};
    e.kind = "p123";
    e.filename = "p123/Np30_Nq30_J2max2_JP1+1__a3f9c2.h5";
    e.key_hash_full = std::string(64, 'a');
    e.key_json = "{\"Np_WP\":30}";
    e.size_bytes = 12345;
    e.schema_version = 1;
    e.created_utc = "2026-05-01T00:00:00Z";
    e.last_used_utc = "2026-05-01T00:00:00Z";
    e.legacy_imported = false;
    e.writer_version = "tictac-test";
    bool inserted = m.upsert(e);
    EXPECT(inserted);

    m.save(root);

    auto m2 = tictac::cache::Manifest::load(root);
    EXPECT_EQ(m2.entries().size(), (size_t)1);
    auto found = m2.find(e.key_hash_full);
    EXPECT(found != nullptr);
    if (found) {
        EXPECT_EQ(found->filename, e.filename);
        EXPECT_EQ(found->size_bytes, e.size_bytes);
    }
    rmrf(root);
}

void test_manifest_missing_returns_empty() {
    auto root = make_tmpdir();
    auto m = tictac::cache::Manifest::load(root);
    EXPECT_EQ(m.entries().size(), (size_t)0);
    rmrf(root);
}

void test_manifest_corrupt_returns_empty() {
    auto root = make_tmpdir();
    {
        std::ofstream out(root + "/manifest.json");
        out << "{ this is not valid json";
    }
    auto m = tictac::cache::Manifest::load(root);
    EXPECT_EQ(m.entries().size(), (size_t)0);
    rmrf(root);
}

void test_manifest_upsert_replaces() {
    tictac::cache::Manifest m;
    tictac::cache::ManifestEntry e{};
    e.kind = "p123"; e.key_hash_full = std::string(64, 'b');
    e.size_bytes = 100; e.filename = "p123/x.h5";
    EXPECT(m.upsert(e));
    e.size_bytes = 200;
    EXPECT(!m.upsert(e));
    EXPECT_EQ(m.find(e.key_hash_full)->size_bytes, (int64_t)200);
}

void test_p123_io_roundtrip() {
    auto root = make_tmpdir();
    auto path = root + "/test_p123.h5";

    auto k = make_p123_key();
    tictac::cache::P123Arrays w{};
    w.nnz = 5;
    w.val_array = new double[5]{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    w.row_array = new int[5]   { 0, 1, 2, 3, 4 };
    w.col_array = new int[5]   { 4, 3, 2, 1, 0 };

    tictac::cache::write_p123_h5(path, k, w, "tictac-test");

    tictac::cache::P123Arrays r{};
    std::string miss;
    bool ok = tictac::cache::read_p123_h5(path, k, &r, &miss);
    EXPECT(ok);
    EXPECT_EQ(r.nnz, (size_t)5);
    if (r.nnz == 5) {
        for (int i = 0; i < 5; ++i) {
            EXPECT_EQ(r.val_array[i], w.val_array[i]);
            EXPECT_EQ(r.row_array[i], w.row_array[i]);
            EXPECT_EQ(r.col_array[i], w.col_array[i]);
        }
    }
    delete[] w.val_array; delete[] w.row_array; delete[] w.col_array;
    delete[] r.val_array; delete[] r.row_array; delete[] r.col_array;
    rmrf(root);
}

void test_p123_io_key_mismatch_rejected() {
    auto root = make_tmpdir();
    auto path = root + "/test_p123.h5";

    auto k1 = make_p123_key();
    tictac::cache::P123Arrays w{}; w.nnz = 1;
    w.val_array = new double[1]{1.0}; w.row_array = new int[1]{0}; w.col_array = new int[1]{0};
    tictac::cache::write_p123_h5(path, k1, w, "tictac-test");

    auto k2 = make_p123_key();
    k2.Np_WP = 50;  // different key

    tictac::cache::P123Arrays r{};
    std::string miss;
    bool ok = tictac::cache::read_p123_h5(path, k2, &r, &miss);
    EXPECT(!ok);
    EXPECT_EQ(miss, std::string("key_mismatch"));

    delete[] w.val_array; delete[] w.row_array; delete[] w.col_array;
    rmrf(root);
}

void test_w1_io_roundtrip() {
    auto root = make_tmpdir();
    auto path = root + "/test_w1.h5";

    auto k = make_w1_key();
    tictac::cache::W1Block b{};
    b.Nq = 2; b.Np = 3; b.a_r = 4; b.a_c = 5;
    b.data.assign(b.Nq * b.Nq * b.Np * b.Np, 0.0f);
    for (size_t i = 0; i < b.data.size(); ++i) b.data[i] = (float)i * 0.5f;

    tictac::cache::write_w1_h5(path, k, b, "tictac-test");

    tictac::cache::W1Block r{};
    std::string miss;
    bool ok = tictac::cache::read_w1_h5(path, k, &r, &miss);
    EXPECT(ok);
    EXPECT_EQ(r.Nq, b.Nq); EXPECT_EQ(r.Np, b.Np);
    EXPECT_EQ(r.a_r, b.a_r); EXPECT_EQ(r.a_c, b.a_c);
    EXPECT_EQ(r.data.size(), b.data.size());
    for (size_t i = 0; i < b.data.size(); ++i)
        EXPECT_EQ(r.data[i], b.data[i]);
    rmrf(root);
}

void test_cache_layer_p123_roundtrip() {
    auto root = make_tmpdir();
    tictac::cache::initialize(root);

    auto k = make_p123_key();
    tictac::cache::P123Arrays w{};
    w.nnz = 3;
    w.val_array = new double[3]{ 7.0, 8.0, 9.0 };
    w.row_array = new int[3]   { 1, 2, 3 };
    w.col_array = new int[3]   { 0, 0, 0 };

    tictac::cache::store_p123(k, w);

    tictac::cache::P123Arrays r{};
    auto res = tictac::cache::lookup_p123(k, &r);
    EXPECT(res.hit);
    EXPECT_EQ(r.nnz, (size_t)3);

    auto stats = tictac::cache::summary();
    EXPECT_EQ(stats.n_p123, (size_t)1);

    delete[] w.val_array; delete[] w.row_array; delete[] w.col_array;
    delete[] r.val_array; delete[] r.row_array; delete[] r.col_array;
    tictac::cache::shutdown();
    rmrf(root);
}

void test_cache_layer_lookup_miss_when_no_file() {
    auto root = make_tmpdir();
    tictac::cache::initialize(root);

    auto k = make_p123_key();
    tictac::cache::P123Arrays r{};
    auto res = tictac::cache::lookup_p123(k, &r);
    EXPECT(!res.hit);
    EXPECT_EQ(res.miss_reason, std::string("not_found"));

    tictac::cache::shutdown();
    rmrf(root);
}

void test_cache_layer_w1_lookup_miss_when_no_file() {
    auto root = make_tmpdir();
    tictac::cache::initialize(root);

    auto k = make_w1_key();
    tictac::cache::W1Block r{};
    auto res = tictac::cache::lookup_w1(k, &r);
    EXPECT(!res.hit);
    EXPECT_EQ(res.miss_reason, std::string("not_found"));

    tictac::cache::shutdown();
    rmrf(root);
}

void test_cache_layer_w1_roundtrip() {
    auto root = make_tmpdir();
    tictac::cache::initialize(root);

    auto k = make_w1_key();
    tictac::cache::W1Block b{};
    b.Nq = 2; b.Np = 2; b.a_r = 0; b.a_c = 0;
    b.data.assign(b.Nq*b.Nq*b.Np*b.Np, 0.0);
    for (size_t i = 0; i < b.data.size(); ++i) b.data[i] = static_cast<double>(i + 1);

    tictac::cache::store_w1(k, b);

    tictac::cache::W1Block r{};
    auto res = tictac::cache::lookup_w1(k, &r);
    EXPECT(res.hit);
    EXPECT_EQ(r.data.size(), b.data.size());
    for (size_t i = 0; i < b.data.size(); ++i)
        EXPECT_EQ(r.data[i], b.data[i]);

    tictac::cache::shutdown();
    rmrf(root);
}

int main() {
    test_sha256_known_vectors();
    test_p123_canonical_json_keys_sorted();
    test_p123_hash_stable();
    test_p123_hash_changes_on_field_change();
    test_p123_grid_hashes_change_identity();
    test_w1_double_quantization();
    test_w1_double_above_quantum_changes_hash();
    test_w1_a_indices_change_hash();
    test_w1_grid_and_channel_fields_change_hash();
    test_filename_prefix_p123();

    test_manifest_roundtrip();
    test_manifest_missing_returns_empty();
    test_manifest_corrupt_returns_empty();
    test_manifest_upsert_replaces();

    test_p123_io_roundtrip();
    test_p123_io_key_mismatch_rejected();

    test_w1_io_roundtrip();

    test_cache_layer_p123_roundtrip();
    test_cache_layer_lookup_miss_when_no_file();
    test_cache_layer_w1_lookup_miss_when_no_file();
    test_cache_layer_w1_roundtrip();

    if (failures > 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "all key tests passed\n";
    return 0;
}
