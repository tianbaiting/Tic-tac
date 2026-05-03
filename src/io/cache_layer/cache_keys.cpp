#include "cache_keys.h"
#include "third_party/sha256.h"
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace tictac::cache {

double quantize_1e9(double x) {
    return std::floor(x * 1e9 + 0.5) / 1e9;
}

static std::string fmt_double_q(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.9f", quantize_1e9(x));
    return std::string(buf);
}

static std::string json_escape(const std::string& s) {
    std::string out; out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if ((unsigned char)c < 0x20) {
                    char hex[8]; std::snprintf(hex, sizeof(hex), "\\u%04x", c);
                    out += hex;
                } else out.push_back(c);
        }
    }
    out.push_back('"');
    return out;
}

bool P123Key::operator==(const P123Key& o) const {
    return schema_version == o.schema_version
        && Np_WP == o.Np_WP && Nq_WP == o.Nq_WP
        && J_2N_max == o.J_2N_max
        && two_J_3N == o.two_J_3N && P_3N == o.P_3N
        && Nphi == o.Nphi && Nx == o.Nx
        && tensor_force == o.tensor_force
        && isospin_breaking_1S0 == o.isospin_breaking_1S0
        && quantize_1e9(chebyshev_s) == quantize_1e9(o.chebyshev_s)
        && quantize_1e9(chebyshev_t) == quantize_1e9(o.chebyshev_t);
}

bool W1Key::operator==(const W1Key& o) const {
    return schema_version == o.schema_version
        && potential_model == o.potential_model
        && tnf_model == o.tnf_model
        && Np_WP == o.Np_WP && Nq_WP == o.Nq_WP
        && J_2N_max == o.J_2N_max && two_J_3N_max == o.two_J_3N_max
        && two_J_3N == o.two_J_3N && P_3N == o.P_3N
        && a_r == o.a_r && a_c == o.a_c
        && quantize_1e9(c_D) == quantize_1e9(o.c_D)
        && quantize_1e9(c_E) == quantize_1e9(o.c_E)
        && quantize_1e9(Lambda_3NF) == quantize_1e9(o.Lambda_3NF)
        && regulator_kind == o.regulator_kind
        && quantize_1e9(chebyshev_s) == quantize_1e9(o.chebyshev_s)
        && quantize_1e9(chebyshev_t) == quantize_1e9(o.chebyshev_t)
        && tensor_force == o.tensor_force
        && isospin_breaking_1S0 == o.isospin_breaking_1S0;
}

std::string canonical_json(const P123Key& k) {
    std::ostringstream os;
    // Keys sorted lexicographically.
    os << "{"
       << "\"J_2N_max\":" << k.J_2N_max
       << ",\"Nphi\":" << k.Nphi
       << ",\"Np_WP\":" << k.Np_WP
       << ",\"Nq_WP\":" << k.Nq_WP
       << ",\"Nx\":" << k.Nx
       << ",\"P_3N\":" << k.P_3N
       << ",\"chebyshev_s\":" << fmt_double_q(k.chebyshev_s)
       << ",\"chebyshev_t\":" << fmt_double_q(k.chebyshev_t)
       << ",\"isospin_breaking_1S0\":" << (k.isospin_breaking_1S0 ? "true" : "false")
       << ",\"schema_version\":" << k.schema_version
       << ",\"tensor_force\":" << (k.tensor_force ? "true" : "false")
       << ",\"two_J_3N\":" << k.two_J_3N
       << "}";
    return os.str();
}

std::string canonical_json(const W1Key& k) {
    std::ostringstream os;
    os << "{"
       << "\"J_2N_max\":" << k.J_2N_max
       << ",\"Lambda_3NF\":" << fmt_double_q(k.Lambda_3NF)
       << ",\"Np_WP\":" << k.Np_WP
       << ",\"Nq_WP\":" << k.Nq_WP
       << ",\"P_3N\":" << k.P_3N
       << ",\"a_c\":" << k.a_c
       << ",\"a_r\":" << k.a_r
       << ",\"c_D\":" << fmt_double_q(k.c_D)
       << ",\"c_E\":" << fmt_double_q(k.c_E)
       << ",\"chebyshev_s\":" << fmt_double_q(k.chebyshev_s)
       << ",\"chebyshev_t\":" << fmt_double_q(k.chebyshev_t)
       << ",\"isospin_breaking_1S0\":" << (k.isospin_breaking_1S0 ? "true" : "false")
       << ",\"potential_model\":" << json_escape(k.potential_model)
       << ",\"regulator_kind\":" << json_escape(k.regulator_kind)
       << ",\"schema_version\":" << k.schema_version
       << ",\"tensor_force\":" << (k.tensor_force ? "true" : "false")
       << ",\"tnf_model\":" << json_escape(k.tnf_model)
       << ",\"two_J_3N\":" << k.two_J_3N
       << ",\"two_J_3N_max\":" << k.two_J_3N_max
       << "}";
    return os.str();
}

std::string hash_full(const P123Key& k) { return sha256::hex_digest(canonical_json(k)); }
std::string hash_full(const W1Key&   k) { return sha256::hex_digest(canonical_json(k)); }

std::string hash_short(const std::string& full_hex) {
    return full_hex.substr(0, 6);
}

static std::string parity_str(int p) { return p > 0 ? "+1" : "-1"; }

std::string filename_prefix(const P123Key& k) {
    std::ostringstream os;
    os << "Np" << k.Np_WP << "_Nq" << k.Nq_WP
       << "_J2max" << k.J_2N_max
       << "_JP" << k.two_J_3N << parity_str(k.P_3N);
    return os.str();
}

std::string filename_prefix(const W1Key& k) {
    std::ostringstream os;
    os << k.potential_model
       << "__cD" << fmt_double_q(k.c_D)
       << "_cE" << fmt_double_q(k.c_E)
       << "_L"  << fmt_double_q(k.Lambda_3NF)
       << "__Np" << k.Np_WP << "_Nq" << k.Nq_WP
       << "_J2max" << k.J_2N_max
       << "_J3Nmax" << k.two_J_3N_max
       << "_JP" << k.two_J_3N << parity_str(k.P_3N)
       << "_a" << k.a_r << "_" << k.a_c;
    return os.str();
}

}  // namespace tictac::cache
