#include "cache_io_w1.h"
#include "cache_manifest.h"
#include "cache_schema.h"
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

#include <atomic>
#include <cstdio>      // std::rename, std::remove
#include <cstring>
#include <iostream>
#include <string>
#include <unistd.h>    // getpid

namespace tictac::cache {

// (Local helpers identical to those in cache_io_p123.cpp; we duplicate them
// here intentionally to keep IO modules self-contained — TUs don't share.)
static void write_str_attr(hid_t loc, const char* name, const std::string& v) {
    hid_t dt = H5Tcopy(H5T_C_S1);
    H5Tset_size(dt, v.size() + 1); H5Tset_strpad(dt, H5T_STR_NULLTERM);
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate2(loc, name, dt, sp, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, dt, v.c_str());
    H5Aclose(at); H5Sclose(sp); H5Tclose(dt);
}
static void write_int_attr(hid_t loc, const char* name, int v) {
    hid_t sp = H5Screate(H5S_SCALAR);
    hid_t at = H5Acreate2(loc, name, H5T_NATIVE_INT, sp, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(at, H5T_NATIVE_INT, &v); H5Aclose(at); H5Sclose(sp);
}
static bool read_str_attr(hid_t loc, const char* name, std::string& out) {
    if (H5Aexists(loc, name) <= 0) return false;
    hid_t at = H5Aopen(loc, name, H5P_DEFAULT);
    hid_t dt = H5Aget_type(at);
    size_t sz = H5Tget_size(dt);
    std::string buf(sz, '\0');
    H5Aread(at, dt, buf.data());
    if (!buf.empty() && buf.back() == '\0') buf.pop_back();
    out = buf; H5Tclose(dt); H5Aclose(at);
    return true;
}
static bool read_int_attr(hid_t loc, const char* name, int& out) {
    if (H5Aexists(loc, name) <= 0) return false;
    hid_t at = H5Aopen(loc, name, H5P_DEFAULT);
    H5Aread(at, H5T_NATIVE_INT, &out); H5Aclose(at);
    return true;
}

void write_w1_h5(const std::string& path,
                 const W1Key&       key,
                 const W1Block&     in,
                 const std::string& writer_version)
{
    // [EN] Atomic publication (Phase E crash safety). The block is written to a
    // unique temp file in the SAME directory (hence the same filesystem, where
    // POSIX rename() is atomic) and only then renamed onto the final path. A
    // process killed mid-write therefore never leaves a partial block visible
    // at the final path; the worst case is an orphaned temp file, which a later
    // reader never consults. Two workers building the same block each get a
    // distinct temp name (pid + per-process counter); the second rename simply
    // overwrites the first with another valid file. read_w1_h5 still re-checks
    // the key hash, so a corrupt file can never be silently accepted.
    // / [CN] 原子发布（Phase E 崩溃安全）：先写同目录唯一临时文件，再 rename 到最终路径。
    static std::atomic<unsigned long long> seq{0};
    std::string tmp = path + ".tmp." + std::to_string(::getpid()) + "."
                    + std::to_string(seq.fetch_add(1));

    hid_t file = H5Fcreate(tmp.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) { std::cerr << "write_w1_h5: failed (tmp create) " << tmp << "\n"; return; }

    // 3NF audit B6 (2026-06-21): changed from H5T_NATIVE_FLOAT to DOUBLE.
    hsize_t dims4[4] = { (hsize_t)in.Nq, (hsize_t)in.Nq, (hsize_t)in.Np, (hsize_t)in.Np };
    H5LTmake_dataset(file, "/data", 4, dims4, H5T_NATIVE_DOUBLE, in.data.data());

    write_int_attr(file, "Nq", in.Nq);
    write_int_attr(file, "Np", in.Np);
    write_int_attr(file, "a_r", in.a_r);
    write_int_attr(file, "a_c", in.a_c);

    write_str_attr(file, "tictac_cache_kind", "w1");
    write_int_attr(file, "tictac_schema_version", key.schema_version);
    write_str_attr(file, "tictac_key_hash_full", hash_full(key));
    write_str_attr(file, "tictac_key_json", canonical_json(key));
    write_str_attr(file, "tictac_created_utc", now_utc_iso8601());
    write_str_attr(file, "tictac_writer_version", writer_version);

    H5Fclose(file);

    // Atomically promote the complete temp file to the final path.
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        // Raced with another writer, or other failure: remove our orphan and
        // rely on whichever writer did succeed. The final path is still valid.
        std::remove(tmp.c_str());
    }
}

bool read_w1_h5(const std::string& path,
                const W1Key&       expected_key,
                W1Block*           out,
                std::string*       miss_reason)
{
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) { if (miss_reason) *miss_reason = "not_found"; return false; }

    std::string kind; int sv = 0; std::string khash;
    if (!read_str_attr(file, "tictac_cache_kind", kind) || kind != "w1") {
        if (miss_reason) *miss_reason = "kind_mismatch";
        H5Fclose(file); return false;
    }
    if (!read_int_attr(file, "tictac_schema_version", sv) || sv != expected_key.schema_version) {
        if (miss_reason) *miss_reason = "schema_mismatch";
        H5Fclose(file); return false;
    }
    if (!read_str_attr(file, "tictac_key_hash_full", khash) || khash != hash_full(expected_key)) {
        if (miss_reason) *miss_reason = "key_mismatch";
        H5Fclose(file); return false;
    }

    int Nq=0, Np=0, ar=0, ac=0;
    read_int_attr(file, "Nq", Nq); read_int_attr(file, "Np", Np);
    read_int_attr(file, "a_r", ar); read_int_attr(file, "a_c", ac);
    out->Nq = Nq; out->Np = Np; out->a_r = ar; out->a_c = ac;
    size_t total = (size_t)Nq * Nq * Np * Np;
    // 3NF audit B6: storage is now double.
    out->data.assign(total, 0.0);
    if (H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, out->data.data()) < 0) {
        if (miss_reason) *miss_reason = "corrupt";
        H5Fclose(file); return false;
    }
    H5Fclose(file);
    return true;
}

}  // namespace tictac::cache
