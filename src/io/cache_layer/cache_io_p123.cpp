#include "cache_io_p123.h"
#include "cache_manifest.h"
#include "cache_schema.h"
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/hdf5_hl.h"

#include <cstring>
#include <iostream>

namespace tictac::cache {

static void write_str_attr(hid_t loc, const char* name, const std::string& v) {
    hid_t dtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(dtype, v.size() + 1);
    H5Tset_strpad(dtype, H5T_STR_NULLTERM);
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr  = H5Acreate2(loc, name, dtype, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, dtype, v.c_str());
    H5Aclose(attr); H5Sclose(space); H5Tclose(dtype);
}

static void write_int_attr(hid_t loc, const char* name, int v) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr  = H5Acreate2(loc, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &v);
    H5Aclose(attr); H5Sclose(space);
}

static bool read_str_attr(hid_t loc, const char* name, std::string& out) {
    if (H5Aexists(loc, name) <= 0) return false;
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    hid_t dtype = H5Aget_type(attr);
    size_t sz = H5Tget_size(dtype);
    std::string buf(sz, '\0');
    H5Aread(attr, dtype, buf.data());
    if (!buf.empty() && buf.back() == '\0') buf.pop_back();
    out = buf;
    H5Tclose(dtype); H5Aclose(attr);
    return true;
}

static bool read_int_attr(hid_t loc, const char* name, int& out) {
    if (H5Aexists(loc, name) <= 0) return false;
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &out);
    H5Aclose(attr);
    return true;
}

void write_p123_h5(const std::string& path,
                   const P123Key&     key,
                   const P123Arrays&  in,
                   const std::string& writer_version)
{
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "write_p123_h5: failed to create " << path << "\n";
        return;
    }

    hsize_t dims[1] = { (hsize_t)in.nnz };
    H5LTmake_dataset(file, "/val", 1, dims, H5T_NATIVE_DOUBLE, in.val_array);
    H5LTmake_dataset(file, "/row", 1, dims, H5T_NATIVE_INT, in.row_array);
    H5LTmake_dataset(file, "/col", 1, dims, H5T_NATIVE_INT, in.col_array);
    hsize_t one[1] = { 1 };
    uint64_t nnz_v = (uint64_t)in.nnz;
    H5LTmake_dataset(file, "/nnz", 1, one, H5T_NATIVE_UINT64, &nnz_v);

    write_str_attr(file, "tictac_cache_kind", "p123");
    write_int_attr(file, "tictac_schema_version", key.schema_version);
    write_str_attr(file, "tictac_key_hash_full", hash_full(key));
    write_str_attr(file, "tictac_key_json", canonical_json(key));
    write_str_attr(file, "tictac_created_utc", now_utc_iso8601());
    write_str_attr(file, "tictac_writer_version", writer_version);

    H5Fclose(file);
}

bool read_p123_h5(const std::string& path,
                  const P123Key&     expected_key,
                  P123Arrays*        out,
                  std::string*       miss_reason)
{
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) { if (miss_reason) *miss_reason = "not_found"; return false; }

    std::string kind; int sv = 0; std::string khash, kjson;
    if (!read_str_attr(file, "tictac_cache_kind", kind) || kind != "p123") {
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

    uint64_t nnz_v = 0;
    if (H5LTread_dataset(file, "/nnz", H5T_NATIVE_UINT64, &nnz_v) < 0) {
        if (miss_reason) *miss_reason = "corrupt";
        H5Fclose(file); return false;
    }
    out->nnz = (size_t)nnz_v;
    out->val_array = new double[out->nnz];
    out->row_array = new int[out->nnz];
    out->col_array = new int[out->nnz];
    if (H5LTread_dataset(file, "/val", H5T_NATIVE_DOUBLE, out->val_array) < 0
     || H5LTread_dataset(file, "/row", H5T_NATIVE_INT, out->row_array) < 0
     || H5LTread_dataset(file, "/col", H5T_NATIVE_INT, out->col_array) < 0) {
        delete[] out->val_array; delete[] out->row_array; delete[] out->col_array;
        out->val_array = nullptr; out->row_array = nullptr; out->col_array = nullptr;
        if (miss_reason) *miss_reason = "corrupt";
        H5Fclose(file); return false;
    }

    H5Fclose(file);
    return true;
}

}  // namespace tictac::cache
