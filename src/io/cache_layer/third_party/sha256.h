// Public-domain SHA-256, adapted from B-Con's sha-2 implementation.
#pragma once
#include <cstdint>
#include <cstddef>
#include <string>

namespace tictac::cache::sha256 {

struct Ctx {
    uint8_t  data[64];
    uint32_t datalen;
    uint64_t bitlen;
    uint32_t state[8];
};

void init(Ctx& ctx);
void update(Ctx& ctx, const uint8_t* data, size_t len);
void finalize(Ctx& ctx, uint8_t hash[32]);

// Convenience: returns 64-char lowercase hex.
std::string hex_digest(const std::string& message);

}  // namespace tictac::cache::sha256
