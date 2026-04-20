#include <cstdio>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <psi_3H_file>\n", argv[0]);
        return 1;
    }
    std::printf("check_3nf_normalization: input=%s\n", argv[1]);
    return 0;
}
