// Utilidades compartilhadas pelas tres versoes (sequencial instrumentada, OpenMP, MPI).
// Nada aqui altera a logica de calculo: so cronometro e verificacao do resultado.
#pragma once
#include <cstdint>
#include <cstddef>
#include <ctime>

static inline double agora() {
    timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

// Hash FNV-1a (64 bits) sobre os bytes da grade, em ordem linha a linha.
// Duas grades com o mesmo hash sao identicas bit a bit (validacao exata,
// mais forte que o checksum, que conserva massa e quase nao discrimina erros).
static inline uint64_t hash_grade(const double* p, size_t n) {
    const unsigned char* b = reinterpret_cast<const unsigned char*>(p);
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n * sizeof(double); i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}
