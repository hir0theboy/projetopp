/* Baseline 05 — Equacao do calor 2D por diferencas finitas (stencil iterativo)
 * IESB 2026/2 — CCO085 — Prof. Rodrigo Goncalves Pinto
 * Aspecto interessante: stencil iterativo com dependencia entre passos;
 * caso classico de MPI com troca de bordas (halo exchange) a cada iteracao.
 *
 * Compilar: g++ -O2 -o 05_calor2d 05_calor2d.cpp
 * Executar: ./05_calor2d 500 500   (lado da grade; num. de passos de tempo)
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
using namespace std;

static double agora() {
    timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(int argc, char** argv) {
    int N = (argc > 1) ? atoi(argv[1]) : 500;
    int T = (argc > 2) ? atoi(argv[2]) : 500;
    vector<double> u(N*N, 0.0), un(N*N, 0.0);

    // condicao inicial: fonte quente no centro
    u[(N/2)*N + N/2] = 100.0;

    double t0 = agora();
    for (int t = 0; t < T; t++) {
        for (int y = 1; y < N-1; y++)
            for (int x = 1; x < N-1; x++)
                un[y*N + x] = 0.25 * (u[(y-1)*N+x] + u[(y+1)*N+x] +
                                      u[y*N+(x-1)] + u[y*N+(x+1)]);
        u.swap(un);
    }
    double t1 = agora();

    double soma = 0;
    for (int i = 0; i < N*N; i++) soma += u[i];
    printf("N=%d passos=%d  checksum=%.6f  tempo=%.6f s\n", N, T, soma, t1 - t0);
    return 0;
}
