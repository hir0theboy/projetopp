/* Versao MPI (memoria distribuida) — equacao do calor 2D.
 *
 * DECOMPOSICAO: 1D por faixas de linhas (row-block). As N linhas sao divididas em P
 * blocos contiguos (os primeiros N%P ranks recebem uma linha a mais). Cada rank guarda
 * suas linhas + 2 linhas FANTASMA (halo): a de cima e a de baixo, copias das linhas de
 * borda dos vizinhos. Linhas sao contiguas na memoria (row-major), entao o halo e um
 * unico bloco de N doubles, sem tipo derivado.
 *
 * TROCA DE BORDAS (halo exchange), a cada passo, com MPI_Sendrecv (sem risco de deadlock):
 *    1) envia minha 1a linha para o vizinho de cima  / recebe a fantasma de baixo
 *    2) envia minha ultima linha ao vizinho de baixo / recebe a fantasma de cima
 *  Nas pontas usa-se MPI_PROC_NULL (a operacao vira no-op).
 *
 * COLETIVAS: MPI_Bcast (parametros), MPI_Reduce (checksum e tempos), MPI_Gatherv (grade
 *  final no rank 0, so para validar por hash — fora da regiao cronometrada).
 *
 * Custo por passo e por rank: computacao ~ (N/P)*N ; comunicacao ~ 2N doubles (constante
 * em P). Por isso a eficiencia cai quando P cresce com N fixo.
 *
 * A logica numerica (stencil, bordas fixas em 0, fonte no centro) e a do baseline.
 */
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "../comum/util.h"
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    // --- parametros: lidos no rank 0 e distribuidos por broadcast
    int par[2] = {0, 0};
    if (rank == 0) { par[0] = (argc > 1) ? atoi(argv[1]) : 500; par[1] = (argc > 2) ? atoi(argv[2]) : 500; }
    MPI_Bcast(par, 2, MPI_INT, 0, MPI_COMM_WORLD);
    const int N = par[0], T = par[1];

    if (P > N) {
        if (rank == 0) fprintf(stderr, "erro: processos (%d) > linhas (%d)\n", P, N);
        MPI_Finalize();
        return 1;
    }

    // --- decomposicao por linhas
    const int base = N / P, rem = N % P;
    const int rows  = base + (rank < rem ? 1 : 0);        // linhas deste rank
    const int start = rank * base + min(rank, rem);       // 1a linha global deste rank
    vector<int> counts(P), displs(P);
    for (int r = 0; r < P; r++) {
        counts[r] = (base + (r < rem ? 1 : 0)) * N;
        displs[r] = (r * base + min(r, rem)) * N;
    }

    // local: (rows+2) linhas; linha 0 e rows+1 sao fantasmas; linha i <-> global start+i-1
    vector<double> a((size_t)(rows + 2) * N, 0.0), b((size_t)(rows + 2) * N, 0.0);
    double* cur = a.data();
    double* nxt = b.data();

    const int gy = N / 2;                                 // fonte quente no centro
    if (gy >= start && gy < start + rows) cur[(size_t)(gy - start + 1) * N + N / 2] = 100.0;

    const int up   = (rank > 0)     ? rank - 1 : MPI_PROC_NULL;
    const int down = (rank < P - 1) ? rank + 1 : MPI_PROC_NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    double tcomm = 0.0;

    for (int t = 0; t < T; t++) {
        double c0 = MPI_Wtime();
        MPI_Sendrecv(cur + (size_t)1 * N,        N, MPI_DOUBLE, up,   0,
                     cur + (size_t)(rows + 1) * N, N, MPI_DOUBLE, down, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(cur + (size_t)rows * N,     N, MPI_DOUBLE, down, 1,
                     cur,                        N, MPI_DOUBLE, up,   1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        tcomm += MPI_Wtime() - c0;

        for (int i = 1; i <= rows; i++) {
            int g = start + i - 1;                         // linha global
            if (g == 0 || g == N - 1) continue;            // bordas fixas em 0
            for (int x = 1; x < N - 1; x++)
                nxt[(size_t)i*N + x] = 0.25 * (cur[(size_t)(i-1)*N+x] + cur[(size_t)(i+1)*N+x] +
                                               cur[(size_t)i*N+(x-1)] + cur[(size_t)i*N+(x+1)]);
        }
        swap(cur, nxt);
    }
    double tloc = MPI_Wtime() - t0;

    // --- checksum: soma local das linhas proprias + MPI_Reduce
    double sloc = 0.0, soma = 0.0;
    for (int i = 1; i <= rows; i++)
        for (int x = 0; x < N; x++) sloc += cur[(size_t)i*N + x];
    MPI_Reduce(&sloc, &soma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // --- tempos: o do programa e o do rank mais lento; comunicacao: max e media
    double tmax = 0, cmax = 0, csum = 0;
    MPI_Reduce(&tloc,  &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tcomm, &cmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tcomm, &csum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // --- validacao: junta a grade no rank 0 e calcula o mesmo hash do baseline
    vector<double> full;
    if (rank == 0) full.resize((size_t)N * N);
    MPI_Gatherv(cur + N, rows * N, MPI_DOUBLE,
                rank == 0 ? full.data() : nullptr, counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        uint64_t h = hash_grade(full.data(), (size_t)N * N);
        printf("N=%d passos=%d  checksum=%.6f  tempo=%.6f s  comm_max=%.6f comm_media=%.6f procs=%d hash=%016llx\n",
               N, T, soma, tmax, cmax, csum / P, P, (unsigned long long)h);
    }
    MPI_Finalize();
    return 0;
}
