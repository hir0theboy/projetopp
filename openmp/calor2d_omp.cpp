/* Versao OpenMP (memoria compartilhada) — equacao do calor 2D.
 *
 * Padrao: laco paralelo (omp for) sobre as linhas da grade, com a regiao
 * paralela ABRANGENDO o laco de tempo (um unico fork/join, em vez de um por passo).
 *
 * Condicoes de corrida:
 *  - Stencil: NAO ha corrida. Cada passo le so `cur` e escreve so `nxt`
 *    (dois buffers), e cada thread escreve linhas distintas. A barreira implicita
 *    ao fim do `omp for` garante que o passo t termina antes do t+1 comecar.
 *  - Checksum: a soma de todas as celulas EH uma corrida em potencial (todos escrevem
 *    em `soma`). Tres solucoes selecionaveis na compilacao:
 *        (padrao)         reduction(+:soma)   -> soma parcial privada por thread
 *        -DSOMA_ATOMIC    #pragma omp atomic
 *        -DSOMA_CRITICAL  #pragma omp critical
 *        -DSOMA_RACE      SEM protecao (so demonstracao: resultado errado/variavel)
 *
 * Escalonamento: schedule(runtime) -> escolha por OMP_SCHEDULE=static|dynamic|guided[,chunk].
 * Primeiro acesso (first-touch): a grade e inicializada em paralelo com schedule(static),
 * para que as paginas de memoria fiquem proximas da thread que as usara (NUMA).
 */
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <omp.h>
#include "../comum/util.h"
using namespace std;

int main(int argc, char** argv) {
    int N = (argc > 1) ? atoi(argv[1]) : 500;
    int T = (argc > 2) ? atoi(argv[2]) : 500;
    size_t n = (size_t)N * N;

    double* A = new double[n];
    double* B = new double[n];

    double ta = agora();
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) { A[(size_t)y*N + x] = 0.0; B[(size_t)y*N + x] = 0.0; }
    A[(size_t)(N/2)*N + N/2] = 100.0;          // fonte quente no centro
    double tb = agora();

    #pragma omp parallel
    {
        double* cur = A;                        // ponteiros PRIVADOS de cada thread
        double* nxt = B;
        for (int t = 0; t < T; t++) {
            #pragma omp for schedule(runtime)   // barreira implicita no fim
            for (int y = 1; y < N-1; y++)
                for (int x = 1; x < N-1; x++)
                    nxt[(size_t)y*N + x] = 0.25 * (cur[(size_t)(y-1)*N+x] + cur[(size_t)(y+1)*N+x] +
                                                   cur[(size_t)y*N+(x-1)] + cur[(size_t)y*N+(x+1)]);
            swap(cur, nxt);
        }
    }
    double tc = agora();

    const double* res = (T % 2 == 0) ? A : B;   // apos T trocas, o resultado esta aqui

    double soma = 0.0;
#if defined(SOMA_ATOMIC)
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)n; i++) {
        #pragma omp atomic
        soma += res[i];
    }
#elif defined(SOMA_CRITICAL)
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)n; i++) {
        #pragma omp critical
        soma += res[i];
    }
#elif defined(SOMA_RACE)
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)n; i++) soma += res[i];   // CORRIDA DE DADOS (demonstracao)
#else
    #pragma omp parallel for schedule(static) reduction(+:soma)
    for (long i = 0; i < (long)n; i++) soma += res[i];
#endif
    double td = agora();

    uint64_t h = hash_grade(res, n);
    printf("N=%d passos=%d  checksum=%.6f  tempo=%.6f s  t_init=%.6f t_soma=%.6f threads=%d hash=%016llx\n",
           N, T, soma, tc - tb, tb - ta, td - tc, omp_get_max_threads(), (unsigned long long)h);
    delete[] A; delete[] B;
    return 0;
}
