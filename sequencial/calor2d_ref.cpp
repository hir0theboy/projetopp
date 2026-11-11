/* Versao sequencial INSTRUMENTADA (copia de 05_calor2d.cpp).
 * A logica de calculo e a entrada sao identicas as do codigo do professor
 * (sequencial/05_calor2d.cpp). Acrescenta apenas:
 *   - cronometro por trecho (init / laco / soma)  -> perfilamento
 *   - hash da grade final                          -> validacao exata das versoes paralelas
 * O campo `tempo` mede so o laco de passos, como no original.
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../comum/util.h"
using namespace std;

int main(int argc, char** argv) {
    int N = (argc > 1) ? atoi(argv[1]) : 500;
    int T = (argc > 2) ? atoi(argv[2]) : 500;

    double ta = agora();
    vector<double> u(N*N, 0.0), un(N*N, 0.0);
    u[(N/2)*N + N/2] = 100.0;
    double tb = agora();

    for (int t = 0; t < T; t++) {
        for (int y = 1; y < N-1; y++)
            for (int x = 1; x < N-1; x++)
                un[y*N + x] = 0.25 * (u[(y-1)*N+x] + u[(y+1)*N+x] +
                                      u[y*N+(x-1)] + u[y*N+(x+1)]);
        u.swap(un);
    }
    double tc = agora();

    double soma = 0;
    for (int i = 0; i < N*N; i++) soma += u[i];
    double td = agora();

    uint64_t h = hash_grade(u.data(), (size_t)N*N);
    printf("N=%d passos=%d  checksum=%.6f  tempo=%.6f s  t_init=%.6f t_soma=%.6f hash=%016llx\n",
           N, T, soma, tc - tb, tb - ta, td - tc, (unsigned long long)h);
    return 0;
}
