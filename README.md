# Equacao do calor 2D — OpenMP e MPI

Paralelizacao do stencil de Jacobi de 5 pontos (`sequencial/05_calor2d.cpp`, fornecido pelo professor).
Todas as versoes imprimem uma linha `chave=valor` com `checksum`, `tempo` (so o laco de passos) e `hash`
(FNV-1a 64 bits da grade final). **O hash igual ao do baseline prova que a grade e identica bit a bit.**

## Estrutura
```
sequencial/05_calor2d.cpp   baseline do professor (INTOCADO)
sequencial/calor2d_ref.cpp  copia instrumentada (cronometro por trecho + hash); mesma logica
openmp/calor2d_omp.cpp      versao OpenMP (+ variantes de soma: atomic, critical, race)
mpi/calor2d_mpi.cpp         versao MPI (decomposicao 1D por linhas, halo exchange, coletivas)
comum/util.h                cronometro e hash
scripts/                    validar.sh  perfil.sh  medir.sh  plotar.py
dados/  graficos/           CSVs de medicao e figuras
```

## Ambiente
Referencia: Linux + g++ + OpenMPI (laboratorio). Em casa, use o contêiner:

```bash
docker build -t calor2d .
docker run -it --rm -v "$PWD":/workspace calor2d      # (Windows PowerShell: -v ${PWD}:/workspace)
```
Ou abra a pasta no VS Code e use **Dev Containers: Reopen in Container** (usa `.devcontainer/`).
No Docker Desktop (Windows/macOS), aloque todos os CPUs em Settings > Resources.
macOS sem Docker: `brew install gcc open-mpi` e `make CXX=g++-14`.

## Uso (do zero)
```bash
make            # compila tudo em bin/
make test       # valida OpenMP e MPI contra o baseline (hash exato) — deve terminar com "TODOS OS ... PASSARAM"
make perfil     # perfilamento do baseline -> dados/perfil_<maquina>.csv
make medir      # medicoes -> dados/medicoes_<maquina>.csv
make graficos   # figuras em graficos/ e tabelas em dados/
```
Use `MAQUINA=notebook make medir` para rotular a maquina (o CSV tem uma coluna `maquina`).
Variaveis de `medir.sh`: `REPS` (5), `T` (500), `SIZES` ("500 1000 2000"), `NW0` (1000), `NCORES`, `MPIRUN_EXTRA`.
Protocolo: 1 execucao de aquecimento descartada + `REPS` registradas; graficos usam a mediana.

Execucao manual:
```bash
bin/calor2d_seq 1000 500
OMP_NUM_THREADS=4 OMP_SCHEDULE=static bin/calor2d_omp 1000 500
mpirun -np 4 bin/calor2d_mpi 1000 500
OMP_NUM_THREADS=8 bin/calor2d_omp_race 1000 20000   # demonstracao de corrida de dados (checksum varia)
```
