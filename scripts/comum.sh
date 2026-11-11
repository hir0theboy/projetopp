# Funcoes compartilhadas pelos scripts (use: source scripts/comum.sh)
# get SAIDA CHAVE -> valor do campo chave=valor da linha de saida dos programas
get() { grep -oE "(^| )$2=[^ ]+" <<<"$1" | head -1 | sed 's/.*=//'; }

# opcoes padrao do mpirun: --oversubscribe permite mais processos que "slots"
# (OpenMPI conta nucleos fisicos); --bind-to none evita erro de binding em contêiner.
MPIRUN_EXTRA=${MPIRUN_EXTRA:-"--bind-to none"}
[ "$(id -u)" = "0" ] && MPIRUN_EXTRA="$MPIRUN_EXTRA --allow-run-as-root"
mpi_run() { local np=$1; shift; mpirun $MPIRUN_EXTRA --oversubscribe -np "$np" "$@"; }

nucleos() { nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1; }
