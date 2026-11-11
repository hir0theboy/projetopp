#!/usr/bin/env bash
# Valida que OpenMP e MPI produzem EXATAMENTE a mesma grade final do baseline
# (hash bit a bit) e o mesmo checksum (tolerancia 1e-9 relativa, por causa da ordem da soma).
# Cobre: varios tamanhos, N impar, N nao divisivel por P, caso em que o calor chega na borda,
# varias threads/processos e as tres politicas de escalonamento.
cd "$(dirname "$0")/.." || exit 1
source scripts/comum.sh
make -s all || { echo "falha na compilacao"; exit 1; }

CASOS=("16 10" "50 100" "64 3000" "101 500" "257 300" "500 500")
THREADS=(1 2 3 4 8)
PROCS=(1 2 3 4 5 7)
falhas=0; total=0

confere() {  # rotulo saida_ref saida_teste
  total=$((total+1))
  local hr=$(get "$2" hash) ht=$(get "$3" hash) cr=$(get "$2" checksum) ct=$(get "$3" checksum)
  local ok_hash=0 ok_chk
  [ -n "$ht" ] && [ "$hr" = "$ht" ] && ok_hash=1
  ok_chk=$(awk -v a="$cr" -v b="$ct" 'BEGIN{d=a-b; if(d<0)d=-d; m=(a<0?-a:a); if(m<1)m=1; print (d/m<1e-9)?1:0}')
  if [ $ok_hash = 1 ] && [ "$ok_chk" = 1 ]; then printf "  OK    %s\n" "$1"
  else printf "  FALHA %s  (ref hash=%s chk=%s | teste hash=%s chk=%s)\n" "$1" "$hr" "$cr" "$ht" "$ct"; falhas=$((falhas+1)); fi
}

for caso in "${CASOS[@]}"; do
  set -- $caso; N=$1; T=$2
  echo "== N=$N passos=$T"
  ref=$(bin/calor2d_ref $N $T)
  orig=$(bin/calor2d_seq $N $T)
  total=$((total+1))
  if [ "$(get "$orig" checksum)" = "$(get "$ref" checksum)" ]; then echo "  OK    instrumentada == original do professor (checksum)"
  else echo "  FALHA instrumentada difere do original"; falhas=$((falhas+1)); fi
  for s in static dynamic guided; do
    for t in "${THREADS[@]}"; do
      out=$(OMP_NUM_THREADS=$t OMP_SCHEDULE=$s bin/calor2d_omp $N $T)
      confere "OpenMP  schedule=$s threads=$t" "$ref" "$out"
    done
  done
  for v in omp_atomic omp_critical; do
    out=$(OMP_NUM_THREADS=4 bin/calor2d_$v $N $T); confere "OpenMP  $v threads=4" "$ref" "$out"
  done
  for p in "${PROCS[@]}"; do
    [ $p -gt $N ] && continue
    out=$(mpi_run $p bin/calor2d_mpi $N $T 2>&1)
    confere "MPI     np=$p" "$ref" "$out"
  done
done
echo
if [ $falhas -eq 0 ]; then echo "TODOS OS $total TESTES PASSARAM"; exit 0
else echo "$falhas FALHA(S) em $total testes"; exit 1; fi
