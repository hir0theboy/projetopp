#!/usr/bin/env bash
# Gera os dados de medicao em CSV: dados/medicoes_<MAQUINA>.csv
# Uso:  MAQUINA=notebook REPS=5 ./scripts/medir.sh
# Variaveis: REPS (repeticoes, padrao 5) T (passos, 500) SIZES ("500 1000 2000")
#            NW0 (tamanho base da escala fraca, 1000) NCORES (padrao: nproc)
#            MPIRUN_EXTRA (padrao "--bind-to none"; no lab tente "--bind-to core")
# Protocolo: 1 execucao de aquecimento descartada + REPS execucoes registradas por configuracao.
cd "$(dirname "$0")/.." || exit 1
source scripts/comum.sh
MAQ=${MAQUINA:-$(hostname)}; REPS=${REPS:-5}; T=${T:-500}
SIZES=${SIZES:-"500 1000 2000"}; NW0=${NW0:-1000}; NCORES=${NCORES:-$(nucleos)}
OUT=dados/medicoes_${MAQ}.csv
make -s all || exit 1

WORKERS=(); p=1
while [ $p -le "$NCORES" ]; do WORKERS+=($p); p=$((p*2)); done
[ "${WORKERS[${#WORKERS[@]}-1]}" -ne "$NCORES" ] && WORKERS+=("$NCORES")
echo "maquina=$MAQ nucleos=$NCORES workers=${WORKERS[*]} reps=$REPS T=$T sizes=$SIZES"

echo "maquina,tipo,versao,variante,N,T,workers,rep,tempo,t_soma,comm_max,comm_media,checksum,hash" > "$OUT"

bench() {   # tipo versao variante N T workers -- comando...
  local tipo=$1 versao=$2 var=$3 n=$4 t=$5 w=$6; shift 7
  "$@" > /dev/null 2>&1                                   # aquecimento
  local rep out
  for rep in $(seq 1 "$REPS"); do
    out=$("$@" 2>&1) || { echo "FALHA: $*" >&2; continue; }
    echo "$MAQ,$tipo,$versao,$var,$n,$t,$w,$rep,$(get "$out" tempo),$(get "$out" t_soma),$(get "$out" comm_max),$(get "$out" comm_media),$(get "$out" checksum),$(get "$out" hash)" >> "$OUT"
  done
  echo "  $tipo $versao/$var N=$n w=$w ok"
}

echo ">> escalabilidade FORTE (N fixo, workers variando)"
for N in $SIZES; do
  bench forte seq - $N $T 1 -- bin/calor2d_seq $N $T
  for w in "${WORKERS[@]}"; do
    bench forte omp static $N $T $w -- env OMP_NUM_THREADS=$w OMP_SCHEDULE=static bin/calor2d_omp $N $T
    bench forte mpi -      $N $T $w -- mpirun $MPIRUN_EXTRA --oversubscribe -np $w bin/calor2d_mpi $N $T
  done
done

echo ">> escalabilidade FRACA (celulas por worker constantes: N = NW0*sqrt(w))"
for w in "${WORKERS[@]}"; do
  Nw=$(awk -v n=$NW0 -v w=$w 'BEGIN{printf "%d", n*sqrt(w)+0.5}')
  bench fraca omp static $Nw $T $w -- env OMP_NUM_THREADS=$w OMP_SCHEDULE=static bin/calor2d_omp $Nw $T
  bench fraca mpi -      $Nw $T $w -- mpirun $MPIRUN_EXTRA --oversubscribe -np $w bin/calor2d_mpi $Nw $T
done

echo ">> politicas de escalonamento (OpenMP, maior N, todos os nucleos)"
NMAX=$(echo $SIZES | tr ' ' '\n' | sort -n | tail -1)
rotulos=(static dynamic dynamic16 guided); valores=("static" "dynamic" "dynamic,16" "guided")
for i in 0 1 2 3; do
  bench sched omp ${rotulos[$i]} $NMAX $T $NCORES -- env OMP_NUM_THREADS=$NCORES OMP_SCHEDULE="${valores[$i]}" bin/calor2d_omp $NMAX $T
done

echo ">> condicao de corrida no checksum: reduction vs atomic vs critical (t_soma)"
NR=$(echo $SIZES | tr ' ' '\n' | sort -n | head -2 | tail -1)
for v in omp omp_atomic omp_critical; do
  bench reducao $v - $NR $T $NCORES -- env OMP_NUM_THREADS=$NCORES bin/calor2d_$v $NR $T
done
echo "gravado: $OUT"
