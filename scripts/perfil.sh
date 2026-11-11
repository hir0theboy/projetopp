#!/usr/bin/env bash
# Perfilamento por trecho (init / laco do stencil / soma) do baseline instrumentado.
# Saida: dados/perfil_<MAQUINA>.csv   Uso: MAQUINA=notebook REPS=5 ./scripts/perfil.sh
cd "$(dirname "$0")/.." || exit 1
source scripts/comum.sh
make -s seq || exit 1
MAQ=${MAQUINA:-$(hostname)}; REPS=${REPS:-5}; T=${T:-500}; SIZES=${SIZES:-"500 1000 2000"}
OUT=dados/perfil_${MAQ}.csv
echo "maquina,N,T,rep,t_init,t_laco,t_soma" > "$OUT"
for N in $SIZES; do
  bin/calor2d_ref $N $T > /dev/null            # aquecimento (descartado)
  for rep in $(seq 1 $REPS); do
    s=$(bin/calor2d_ref $N $T)
    echo "$MAQ,$N,$T,$rep,$(get "$s" t_init),$(get "$s" tempo),$(get "$s" t_soma)" >> "$OUT"
  done
done
echo "gravado: $OUT"
