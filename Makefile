# Equacao do calor 2D — baseline, OpenMP e MPI
# Ambiente de referencia: Linux + g++ + OpenMPI (laboratorio / contêiner Docker).
# Em outros ambientes: make CXX=g++-14   (macOS com Homebrew)
CXX      ?= g++
MPICXX   ?= mpicxx
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
BIN      := bin

.PHONY: all seq omp mpi test perfil medir graficos clean
all: seq omp mpi

seq: $(BIN)/calor2d_seq $(BIN)/calor2d_ref
omp: $(BIN)/calor2d_omp $(BIN)/calor2d_omp_atomic $(BIN)/calor2d_omp_critical $(BIN)/calor2d_omp_race
mpi: $(BIN)/calor2d_mpi

$(BIN):
	mkdir -p $(BIN)

# baseline do professor, intocado (mesma linha de compilacao do cabecalho do arquivo)
$(BIN)/calor2d_seq: sequencial/05_calor2d.cpp | $(BIN)
	$(CXX) -O2 -o $@ $<

$(BIN)/calor2d_ref: sequencial/calor2d_ref.cpp comum/util.h | $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ $<

OMPSRC = openmp/calor2d_omp.cpp comum/util.h
$(BIN)/calor2d_omp: $(OMPSRC) | $(BIN)
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $<
$(BIN)/calor2d_omp_atomic: $(OMPSRC) | $(BIN)
	$(CXX) $(CXXFLAGS) -fopenmp -DSOMA_ATOMIC -o $@ $<
$(BIN)/calor2d_omp_critical: $(OMPSRC) | $(BIN)
	$(CXX) $(CXXFLAGS) -fopenmp -DSOMA_CRITICAL -o $@ $<
$(BIN)/calor2d_omp_race: $(OMPSRC) | $(BIN)
	$(CXX) $(CXXFLAGS) -fopenmp -DSOMA_RACE -o $@ $<

$(BIN)/calor2d_mpi: mpi/calor2d_mpi.cpp comum/util.h | $(BIN)
	$(MPICXX) $(CXXFLAGS) -o $@ $<

test: all
	./scripts/validar.sh
perfil: seq
	./scripts/perfil.sh
medir: all
	./scripts/medir.sh
graficos:
	python3 scripts/plotar.py
clean:
	rm -rf $(BIN)
