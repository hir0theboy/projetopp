# Ambiente de referencia da disciplina: Linux + g++ + OpenMPI (+ Python p/ graficos)
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential g++ make gdb git ca-certificates \
        openmpi-bin libopenmpi-dev \
        python3 python3-pandas python3-matplotlib python3-numpy \
        time less nano zip unzip \
    && rm -rf /var/lib/apt/lists/*

# OpenMPI recusa rodar como root; cria usuario comum (UID 1000, sem conflito com o 'ubuntu' da imagem)
RUN userdel -r ubuntu 2>/dev/null || true \
    && useradd -m -s /bin/bash -u 1000 dev

# Em contêiner: permite mais processos que "slots" e evita aviso do mecanismo de memoria compartilhada
ENV OMPI_MCA_rmaps_base_oversubscribe=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

USER dev
WORKDIR /workspace
CMD ["bash"]
