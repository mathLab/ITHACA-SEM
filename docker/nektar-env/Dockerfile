FROM debian:buster-slim

LABEL maintainer="David Moxey <d.moxey@exeter.ac.uk>"

RUN DEBIAN_RELEASE=buster && \
    echo "deb http://deb.debian.org/debian ${DEBIAN_RELEASE} non-free" > \
        /etc/apt/sources.list.d/debian-non-free.list && \
    apt-get update && \
    apt-get install -y build-essential cmake git \
    libboost-thread-dev libboost-iostreams-dev libboost-filesystem-dev libboost-system-dev \
    libboost-program-options-dev libboost-python-dev libboost-numpy-dev libboost-regex-dev \
    libfftw3-dev libhdf5-dev libhdf5-mpi-dev liblapack-dev liboce-foundation-dev \
    liboce-modeling-dev liboce-ocaf-dev libtet1.5-dev libtriangle-dev \
    libtinyxml-dev mpi-default-dev petsc-dev zlib1g-dev python-numpy libptscotch-dev \
    libarpack2-dev \
    texlive-base texlive-latex-extra texlive-science texlive-fonts-recommended texlive-pstricks \
    doxygen graphviz imagemagick

RUN groupadd nektar && useradd -m -g nektar nektar
USER nektar:nektar
WORKDIR /home/nektar
