FROM nektarpp/nektar-env:latest AS build

LABEL maintainer="David Moxey <d.moxey@exeter.ac.uk>"

COPY --chown=nektar:nektar . /home/nektar

ARG BUILD_DOCS=ON
ARG BUILD_DEMOS=OFF
ARG BUILD_SOLVERS=ON
ARG INSTALL_PREFIX=/usr/local

RUN mkdir build && cd build &&                                  \
    cmake                                                       \
    -DNEKTAR_BUILD_DEMOS=${BUILD_DEMOS}                         \
    -DNEKTAR_BUILD_DOC=${BUILD_DOCS}                            \
    -DNEKTAR_BUILD_TESTS=OFF                                    \
    -DNEKTAR_BUILD_UNIT_TESTS=OFF                               \
    -DNEKTAR_BUILD_SOLVERS=${BUILD_SOLVERS}                     \
    -DNEKTAR_BUILD_PYTHON=ON                                    \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}                    \
    -DNEKTAR_USE_MESHGEN=ON                                     \
    -DNEKTAR_USE_MPI=ON                                         \
    -DNEKTAR_USE_ARPACK=ON                                      \
    -DNEKTAR_USE_FFTW=ON                                        \
    -DNEKTAR_USE_HDF5=ON                                        \
    -DNEKTAR_USE_CCM=ON                                         \
    .. &&                                                       \
    make -j$(nproc) &&                                          \
    make -j$(nproc) user-guide-pdf developer-guide-pdf          \
        tutorials-pdf doc

USER root
RUN cd build && make install && make nekpy-install-system && cd ..
USER nektar

FROM debian:buster-slim
RUN DEBIAN_RELEASE=buster && \
    echo "deb http://deb.debian.org/debian ${DEBIAN_RELEASE} non-free" > \
        /etc/apt/sources.list.d/debian-non-free.list && \
    apt-get update && \
    apt-get install -y \
    openssh-client nano vim binutils \
    libboost-thread1.67 libboost-iostreams1.67 libboost-filesystem1.67 libboost-system1.67 \
    libboost-program-options1.67 libboost-python1.67 libboost-regex1.67 libboost-numpy1.67 \
    libfftw3-3 libhdf5-103 libhdf5-openmpi-103 liblapack3 liboce-foundation11 \
    liboce-modeling11 liboce-ocaf11 libtet1.5 libtriangle-1.6 python-numpy \
    openmpi-bin libpetsc3.10 libptscotch-6.0 libarpack2 libtinyxml2.6.2v5 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
ARG INSTALL_PREFIX=/usr/local
COPY --from=build ${INSTALL_PREFIX} ${INSTALL_PREFIX}

# Set up entrypoint for copying test files.
COPY docker/nektar/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && \
    ln -s /usr/local/bin/docker-entrypoint.sh /
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

RUN groupadd nektar && useradd -m -g nektar nektar
USER nektar
WORKDIR /home/nektar
CMD ["/bin/bash"]
