FROM jupyter/scipy-notebook:ae5f7e104dd5 AS build

USER root
RUN apt-get update && \
    apt-get install -y libtinyxml-dev && \
    rm -rf /var/lib/apt/lists/*

USER $NB_UID
RUN conda install --quiet --yes cmake==3.14.3 boost==1.69 jupyter_contrib_nbextensions && \
    jupyter contrib nbextension install --sys-prefix && \
    jupyter nbextension enable --sys-prefix hide_input/main && \
    jupyter nbextension enable --sys-prefix hide_input_all/main && \
    jupyter nbextension enable --sys-prefix init_cell/main && \
    jupyter nbextension enable --sys-prefix equation-numbering/main

USER root
COPY . /tmp/nektar
RUN cd /tmp/nektar && mkdir build && cd build && \
    cmake -DNEKTAR_BUILD_SOLVERS=OFF -DNEKTAR_BUILD_UTILITIES=ON -DNEKTAR_BUILD_PYTHON=ON \
          -DNEKTAR_USE_PYTHON3=ON -DNEKTAR_BUILD_TESTS=OFF -DNEKTAR_BUILD_UNIT_TESTS=OFF \
          -DNEKTAR_BUILD_DEMOS=OFF -DNEKTAR_USE_SCOTCH=OFF -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) install && \
    chown -R $NB_UID:$NB_UID /tmp/nektar && \
    cd /tmp/nektar/build && sudo -u $NB_USER make nekpy-install-system && cd && rm -R /tmp/nektar

USER $NB_UID