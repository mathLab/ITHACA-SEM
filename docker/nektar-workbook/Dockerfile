FROM jupyter/scipy-notebook:42f4c82a07ff AS build

LABEL maintainer="Nektar++ Development Team <nektar-users@imperial.ac.uk>"

USER root
RUN apt-get update && \
    apt-get install -y libgl1-mesa-dev xvfb && \
    rm -rf /var/lib/apt/lists/*

USER $NB_UID
RUN conda install --quiet --yes cmake=3.19.1 boost=1.74.0 occt=7.4.0 gmsh=4.6.0 pyvista=0.27.2 itkwidgets=0.32.0 arpack=3.7.0 fftw=3.3.8 tinyxml jupyter_contrib_nbextensions && \
    jupyter contrib nbextension install --sys-prefix && \
    jupyter nbextension enable --sys-prefix hide_input/main && \
    jupyter nbextension enable --sys-prefix hide_input_all/main && \
    jupyter nbextension enable --sys-prefix init_cell/main && \
    jupyter nbextension enable --sys-prefix equation-numbering/main

USER root
COPY . /tmp/nektar
RUN cd /tmp/nektar && mkdir build && cd build && \
    cmake -DNEKTAR_BUILD_SOLVERS=ON -DNEKTAR_SOLVER_DIFFUSION=OFF \
          -DNEKTAR_SOLVER_DUMMY=OFF -DNEKTAR_SOLVER_ELASTICITY=OFF \
          -DNEKTAR_SOLVER_MMF=OFF -DNEKTAR_SOLVER_CARDIAC_EP=OFF \
          -DNEKTAR_BUILD_UTILITIES=ON -DNEKTAR_BUILD_PYTHON=ON \
          -DNEKTAR_SOLVER_VORTEXWAVE=OFF -DNEKTAR_LIB_DIR=lib/x86_64-linux-gnu \
          -DNEKTAR_BUILD_TESTS=OFF -DNEKTAR_BUILD_UNIT_TESTS=OFF \
          -DNEKTAR_BUILD_DEMOS=OFF -DNEKTAR_USE_SCOTCH=OFF -DNEKTAR_USE_MESHGEN=ON \
          -DNEKTAR_USE_ARPACK=ON -DNEKTAR_USE_FFTW=ON \
          -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) install && \
    chown -R $NB_UID:$NB_UID /tmp/nektar && \
    cd /tmp/nektar/build && sudo -u $NB_USER make nekpy-install-system && cd && ldconfig && rm -R /tmp/nektar

USER $NB_UID
