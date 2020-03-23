# `nektar-env` image

This image is designed to provide a build environment for Nektar++ based on the
Debian 10 (buster) build image. It installs all libraries to enable Nektar++ to
be compiled with most/all third-party dependencies turned on. In particular we
install the following development libraries:

- Boost
- TinyXML
- LAPACK/BLAS
- OpenMPI (enables `NEKTAR_USE_MPI`)
- FFTW (enables `NEKTAR_USE_FFT`)
- Python and Boost.Python (enables `NEKTAR_BUILD_PYTHON`)
- HDF5 (enables `NEKTAR_USE_HDF5`)
- OCE, Triangle and TetGen (enables `NEKTAR_USE_MESHGEN`)
- PETSc (enables `NEKTAR_USE_PETSC`)
- ARPACK (enables `NEKTAR_USE_ARPACK`)

## Building

No particular context is required to build this image. Building using the below
or similar.

```sh
docker build -t nektarpp/nektar-env -f Dockerfile .
```

# Other environment images

The other dockerfiles for different operating systems and package lists are used
to provide environment images for the CI system. See .gitlab-ci.yml for details.
