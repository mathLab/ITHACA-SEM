# `nektar` image

This image provides a full installation of Nektar++ with the following options
enabled:

- `NEKTAR_USE_MPI`
- `NEKTAR_USE_FFT`
- `NEKTAR_BUILD_PYTHON`
- `NEKTAR_USE_HDF5`
- `NEKTAR_USE_MESHGEN`
- `NEKTAR_USE_PETSC`
- `NEKTAR_USE_CCM`

## Building

The image is built using `nektarpp/nektar-env` and requires as build context the
path to the Nektar++ source tree. It supports three additional build arguments:

- `BUILD_DEMOS` can be set to `ON` to build demos, which are disabled by
  default;
- `BUILD_SOLVERS` can be set to `OFF` to disable build of solvers, which are
  enabled by default;
- `BUILD_ARGS` can be set to `ON` to build doxygen documentation, which are
  enabled by default;
- `INSTALL_PREFIX` can be set to adjust the install prefix, which is
  `/usr/local` by default.
  
Then build the image using a command similar to:

```sh
docker build -t nektarpp/nektar -f Dockerfile ~/nektar++
```

Note that this is a multi-stage build: after the initial build phase is
completed, the build tree is erased along with development headers, and the
libraries/executables installed into a fresh Debian 10 container (with
appropriate packages for runtime libraries). If you want to keep the build
files, you should instead only build up to the end of the build stage using
`--target build`:

```sh
docker build --target build -t nektar-build -f Dockerfile ~/nektar++
```
