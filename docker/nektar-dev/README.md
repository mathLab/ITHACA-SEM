# `nektar-dev` image

This image provides a full installation of Nektar++, based on the `nektar`
image, but includes the full set of development headers from the `nektar-env`
image for code development purposes.

## Building

The image is built using `nektarpp/nektar-env` or similar and requires as build
context the path to the Nektar++ source tree. It supports the build environment
variables:

- `ENV_IMAGE` and `NEKTAR_IMAGE`, which are used to select the environment image
  and Nektar image to build against. This is used by the CI to e.g. consistently
  build against the correct commits. By default these are set to
  `nektarpp/nektar-env:latest` and `nektarpp/nektar:latest`.
- `INSTALL_PREFIX` can be set to adjust the install prefix, which is
  `/usr/local` by default.
  
Then build the image using a command similar to:

```sh
docker build \
    -t nektarpp/nektar-dev \
    -f ~/nektar++/docker/nektar-dev/Dockerfile \
    --build-arg ENV_IMAGE=nektarpp/nektar-env:latest \
    --build-arg NEKTAR_IMAGE=nektarpp/nektar:latest \
    ~/nektar++
```
