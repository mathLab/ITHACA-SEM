# `nektar-dev` image

This image provides a full installation of Nektar++, based on the `nektar`
image, but includes the full set of development headers from the `nektar-env`
image for code development purposes.

## Building

The image is built using `nektarpp/nektar-env` and requires as build context the
path to the Nektar++ source tree. It supports the build environment variable:

- `INSTALL_PREFIX` can be set to adjust the install prefix, which is
  `/usr/local` by default.
  
Then build the image using a command similar to:

```sh
docker build -t nektarpp/nektar-dev -f ~/nektar++/docker/nektar-dev/Dockerfile ~/nektar++
```
