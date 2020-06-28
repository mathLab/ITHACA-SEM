# `nektar-workbook` image

This image is based on the `scipy-notebook` image but augments it with the
`NekPy` Python interface. It also includes the `NekMesh` and `FieldConvert`
utilities for post-processing.

## Building

The image requires as build context the path to the Nektar++ source tree.  Build
the image using a command similar to:

```sh
docker build -t nektarpp/nektar-notebook -f Dockerfile ~/nektar++
```
