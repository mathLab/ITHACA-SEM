# About

This repository contains the `NekPy` Python wrappers for the Nektar++
spectral/_hp_ element framework. **As a disclaimer, thhese wrappings are
_experimental_ and _incomplete_.** You should not rely on their current
structure and API remaining unchanged.

Currently, representative classes from the `LibUtilities`, `StdRegions`,
`SpatialDomains` and `LocalRegions` libraries have been wrapped in order to show
the proof-of-concept.

# Features and functionality

`NekPy` uses the `Boost.Python` library to provide a set of high-quality,
hand-written Python bindings for selected functions and classes in Nektar++. A
typical snippet could look something like:

```python
from NekPy.LibUtilities import PointsKey, PointsType, BasisKey, BasisType
from NekPy.StdRegions import StdQuadExp
import numpy as np
numModes = 8
numPts   = 9
ptsKey   = PointsKey(numPts, PointsType.GaussLobattoLegendre)
basisKey = BasisKey(BasisType.Modified_A, numModes, ptsKey)
quadExp  = StdQuadExp(basisKey, basisKey)
x, y     = quadExp.GetCoords()
fx       = np.sin(x) * np.cos(y)
proj     = quadExp.FwdTrans(fx)
```

`NekPy` uses the `Boost.NumPy` library, contained in Boost 1.63+, to
automatically convert C++ `Array<OneD, >` objects to and from the commonly-used
`numpy.ndarray` object, which makes the integration more seamless between Python
and C++.

# How do I wrap things?

`Boost.Python` is a pretty comprehensive package and an extended discussion is
really beyond the scope of this project. See `doc/wrapping-guide.md` for some
basic concepts and frequently-encountered issues.

# Compiling

NekPy has the following list of requirements:

- Boost with Python support
- Nektar++ `master` branch compiled from source (i.e. not from packages)
- Python 2.7+ (note that examples rely on Python 2.7)
- NumPy

Most of these can be installed using package managers on various operating
systems, as we describe below. We also have a requirement on the `Boost.NumPy`
package, which is available in boost 1.63 or later. If this isn't found on your
system, it will be automatically downloaded and compiled.

## Compiling and installing Nektar++

Nektar++ should be compiled as per the user guide instructions and installed
into a directory which we will refer to as `$NEKDIR`. By default this is the
`dist` directory inside the Nektar++ build directory.

Note that Nektar++ must, at a minimum, be compiled with `NEKTAR_BUILD_LIBRARY`,
`NEKTAR_BUILD_UTILITIES` and `NEKTAR_BUILD_SOLVERS`. Note that all solvers may
be disabled as long as the `NEKTAR_BUILD_SOLVERS` option is set.

## macOS

### Homebrew 
Users of Homebrew should make sure their installation is up-to-date with `brew
upgrade`. Then run

```
brew install python boost-python
```

To install the NumPy package, use the `pip` package manager:

```
pip install numpy
```

### MacPorts

Users of MacPorts should sure their installation is up-to-date with `sudo port
selfupdate && sudo port upgrade outdated`. Then run

```
sudo port install python27 py27-numpy
sudo port select --set python python27
```

### Compiling the wrappers

Clone and compile the wrappers as follows:

```
git clone https://gitlab.nektar.info/nektar/nektar-python.git
cd nektar-python && mkdir builds && cd builds
cmake -DNektar++_DIR=$NEKDIR/lib/nektar++-4.5.0/cmake ..
make install
```

This will automatically download and install `Boost.NumPy` if required.

### Linux: Ubuntu/Debian

Users of Debian and Ubuntu Linux systems should sure their installation is
up-to-date with `sudo apt-get update && sudo apt-get upgrade`

```
sudo apt-get install libboost-python-dev python-numpy
```

# Using the bindings

By default, the bindings will install into the `dist` directory, along with a
number of examples that are stored in the `example` directory. To test your
installation, you can for example run one of these (e.g. `python Basis.py`) or
launch an interactive session:

```
$ cd builds
$ python
Python 2.7.13 (default, Apr  4 2017, 08:47:57) 
[GCC 4.2.1 Compatible Apple LLVM 8.1.0 (clang-802.0.38)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from NekPy.LibUtilities import PointsKey, PointsType
>>> PointsKey(10, PointsType.GaussLobattoLegendre)
<NekPy.LibUtilities._LibUtilities.PointsKey object at 0x11005c310>
```

## Examples

A number of examples of the wrappers can be found in the `dist` directory, along
with a sample mesh `newsquare_2x2.xml`:

- `SessionReader.py` is the simplest example and shows how to construct a
  session reader object. Run it as `python SessionReader.py mesh.xml`.
- `Basis.py` shows functionality of basic `LibUtilities` points and basis
  classes. Run this as `python Basis.py`.
- `StdProject.py` shows how to use some of the `StdRegions` wrappers and
  duplicates the functionality of `Basis.py` using the `StdExpansion` class. Run
  this as `python StdProject.py`.
- `MeshGraph.py` loads a mesh and prints out some basic properties of its
  quadrilateral elements. Run it as `python MeshGraph.py newsquare_2x2.xml`.

If you want to modify the source files, it's advisable to edit them in the
`examples` directory and re-run `make install`, otherwise local changes will be
overwritten by the next `make install`.
