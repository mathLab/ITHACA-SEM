# Wrapping guide

This document attempts to outline some of the basic principles of the NekPy
wrapper, which relies entirely on the excellent `Boost.Python` library. An
extensive documentation is therefore beyond the scope of this document, but we
highlight aspects that are important for the NekPy wrappers.

In general, note that when things go wrong with `Boost.Python`, it'll be
indicated either by an extensive compiler error, or a runtime error in the
Python interpreter when you try to use your wrapper. Judicious use of Google is
therefore recommended to track down these issues!

You may also find the following resources useful:

- [The `Boost.Python` tutorial](http://www.boost.org/doc/libs/1_61_0/libs/python/doc/html/tutorial/index.html)
- The `Boost.Python`
  [entry on the Python wiki](https://wiki.python.org/moin/boost.python/HowTo)
- Examples [on GitHub](https://github.com/TNG/boost-python-examples)
- The
  [OpenStreetGraph cookbook](https://github.com/Skylark13/osgboostpython/blob/wiki/WrappingCookbook.md)
  and
  [rationale for using manual-wrapping](https://github.com/Skylark13/osgboostpython/blob/wiki/ManualWrappingRationale.md)
  which served as a starting point for this project.

# Basic structure

The NekPy wrapper is designed to mimic the library structure of Nektar++, with
directories for the `LibUtilities`, `SpatialDomains` and `StdRegions`
libraries. This is a deliberate design decision, so that classes and definitions
in Nektar++ can be easily located inside NekPy.

There are also some other directories and files:

- `NekPyConfig.hpp` is a convenience header that all `.cpp` files should
  import. It sets appropriate namespaces for `boost::python` and
  `boost::python::numpy`, depending on whether the `Boost.NumPy` library was
  compiled or is included in Boost.
- `lib` is a skeleton Python directory, which will be installed by CMake into
  the `dist` directory and automatically import the compiled binary libraries.
- `example` contains some basic Python examples.
- `cmake` has some CMake configuration files.

To demonstrate how to wrap classes, we'll refer to a number of existing parts of
the code below.

# Defining a library

First consider `LibUtilities`. An abbreviated version of the base file,
`LibUtilities.cpp` has the following structure:

```c++
#include <LibUtilities/Python/NekPyConfig.hpp>

void export_Basis();
void export_SessionReader();

BOOST_PYTHON_MODULE(_LibUtilities)
{
    // Initialise Boost.NumPy.
    np::initialize();

    // Export classes.
    export_Basis();
    export_SessionReader();
}
```

The `BOOST_PYTHON_MODULE(name)` macro allows us to define a Python module inside
C++. Note that in this case, the leading underscore in the name
(i.e. `_LibUtilities`) is deliberate. To define the contents of the module, we
call a number of functions that are prefixed by `export_`, which will define one
or more Python classes that live in this module. These Python classes correspond
with our Nektar++ classes. We adopt this approach simply so that we can split up
the different classes into different files, because it is not possible to call
`BOOST_PYTHON_MODULE` more than once. These functions are defined in
appropriately named files, for example:

- `export_Basis()` lives in the file
  `library/LibUtilities/Python/Foundations/Basis.cpp`
- This corresponds to the Nektar++ file
  `library/LibUtilities/Foundations/Basis.cpp` and the classes defined therein.

# Basic class wrapping

As a very basic example of wrapping a class, let's consider the `SessionReader`
wrapper.

```c++
void export_SessionReader()
{
    py::class_<SessionReader,
           std::shared_ptr<SessionReader>,
           boost::noncopyable>(
               "SessionReader", py::no_init)

        .def("CreateInstance", SessionReader_CreateInstance)
        .staticmethod("CreateInstance")

        .def("GetSessionName", &SessionReader::GetSessionName,
             py::return_value_policy<py::copy_const_reference>())

        .def("Finalise", &SessionReader::Finalise)
        ;
}
```

## `py::class_<>`

This `Boost.Python` object defines a Python class in C++. It is templated, and
in this case we have the following template arguments:

- `SessionReader` is the class that will be wrapped
- `std::shared_ptr<SessionReader>` indicates that this object should be stored
  inside a shared (or smart) pointer, which we frequently use throughout the
  library, as can be seen by the frequent use of `SessionReaderSharedPtr`
- `boost::noncopyable` indicates that `Boost.Python` shouldn't try to
  automatically wrap the copy constructor of `SessionReader`. We add this here
  because of compiler errors due to subclasses used inside `SessionReader`, but
  generally, this should be used for abstract classes which can't be copied.

We then have two arguments to the:
- `"SessionReader"` is the name of the class in Python.
- `py::no_init` indicates this object has no publically-accessible
  initialiser. This is because for `SessionReader`, we define a factory-type
  function called `CreateInstance` instead.

## Wrapping member functions

We then call the `.def` function on the `class_<>`, which allows us to define
member functions on our class. This is equivalent to `def`-ing a function in
Python. `def` has two required parameters, and one optional parameter:

- The function name as a string, e.g. `"GetSessionName"`
- A function pointer that defines the C++ function that will be called
- An optional return policy, which we need to define when the C++ function
  returns a reference.

`Boost.Python` is very smart and can convert many Python objects to their
equivalent C++ function arguments, and C++ return types of the function to their
respective Python object. Many times therefore, one only needs to define the
`.def()` call.

However, there are some instances where we need to do some additional
conversion, mask some C++ complexity from the Python interface, or deal with
functions that return references. We describe ways to deal with this below.

### Thin wrappers

Instead of defining a function pointer to a member of the C++ class, we can
define a function pointer to a separate function that defines some extra
functionality. This is called a **thin wrapper**.

As an example, consider the `CreateInstance` function. In C++ we pass this
function the command line arguments in the usual `argc`, `argv` format. In
Python, command line arguments are defined as a list of strings inside
`sys.argv`. However, `Boost.Python` does not know how to convert this list to
`argc, argv`, so we need some additional code.

```c++
SessionReaderSharedPtr SessionReader_CreateInstance(py::list &ns)
{
    // ... some code here that converts a Python list to the standard
    // c/c++ (int argc, char **argv) format for command line arguments.
    // Then use this to construct a SessionReader and return it.
    SessionReaderSharedPtr sr = SessionReader::CreateInstance(argc, argv);
    return sr;
}
```

In Python, we can then simply call 

```python
session = SessionReader.CreateInstance(sys.argv)
```

In NekPy, you should make 

### Dealing with references

When dealing with functions in C++ that return references, e.g.

```c++
const NekDouble &GetFactor()
```

we need to supply an additional argument to `.def()`, since Python immutable
types such as strings and integers cannot be passed by reference. For a full
list of options, consult the `Boost.Python` guide. However a good rule of thumb
is to use `copy_const_reference` as highlighted above, which will create a copy
of the const reference and return this.

### Dealing with `Array<OneD, >`

The `LibUtilities/SharedArray.cpp` file contains a number of functions that
allow for the automatic conversion of Nektar++ `Array<OneD, >` storage to and
from NumPy `ndarray` objects. This means that you can wrap functions that take
these as parameters and return arrays very easily. However bear in mind the
following caveats:

- NumPy arrays created from Array objects will share their memory, so that
  changing the C++ array changes the contents of the NumPy array. Likewise, C++
  arrays created from NumPy arrays will share memory. Reference counting and
  capsules are used to ensure that memory should be persistent whilst the arrays
  are in use, either within Python or C++.
- Many functions in Nektar++ return Arrays through argument parameters. In
  Python this is a very unnatural way to write functions. For example:
  ```python
  # This is good
  x, y, z = exp.GetCoords()
  # This is bad
  x, y, z = np.zeros(10), np.zeros(10), np.zeros(10)
  exp.GetCoords(x,y,z)
  ```
  Use thin wrappers to overcome this problem. For examples of how to do this,
  particularly in returning tuples, consult the `StdRegions/StdExpansion.cpp`
  wrapper.
- `TwoD` and `ThreeD` arrays are not presently supported.

## Inheritance

Nektar++ makes heavy use of inheritance, which can be translated to Python quite
easily using `Boost.Python`. For a good example of how to do this, you can
examine the `StdRegions` wrapper for `StdExpansion` and its elements such as
`StdQuadExp`. In a cut-down form, these look like the following:

```c++
void export_StdExpansion()
{
    py::class_<StdExpansion,
               std::shared_ptr<StdExpansion>,
               boost::noncopyable>(
                   "StdExpansion", py::no_init);
}
void export_StdQuadExp()
{
    py::class_<StdQuadExp, py::bases<StdExpansion>,
               std::shared_ptr<StdQuadExp> >(
                   "StdQuadExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
}
```

Note the following:
- `StdExpansion` is an abstract class, so it has no initialiser and is
  non-copyable.
- We use `py::bases<StdExpansion>` in the definition of `StdQuadExp` to define
  its parent class. This does not necessarily need to include the full hierarchy
  of C++ inheritance: in `StdRegions` the inheritance graph for `StdQuadExp`
  looks like
  ```
  StdExpansion -> StdExpansion2D -> StdQuadExp
  ```
  In the above wrapper, we omit the StdExpansion2D call entirely.
- `py::init<>` is used to show how to wrap a C++ constructor. This can accept
  any arguments for which you have either written explicit wrappers or
  `Boost.Python` already knows how to convert.

## Wrapping enums

Most Nektar++ enumerators come in the form:

```c++
enum MyEnum {
    eItemOne,
    eItemTwo,
    SIZE_MyEnum
};
static const char *MyEnumMap[] = {
    "ItemOne"
    "ItemTwo"
    "ItemThree"
};
```

To wrap this, you can use the `NEKPY_WRAP_ENUM` macro defined in
`NekPyConfig.hpp`, which in this case can be used as `NEKPY_WRAP_ENUM(MyEnum,
MyEnumMap)`. Note that if instead of `const char *` the map is defined as a
`const std::string`, you can use the `NEKPY_WRAP_ENUM_STRING` macro.
