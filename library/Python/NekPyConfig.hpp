#ifdef BOOST_HAS_NUMPY

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace py = boost::python;
namespace np = boost::python::numpy;

#else

#include <boost/python.hpp>
#include <boost/numpy.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#endif

#define SIZENAME(s) SIZE_##s
#define NEKPY_WRAP_ENUM(ENUMNAME,MAPNAME)                          \
    {                                                              \
        py::enum_<ENUMNAME> tmp(#ENUMNAME);                        \
        for (int a = 0; a < (int)SIZENAME(ENUMNAME); ++a)          \
        {                                                          \
            tmp.value(MAPNAME[a], (ENUMNAME)a);                    \
        }                                                          \
        tmp.export_values();                                       \
    }
#define NEKPY_WRAP_ENUM_STRING(ENUMNAME,MAPNAME)                   \
    {                                                              \
        py::enum_<ENUMNAME> tmp(#ENUMNAME);                        \
        for (int a = 0; a < (int)SIZENAME(ENUMNAME); ++a)          \
        {                                                          \
            tmp.value(MAPNAME[a].c_str(), (ENUMNAME)a);            \
        }                                                          \
        tmp.export_values();                                       \
    }
#define NEKPY_WRAP_ENUM_STRING_DOCS(ENUMNAME,MAPNAME, DOCSTRING)   \
    {                                                              \
        py::enum_<ENUMNAME> tmp(#ENUMNAME);                        \
        for (int a = 0; a < (int)SIZENAME(ENUMNAME); ++a)          \
        {                                                          \
            tmp.value(MAPNAME[a].c_str(), (ENUMNAME)a);            \
        }                                                          \
        tmp.export_values();                                       \
        PyTypeObject * pto =                                       \
            reinterpret_cast<PyTypeObject*>(tmp.ptr());            \
        PyDict_SetItemString(pto->tp_dict, "__doc__",              \
            PyString_FromString(DOCSTRING));                       \
    }

