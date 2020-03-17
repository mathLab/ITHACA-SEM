///////////////////////////////////////////////////////////////////////////////
//
// File NekPyConfig.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/version.hpp>

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
#if PY_MAJOR_VERSION == 2
#define NEKPY_WRAP_ENUM_STRING_DOCS(ENUMNAME,MAPNAME,DOCSTRING)    \
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
#else
#define NEKPY_WRAP_ENUM_STRING_DOCS(ENUMNAME,MAPNAME,DOCSTRING)    \
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
                             PyUnicode_FromString(DOCSTRING));     \
    }
#endif

// Boost 1.62 and earlier don't have native support for std::shared_ptr.
#if BOOST_VERSION < 106300
#define NEKPY_SHPTR_FIX(SOURCE,TARGET)                             \
    py::implicitly_convertible<std::shared_ptr<SOURCE>,            \
                               std::shared_ptr<TARGET>>();
#else
#define NEKPY_SHPTR_FIX(SOURCE,TARGET)
#endif

