///////////////////////////////////////////////////////////////////////////////
//
// File: NekMesh.cpp
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
// Description: Python wrapper for NekMesh classes.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMesh/Module/Log.hpp>

void export_Element();
void export_Mesh();
void export_Module();
void export_Node();

PyObject *NekMeshErrorType = nullptr;

using NekMeshError = Nektar::NekMesh::NekMeshError;

PyObject* CreateExceptionClass(const char* name,
                               PyObject* baseTypeObj = PyExc_Exception)
{
    std::string qualifiedName0 = std::string("NekPy.NekMesh.") + name;

    PyObject* typeObj = PyErr_NewException(
        const_cast<char*>(qualifiedName0.c_str()), baseTypeObj, 0);

    if (!typeObj)
    {
        py::throw_error_already_set();
    }

    py::scope().attr(name) = py::handle<>(py::borrowed(typeObj));
    return typeObj;
}

void TranslateNekMeshError(NekMeshError const &e)
{
    PyErr_SetString(NekMeshErrorType, e.what());
}

BOOST_PYTHON_MODULE(_NekMesh)
{
    np::initialize();

    // Register the NekMeshError exception.
    NekMeshErrorType = CreateExceptionClass("NekMeshError");
    py::register_exception_translator<NekMeshError>(&TranslateNekMeshError);

    export_Element();
    export_Mesh();
    export_Module();
    export_Node();
}
