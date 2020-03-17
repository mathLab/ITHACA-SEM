///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrix.cpp
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
// Description: Python wrapper for NekMatrix.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;

#if PY_MAJOR_VERSION == 2
template<typename T, typename F>
void NekMatrixCapsuleDestructor(void *ptr)
{
    std::shared_ptr<NekMatrix<T, F>> *mat =
        (std::shared_ptr<NekMatrix<T, F>> *)ptr;
    delete mat;
}
#else
template<typename T, typename F>
void NekMatrixCapsuleDestructor(PyObject *ptr)
{
    std::shared_ptr<NekMatrix<T, F>> *mat =
        (std::shared_ptr<NekMatrix<T, F>> *)PyCapsule_GetPointer(ptr, 0);
    delete mat;
}
#endif

template<typename T>
struct NekMatrixToPython
{
    static PyObject *convert(
        std::shared_ptr<NekMatrix<T, StandardMatrixTag>> const &mat)
    {
        // Create a Python capsule to hold a pointer that contains a lightweight
        // copy of arr. That way we guarantee Python will still have access to
        // the memory allocated inside arr even if arr is deallocated in C++.
#if PY_MAJOR_VERSION == 2
        py::object capsule(
            py::handle<>(PyCObject_FromVoidPtr(
                             new std::shared_ptr<NekMatrix<T, StandardMatrixTag>>(mat),
                             NekMatrixCapsuleDestructor<T, StandardMatrixTag>)));
#else
        py::object capsule(
            py::handle<>(PyCapsule_New(
                             (void *)new std::shared_ptr<NekMatrix<T, StandardMatrixTag>>(mat), NULL,
                             (PyCapsule_Destructor)&NekMatrixCapsuleDestructor<T, StandardMatrixTag>)));
#endif

        int nRows = mat->GetRows(), nCols = mat->GetColumns();
        MatrixStorage storage = mat->GetStorageType();

        ASSERTL0(storage == eFULL, "Only full storage matrices are currently supported.");

        return py::incref(
            np::from_data(
                mat->GetRawPtr(), np::dtype::get_builtin<T>(),
                py::make_tuple(nRows, nCols),
                py::make_tuple(sizeof(T), nRows * sizeof(T)),
                capsule).ptr());
    }
};

template<typename T>
void export_NekMatrix()
{
    py::to_python_converter<std::shared_ptr<NekMatrix<T, StandardMatrixTag>>,
                            NekMatrixToPython<T>>();
}

template void export_NekMatrix<double>();
