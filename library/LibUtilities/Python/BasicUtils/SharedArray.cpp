///////////////////////////////////////////////////////////////////////////////
//
// File: SharedArray.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
//
// License for the specific language governing rights and limitations under
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
// Description: Python wrapper for ShareArray.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <type_traits>

using namespace Nektar;
using namespace Nektar::LibUtilities;

#if PY_MAJOR_VERSION == 2
template<typename T>
void CapsuleDestructor(void *ptr)
{
    Array<OneD, T> *tmp = (Array<OneD, T> *)ptr;
    delete tmp;
}
#else
template<typename T>
void CapsuleDestructor(PyObject *ptr)
{
    Array<OneD, T> *tmp = (Array<OneD, T> *)PyCapsule_GetPointer(ptr, 0);
    delete tmp;
}
#endif

template<typename T>
struct OneDArrayToPython
{
    static PyObject *convert(Array<OneD, T> const &arr)
    {
        // Create a Python capsule to hold a pointer that contains a lightweight
        // copy of arr. That way we guarantee Python will still have access to
        // the memory allocated inside arr even if arr is deallocated in C++.
#if PY_MAJOR_VERSION == 2
        py::object capsule(
            py::handle<>(PyCObject_FromVoidPtr(
                             new Array<OneD, T>(arr), CapsuleDestructor<T>)));
#else
        py::object capsule(
            py::handle<>(PyCapsule_New(
                             new Array<OneD, T>(arr), 0,
                             (PyCapsule_Destructor)&CapsuleDestructor<T>)));
#endif
        PyObject *tmp = py::incref(
            np::from_data(
                arr.data(), np::dtype::get_builtin<T>(),
                py::make_tuple(arr.num_elements()), py::make_tuple(sizeof(T)),
                capsule).ptr());

        return tmp;
    }
};

template <typename T>
struct PythonToOneDArray
{
    PythonToOneDArray()
    {
        py::converter::registry::push_back(
            &convertible, &construct, py::type_id<Array<OneD, T> >());
    }
    static void *convertible(PyObject *objPtr)
    {
        try
        {
            py::object obj((py::handle<>(py::borrowed(objPtr))));
            np::ndarray array = py::extract<np::ndarray>(obj);

            // Check data types match
            np::dtype dtype = np::dtype::get_builtin<
                typename boost::remove_const<T>::type>();
            if (dtype != array.get_dtype())
            {
                return 0;
            }

            // Check shape is 1D
            if (array.get_nd() != 1)
            {
                return 0;
            }
        }
        catch (boost::python::error_already_set)
        {
            py::handle_exception();
            PyErr_Clear();
            return 0;
        }

        return objPtr;
    }

    static void decrement(void *objPtr) 
    {
        //PyObject *pyObjPtr = (PyObject *)objPtr;
        //Py_DECREF(pyObjPtr);
        py::decref((PyObject *)objPtr);
    }

    static void construct(
        PyObject *objPtr,
        py::converter::rvalue_from_python_stage1_data* data)
    {
        // This has to be a _borrowed_ reference, otherwise at the end of this
        // scope it seems memory gets deallocated
        py::object obj((py::handle<>(py::borrowed(objPtr))));
        np::ndarray array = py::extract<np::ndarray>(obj);
        
        // Construct OneD array from numpy array
        void *storage = (
            (py::converter::rvalue_from_python_storage<Array<OneD, T> >*)data)
            ->storage.bytes;
        data->convertible = storage;
        void *memory_pointer = objPtr;
        using nonconst_t = typename std::remove_const<T>::type;
        new (storage) Array<OneD, T>(array.shape(0), (nonconst_t *)array.get_data(), memory_pointer, &decrement);
        py::incref(objPtr);
    }

};

template<typename T>
void export_SharedArray()
{
    py::to_python_converter<Array<OneD, const T>, OneDArrayToPython<const T> >();
    py::to_python_converter<Array<OneD, T>, OneDArrayToPython<T> >();

    PythonToOneDArray<const T>();
    PythonToOneDArray<T>();
}

template void export_SharedArray<double>();
