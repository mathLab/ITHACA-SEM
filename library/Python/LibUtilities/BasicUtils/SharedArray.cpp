#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <NekPyConfig.hpp>

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
        return py::incref(
            np::from_data(
                arr.data(), np::dtype::get_builtin<T>(),
                py::make_tuple(arr.num_elements()), py::make_tuple(sizeof(T)),
                capsule).ptr());
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
        new (storage) Array<OneD, T>(array.shape(0), (T *)array.get_data());
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
