#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <NekPyConfig.hpp>

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
                             (void *)new std::shared_ptr<NekMatrix<T, StandardMatrixTag>>(mat),
                             (PyCapsule_Destructor)&CapsuleDestructor<T, StandardMatrixTag>)));
#endif

        int nRows = mat->GetRows(), nCols = mat->GetColumns();
        //StorageType storage = mat->GetStorageType();

        // if (storage != eFULL)
        // {
        //     // Only support full storage types at the moment - ignore symmetric,
        //     // banded, etc.
        // }

        return py::incref(
            np::from_data(
                mat->GetRawPtr(), np::dtype::get_builtin<T>(),
                py::make_tuple(nRows, nCols),
                py::make_tuple(sizeof(T), nCols * sizeof(T)),
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
