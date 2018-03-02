#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

#include <NekPyConfig.hpp>

using namespace Nektar::LibUtilities;

BasisSharedPtr Basis_Create(const BasisKey &pts)
{
    return BasisManager()[pts];
}

py::tuple Basis_GetZW(BasisSharedPtr pts)
{
    return py::make_tuple(pts->GetZ(), pts->GetW());
}

/**
 * @brief Basis exports.
 */
void export_Basis()
{
    // Enumerator for basis type
    NEKPY_WRAP_ENUM(BasisType, BasisTypeMap);

    py::class_<BasisKey>("BasisKey", py::init<const BasisType&, const int,
                         const PointsKey&>())

        .def("GetNumModes", &BasisKey::GetNumModes)
        .def("GetTotNumModes", &BasisKey::GetTotNumModes)
        .def("GetNumPoints", &BasisKey::GetNumPoints)
        .def("GetTotNumPoints", &BasisKey::GetTotNumPoints)
        .def("GetBasisType", &BasisKey::GetBasisType)
        .def("GetPointsKey", &BasisKey::GetPointsKey)
        .def("GetPointsType", &BasisKey::GetPointsType)
        .def("Collocation", &BasisKey::Collocation)
        ;

    py::class_<Basis,
               std::shared_ptr<Basis> >(
               "Basis", py::no_init)

        .def("Create", &Basis_Create)
        .staticmethod("Create")

        .def("GetNumModes", &Basis::GetNumModes)
        .def("GetTotNumModes", &Basis::GetTotNumModes)
        .def("GetNumPoints", &Basis::GetNumPoints)
        .def("GetTotNumPoints", &Basis::GetTotNumPoints)
        .def("GetBasisType", &Basis::GetBasisType)
        .def("GetPointsKey", &Basis::GetPointsKey)
        .def("GetBasisKey", &Basis::GetBasisKey)
        .def("GetPointsType", &Basis::GetBasisType)
        .def("Initialize", &Basis::Initialize)

        .def("GetZ", &Basis::GetZ,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetW", &Basis::GetZ,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetZW", &Basis_GetZW)

        .def("GetBdata", &Basis::GetBdata,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetDbdata", &Basis::GetDbdata,
             py::return_value_policy<py::copy_const_reference>())
        ;
}
