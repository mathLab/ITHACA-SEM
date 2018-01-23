#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/PointsType.h>

#include <NekPyConfig.hpp>

using namespace Nektar::LibUtilities;

typedef std::shared_ptr<Points<double> > PointsSharedPtr;

PointsSharedPtr Points_Create(const PointsKey &pts)
{
    return PointsManager()[pts];
}

py::tuple Points_GetZW(PointsSharedPtr pts)
{
    return py::make_tuple(pts->GetZ(), pts->GetW());
}

/**
 * @brief Points exports.
 */
void export_Points()
{
    NEKPY_WRAP_ENUM_STRING(PointsType, kPointsTypeStr);

    py::class_<PointsKey>("PointsKey", py::init<const int, const PointsType&>())

        .def("GetNumPoints", &PointsKey::GetNumPoints)
        .def("GetPointsType", &PointsKey::GetPointsType)
        .def("GetPointsDim", &PointsKey::GetPointsDim)
        .def("GetTotNumPoints", &PointsKey::GetTotNumPoints)
        ;

    py::class_<Points<double>,
               std::shared_ptr<Points<double> >,
               boost::noncopyable>(
                   "Points", py::no_init)
        .def("Create", &Points_Create)
        .staticmethod("Create")

        .def("Initialise", &Points<double>::Initialize)
        .def("GetPointsDim", &Points<double>::GetPointsDim)
        .def("GetPointsType", &Points<double>::GetPointsType)
        .def("GetNumPoints", &Points<double>::GetNumPoints)
        .def("GetTotNumPoints", &Points<double>::GetTotNumPoints)

        .def("GetZ", &Points<double>::GetZ,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetW", &Points<double>::GetZ,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetZW", &Points_GetZW)
        ;
}
