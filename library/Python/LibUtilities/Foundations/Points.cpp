#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/PointsType.h>

#include <NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;

typedef std::shared_ptr<Points<double> > PointsSharedPtr;
typedef std::shared_ptr<NekMatrix<NekDouble> > MatrixSharedPtrType;

PointsSharedPtr Points_Create(const PointsKey &pts)
{
    return PointsManager()[pts];
}

py::tuple Points_GetZW(PointsSharedPtr pts)
{
    return py::make_tuple(pts->GetZ(), pts->GetW());
}

MatrixSharedPtrType Points_GetD(PointsSharedPtr pts) 
{
    return pts->GetD();
}

MatrixSharedPtrType Points_GetD2(PointsSharedPtr pts, 
    Direction dir)
{
    return pts->GetD(dir);
}


/**
 * @brief Points exports.
 */
void export_Points()
{
    NEKPY_WRAP_ENUM_STRING(PointsType, kPointsTypeStr);

    py::class_<PointsKey>("PointsKey", 
        "Create a PointsKey which uniquely defines quadrature points.\n"
        "\n"
        "Args:\n"
        "\tnQuadPoints (integer): The number of quadrature points.\n"
        "\tpointsType (PointsType object): The type of quadrature points.\n"
        "Returns:\n"
        "\tPointsKey object defining quadrature points.",
        py::init<const int, const PointsType&>())

        .def("GetNumPoints", &PointsKey::GetNumPoints)
        .def("GetPointsType", &PointsKey::GetPointsType)
        .def("GetPointsDim", &PointsKey::GetPointsDim)
        .def("GetTotNumPoints", &PointsKey::GetTotNumPoints)
        ;

    py::class_<Points<double>, 
               std::shared_ptr<Points<double> >, 
               boost::noncopyable>(
                   "Points", 
                   "docstring", 
                   py::no_init)
        .def("Create", &Points_Create)
        .staticmethod("Create")

        .def("Initialise", &Points<double>::Initialize)
        .def("GetPointsDim", &Points<double>::GetPointsDim)
        .def("GetPointsType", &Points<double>::GetPointsType)
        .def("GetNumPoints", &Points<double>::GetNumPoints)
        .def("GetTotNumPoints", &Points<double>::GetTotNumPoints)

        .def("GetZ", &Points<double>::GetZ,
             py::return_value_policy<py::copy_const_reference>(), 
             "Get quadrature zeros.\n"
             "\n"
             "Returns:\n"
             "\tNumPy ndarray of length equal to the number of quadrature\n"
             "\tpoints, containing quadrature zeros.")
        .def("GetW", &Points<double>::GetZ,
             py::return_value_policy<py::copy_const_reference>(),
             "Get quadrature weights.\n"
             "\n"
             "Returns:\n"
             "\tNumPy ndarray of length equal to the number of quadrature\n"
             "\tpoints, containing quadrature weights.")
        .def("GetZW", &Points_GetZW,
            "Get quadrature zeros and weights.\n"
            "\n"
            "Returns:\n"
            "\tTuple containing the quadrature zeros and quadrature weights,\n"
            "\ti.e. (Points.GetZ(), Points.GetW()).")
        .def("GetD", &Points_GetD,
            "Get the differentiation matrix.\n"
             "\n"
             "Returns:\n"
             "\tNumPy ndarray of n x n dimensions, where n is equal to\n"
             "\tthe number of quadrature points, containing the\n"
             "\tdifferentiation matrix.")
        // to-do: fill in Args
        .def("GetD", &Points_GetD2,
            "Get the differentiation matrix.\n"
             "\n"
             "Args:\n"
             "\t WIP\n"
             "Returns:\n"
             "\tNumPy ndarray of n x n dimensions, where n is equal to\n"
             "\tthe number of quadrature points containing the\n"
             "\tdifferentiation matrix.")
        ;
}
