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
        "\tpointsType (PointsType object): The type of quadrature points.",
        py::init<const int, const PointsType&>())

        .def("GetNumPoints", &PointsKey::GetNumPoints,
            "Get the number of quadrature points in PointsKey.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tInteger defining the number of quadrature points in PointsKey.")
        .def("GetPointsType", &PointsKey::GetPointsType,
            "Get the type of points in PointsKey.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tPointsType object specifying the type of points in PointsKey.")
        .def("GetPointsDim", &PointsKey::GetPointsDim,
            "Get the dimension of the PointsKey.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tInteger characterising the dimension of the PointsKey (e.g. 2\n"
            "\tfor 2D PointsKey).")
        .def("GetTotNumPoints", &PointsKey::GetTotNumPoints,
            "Get the total number of points in PointsKey.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tInteger defining the total number of points in PointsKey.")
        ;

    py::class_<Points<double>, 
               std::shared_ptr<Points<double> >, 
               boost::noncopyable>(
                   "Points", 
                   "Create a set of points which can be used to calculate\n"
                   "quadrature zeros, weights etc."
                   "\n"
                   "Args:\n"
                   "\tNone",
                   py::no_init)
        .def("Create", &Points_Create,
            "Create a Points object using PointsKey.\n"
            "\n"
            "Args:\n"
            "\tpointsKey (PointsKey object): The PointsKey to be used to\n"
            "\tcreate points.\n"
            "Returns:\n"
            "\tPoints object created with the given PointsKey.")
        .staticmethod("Create")

        .def("Initialise", &Points<double>::Initialize,
            "Initialise Points object by calculating points, weights\n"
            "and differentiation matrix.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tNone")
        .def("GetPointsDim", &Points<double>::GetPointsDim,
            "Get the dimension of the Points object.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tInteger characterising the dimension of the Points object\n"
            "\t(e.g. 2 for 2D Points object).")
        .def("GetPointsType", &Points<double>::GetPointsType,
            "Get the type of points in Points object.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tPointsType object specifying the type of points in Points\n"
            "\tobject.")
        .def("GetNumPoints", &Points<double>::GetNumPoints,
            "Get the number of quadrature points in Points object.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tInteger defining the number of quadrature points in Points\n"
            "\tobject.")
        .def("GetTotNumPoints", &Points<double>::GetTotNumPoints,
            "Get the total number of points in Points object.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tInteger defining the total number of points in Points\n"
            "\tobject.")

        .def("GetZ", &Points<double>::GetZ,
             py::return_value_policy<py::copy_const_reference>(), 
             "Get quadrature zeros.\n"
             "\n"
             "Args:\n"
             "\tNone\n"
             "Returns:\n"
             "\tNumPy ndarray of length equal to the number of quadrature\n"
             "\tpoints, containing quadrature zeros.")
        .def("GetW", &Points<double>::GetZ,
             py::return_value_policy<py::copy_const_reference>(),
             "Get quadrature weights.\n"
             "\n"
             "Args:\n"
             "\tNone\n"
             "Returns:\n"
             "\tNumPy ndarray of length equal to the number of quadrature\n"
             "\tpoints, containing quadrature weights.")
        .def("GetZW", &Points_GetZW,
            "Get quadrature zeros and weights.\n"
            "\n"
            "Args:\n"
            "\tNone\n"
            "Returns:\n"
            "\tTuple containing the quadrature zeros and quadrature weights,\n"
            "\ti.e. (Points.GetZ(), Points.GetW()).")
        .def("GetD", &Points_GetD,
            "Get the differentiation matrix.\n"
             "\n"
             "Args:\n"
             "\tNone\n"
             "Returns:\n"
             "\tNumPy ndarray of n x n dimensions, where n is equal to\n"
             "\tthe number of quadrature points, containing the\n"
             "\tdifferentiation matrix.")
        .def("GetD", &Points_GetD2,
            "Get the differentiation matrix.\n"
             "\n"
             "Args:\n"
             "\tdir (Direction object): The direction of the desired\n"
             "\tdifferentiation matrix.\n"
             "Returns:\n"
             "\tNumPy ndarray of n x n dimensions, where n is equal to\n"
             "\tthe number of quadrature points containing the\n"
             "\tdifferentiation matrix.")
        ;
}
