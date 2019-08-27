///////////////////////////////////////////////////////////////////////////////
//
// File: Points.cpp
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
// Description: Python wrapper for Points.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/PointsType.h>

#include <LibUtilities/Python/NekPyConfig.hpp>

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
    NEKPY_WRAP_ENUM_STRING_DOCS(PointsType, kPointsTypeStr, 
        "Characterise the type of points.\n"
        "\n"
        "Sample usage: PointsType.GaussLobattoLegendre\n"
        "\n"
        "Available point types:\n"
        "\tNoPointsType: No points type.\n"
        "\tGaussGaussLegendre: 1D Gauss-Gauss-Legendre quadrature points.\n"
        "\tGaussRadauMLegendre: 1D Gauss-Radau-Legendre quadrature points,\n"
            "\t\tpinned at x=-1.\n"
        "\tGaussRadauPLegendre: 1D Gauss-Radau-Legendre quadrature points,\n"
            "\t\tpinned at x=1.\n"
        "\tGaussLobattoLegendre: 1D Gauss-Lobatto-Legendre quadrature points.\n"
        "\tGaussGaussChebyshev: 1D Gauss-Gauss-Chebyshev quadrature points.\n"
        "\tGaussRadauMChebyshev: 1D Gauss-Radau-Chebyshev quadrature points,\n"
            "\t\tpinned at x=-1.\n"
        "\tGaussRadauPChebyshev: 1D Gauss-Radau-Chebyshev quadrature points,\n"
            "\t\tpinned at x=1.\n"
        "\tGaussLobattoChebyshev: 1D Gauss-Lobatto-Legendre quadrature points\n"
        "\tGaussRadauMAlpha0Beta1: Gauss Radau pinned at x=-1,\n"
            "\t\talpha = 0, beta = 1\n"
        "\tGaussRadauMAlpha0Beta2: Gauss Radau pinned at x=-1,\n"
            "\t\talpha = 0, beta = 2\n"
        "\tGaussRadauMAlpha1Beta0: Gauss Radau pinned at x=-1,\n"
            "\t\talpha = 1, beta = 0\n"
        "\tGaussRadauMAlpha2Beta0: Gauss Radau pinned at x=-1,\n"
            "\t\talpha = 2, beta = 0\n"
        "\tGaussKronrodLegendre: 1D Gauss-Kronrod-Legendre quadrature points.\n"
        "\tGaussRadauKronrodMLegendre: 1D Gauss-Radau-Kronrod-Legendre quadrature\n"
            "\t\tpoints, pinned at x=-1.\n"
        "\tGaussRadauKronrodMAlpha1Beta0: 1D Gauss-Radau-Kronrod-Legendre\n"
            "\t\tpinned at x=-1, alpha = 1, beta = 0\n"
        "\tGaussLobattoKronrodLegendre: 1D Lobatto Kronrod quadrature points.\n"
        "\tPolyEvenlySpaced: 1D Evenly-spaced points using Lagrange polynomial.\n"
        "\tFourierEvenlySpaced: 1D Evenly-spaced points using Fourier Fit.\n"
        "\tFourierSingleModeSpaced: 1D Non Evenly-spaced points for Single Mode\n"
            "\t\tanalysis.\n"
        "\tBoundaryLayerPoints: 1D power law distribution for boundary layer points.\n"
        "\tBoundaryLayerPointsRev: 1D power law distribution for boundary layer\n"
            "\t\tpoints.\n"
        "\tNodalTriElec: 2D Nodal Electrostatic Points on a Triangle.\n"
        "\tNodalTriFekete: 2D Nodal Fekete Points on a Triangle.\n"
        "\tNodalTriEvenlySpaced: 2D Evenly-spaced points on a Triangle.\n"
        "\tNodalTetEvenlySpaced: 3D Evenly-spaced points on a Tetrahedron.\n"
        "\tNodalTetElec: 3D Nodal Electrostatic Points on a Tetrahedron.\n"
        "\tNodalPrismEvenlySpaced: 3D Evenly-spaced points on a Prism.\n"
        "\tNodalPrismElec: 3D electrostatically spaced points on a Prism.\n"
        "\tNodalTriSPI: 2D Nodal Symmetric positive internal triangle\n"
            "\t\t(Whitherden, Vincent).\n"
        "\tNodalTetSPI: 3D Nodal Symmetric positive internal tet\n"
            "\t\t(Whitherden, Vincent).\n"
        "\tNodalPrismSPI: 3D prism SPI.\n"
        "\tNodalQuadElec: 2D GLL for quad.\n"
        "\tNodalHexElec: 3D GLL for hex");

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
        .def("GetW", &Points<double>::GetW,
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
