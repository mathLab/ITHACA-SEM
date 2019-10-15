///////////////////////////////////////////////////////////////////////////////
//
// File: StdExpansion.cpp
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
// Description: Python wrapper for StdExpansion.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <StdRegions/StdExpansion.h>

using namespace Nektar;
using namespace Nektar::StdRegions;

Array<OneD, NekDouble> StdExpansion_FwdTrans(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->FwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> StdExpansion_BwdTrans(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetTotPoints());
    exp->BwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> StdExpansion_IProductWRTBase(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->BwdTrans(in, out);
    return out;
}

NekDouble StdExpansion_PhysEvaluate(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &physvals)
{
    return exp->PhysEvaluate(coords, physvals);
}

NekDouble StdExpansion_L2(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    return exp->L2(in);
}

NekDouble StdExpansion_L2_Error(
    StdExpansionSharedPtr exp,
    const Array<OneD, const NekDouble> &in,
    const Array<OneD, const NekDouble> &err)
{
    return exp->L2(in, err);
}

py::tuple StdExpansion_GetCoords(StdExpansionSharedPtr exp)
{
    int nPhys = exp->GetTotPoints();
    int coordim = exp->GetCoordim();

    std::vector<Array<OneD, NekDouble> > coords(coordim);
    for (int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(nPhys);
    }

    switch (coordim)
    {
        case 1:
            exp->GetCoords(coords[0]);
            return py::make_tuple(coords[0]);
            break;
        case 2:
            exp->GetCoords(coords[0], coords[1]);
            return py::make_tuple(coords[0], coords[1]);
            break;
        case 3:
            exp->GetCoords(coords[0], coords[1], coords[2]);
            return py::make_tuple(coords[0], coords[1], coords[2]);
            break;
    }

    return py::tuple();
}

py::tuple StdExpansion_PhysDeriv(StdExpansionSharedPtr exp,
                                 const Array<OneD, const NekDouble> &inarray)
{
    int nPhys = exp->GetTotPoints();
    int coordim = exp->GetCoordim();

    std::vector<Array<OneD, NekDouble> > derivs(coordim);
    for (int i = 0; i < coordim; ++i)
    {
        derivs[i] = Array<OneD, NekDouble>(nPhys);
    }

    switch (coordim)
    {
        case 1:
            exp->PhysDeriv(inarray, derivs[0]);
            return py::make_tuple(derivs[0]);
            break;
        case 2:
            exp->PhysDeriv(inarray, derivs[0], derivs[1]);
            return py::make_tuple(derivs[0], derivs[1]);
            break;
        case 3:
            exp->PhysDeriv(inarray, derivs[0], derivs[1], derivs[2]);
            return py::make_tuple(derivs[0], derivs[1], derivs[2]);
            break;
    }

    return py::tuple();
}

void export_StdExpansion()
{
    py::class_<StdExpansion,
               std::shared_ptr<StdExpansion>,
               boost::noncopyable>(
                   "StdExpansion", py::no_init)

        .def("GetNcoeffs", &StdExpansion::GetNcoeffs)
        .def("GetTotPoints", &StdExpansion::GetTotPoints)
        .def("GetBasisType", &StdExpansion::GetBasisType)
        .def("GetPointsType", &StdExpansion::GetPointsType)
        .def("GetNverts", &StdExpansion::GetNverts)
        .def("GetNedges", &StdExpansion::GetNedges)
        .def("GetNfaces", &StdExpansion::GetNfaces)
        .def("DetShapeType", &StdExpansion::DetShapeType)
        .def("GetShapeDimension", &StdExpansion::GetShapeDimension)
        .def("Integral", &StdExpansion::Integral)

        .def("GetBasis", &StdExpansion::GetBasis,
             py::return_value_policy<py::copy_const_reference>())

        .def("GenMatrix", &StdExpansion::GenMatrix)
        .def("GetStdMatrix", &StdExpansion::GetStdMatrix)

        .def("FwdTrans", &StdExpansion_FwdTrans)
        .def("BwdTrans", &StdExpansion_BwdTrans)
        .def("IProductWRTBase", &StdExpansion_IProductWRTBase)

        .def("PhysEvaluate", &StdExpansion_PhysEvaluate)
        .def("L2", &StdExpansion_L2)
        .def("L2", &StdExpansion_L2_Error)

        .def("GetCoords", &StdExpansion_GetCoords)

        .def("PhysDeriv", &StdExpansion_PhysDeriv)
        ;
}
