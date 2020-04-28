///////////////////////////////////////////////////////////////////////////////
//
// File: ExpList.cpp
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
// Description: Python wrapper for ExpList.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList.h>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <fstream>
#include <sstream>
#include <string>

using namespace Nektar;
using namespace Nektar::StdRegions;
using namespace Nektar::LocalRegions;
using namespace Nektar::SpatialDomains;
using namespace Nektar::MultiRegions;

ExpansionSharedPtr ExpList_GetExp(ExpListSharedPtr exp, int i)
{
    return exp->GetExp(i);
}

void ExpList_WriteVTK(ExpListSharedPtr exp, std::string filename)
{
    std::ofstream out(filename.c_str());
    exp->WriteVtkHeader(out);
    for (int i = 0; i < exp->GetExpSize(); ++i)
    {
        exp->WriteVtkPieceHeader(out, i);
        exp->WriteVtkPieceFooter(out, i);
    }
    exp->WriteVtkFooter(out);
}

Array<OneD, NekDouble> ExpList_FwdTrans(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->FwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> ExpList_BwdTrans(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNpoints());
    exp->BwdTrans(in, out);
    return out;
}

Array<OneD, NekDouble> ExpList_IProductWRTBase(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs());
    exp->IProductWRTBase(in, out);
    return out;
}

Array<OneD, NekDouble> ExpList_MultiplyByInvMassMatrix(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs(), 0.0);
    exp->MultiplyByInvMassMatrix(in, out);
    return out;
}

Array<OneD, NekDouble> ExpList_HelmSolve(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in,
    const py::object constFactorMap,
    const py::object varCoeffMap)
{
    Array<OneD, NekDouble> out(exp->GetNcoeffs(), 0.0);

    StdRegions::ConstFactorMap facMap = StdRegions::NullConstFactorMap;
    StdRegions::VarCoeffMap coeffMap = StdRegions::NullVarCoeffMap;

    if (!constFactorMap.is_none())
    {
        facMap = py::extract<StdRegions::ConstFactorMap>(constFactorMap);
    }
    if (!varCoeffMap.is_none())
    {
        coeffMap = py::extract<StdRegions::VarCoeffMap>(varCoeffMap);
    }

    exp->HelmSolve(in, out, facMap, coeffMap);
    return out;
}

NekDouble ExpList_L2(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    return exp->L2(in);
}

NekDouble ExpList_L2_Error(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in,
    const Array<OneD, const NekDouble> &err)
{
    return exp->L2(in, err);
}

NekDouble ExpList_Linf(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in)
{
    return exp->Linf(in);
}

NekDouble ExpList_Linf_Error(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &in,
    const Array<OneD, const NekDouble> &err)
{
    return exp->Linf(in, err);
}

py::tuple ExpList_GetCoords(ExpListSharedPtr exp)
{
    int nPhys = exp->GetNpoints();
    int coordim = exp->GetCoordim(0);

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

void ExpList_SetPhysArray(
    ExpListSharedPtr exp,
    Array<OneD, NekDouble> inarray)
{
    exp->SetPhysArray(inarray);
}

void ExpList_SetPhys(
    ExpListSharedPtr exp,
    const Array<OneD, const NekDouble> &inarray)
{
    exp->SetPhys(inarray);
}

const Array<OneD, const NekDouble> ExpList_GetPhys(ExpListSharedPtr exp)
{
    return exp->GetPhys();
}

NekDouble ExpList_PhysIntegral(ExpListSharedPtr exp)
{
    return exp->PhysIntegral();
}

std::string ExpList_GetPhysAddress(ExpListSharedPtr exp)
{
    std::stringstream ss;
    ss << static_cast<const void*>(&(exp->GetPhys()[0]));
    return ss.str();
}

void ExpList_ResetManagers(ExpListSharedPtr exp)
{
    exp->ClearGlobalLinSysManager();
    LibUtilities::NekManager<
        LocalRegions::MatrixKey, DNekScalMat,
        LocalRegions::MatrixKey::opLess>::ClearManager();
    LibUtilities::NekManager<
        LocalRegions::MatrixKey, DNekScalBlkMat,
        LocalRegions::MatrixKey::opLess>::ClearManager();
}

void export_ExpList()
{
    int (ExpList::*GetNcoeffs)() const = &ExpList::GetNcoeffs;

    py::class_<ExpList,
               std::shared_ptr<ExpList>,
               boost::noncopyable>(
                   "ExpList", py::no_init)

        // Expansions
        .def("GetExp", ExpList_GetExp)
        .def("GetExpSize", &ExpList::GetExpSize)

        // Query points and offset information
        .def("GetNpoints", &ExpList::GetNpoints)
        .def("GetNcoeffs", GetNcoeffs)
        .def("GetCoords", &ExpList_GetCoords)
        .def("GetPhys_Offset", &ExpList::GetPhys_Offset)
        .def("GetCoeff_Offset", &ExpList::GetCoeff_Offset)

        // Evaluations
        .def("PhysEvaluate", &ExpList::PhysEvaluate)

        // Operators
        .def("FwdTrans", &ExpList_FwdTrans)
        .def("BwdTrans", &ExpList_BwdTrans)
        .def("IProductWRTBase", &ExpList_IProductWRTBase)
        .def("MultiplyByInvMassMatrix", &ExpList_MultiplyByInvMassMatrix)
        .def("HelmSolve", &ExpList_HelmSolve, (
                 py::arg("in"),
                 py::arg("constFactorMap") = py::object(),
                 py::arg("varCoeffMap") = py::object()
                 ))

        // Error norms
        .def("L2", &ExpList_L2)
        .def("L2", &ExpList_L2_Error)
        .def("Linf", &ExpList_Linf)
        .def("Linf", &ExpList_Linf_Error)

        // Storage setups
        .def("SetPhysArray", &ExpList_SetPhysArray)
        .def("SetPhys", &ExpList_SetPhys)
        .def("GetPhys", &ExpList_GetPhys)
        .def("SetPhysState", &ExpList::SetPhysState)
        .def("GetPhysState", &ExpList::GetPhysState)
        .def("PhysIntegral", &ExpList_PhysIntegral)
        .def("GetPhysAddress", &ExpList_GetPhysAddress)

        // Misc functions
        .def("WriteVTK", &ExpList_WriteVTK)
        .def("ResetManagers", &ExpList_ResetManagers)
        ;
}
