///////////////////////////////////////////////////////////////////////////////
//
// File: ContField.cpp
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
// Description: Python wrapper for ContField.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::MultiRegions;

std::shared_ptr<ContField1D> CreateContField1D(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,
    const std::string &var)
{
    return std::make_shared<ContField1D>(session, graph, var);
}

std::shared_ptr<ContField2D> CreateContField2D(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,
    const std::string &var,
    const bool checkSingular)
{
    return std::make_shared<ContField2D>(
        session, graph, var, true, checkSingular);
}

std::shared_ptr<ContField3D> CreateContField3D(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,
    const std::string &var,
    const bool checkSingular)
{
    return std::make_shared<ContField3D>(session, graph, var, checkSingular);
}

void export_ContField()
{
    py::class_<ContField1D, py::bases<ExpList1D>, std::shared_ptr<ContField1D>>(
        "ContField1D", py::no_init)
        .def("__init__", py::make_constructor(
                 &CreateContField1D,
                 py::default_call_policies(),
                 (py::arg("session"), py::arg("graph"), py::arg("var"))));

    NEKPY_SHPTR_FIX(ContField1D, ExpList);
    NEKPY_SHPTR_FIX(ContField1D, ExpList1D);
    NEKPY_SHPTR_FIX(ContField1D, DisContField1D);

    py::class_<ContField2D, py::bases<ExpList2D>, std::shared_ptr<ContField2D>>(
        "ContField2D", py::no_init)
        .def("__init__", py::make_constructor(
                 &CreateContField2D,
                 py::default_call_policies(),
                 (py::arg("session"), py::arg("graph"), py::arg("var"),
                  py::arg("checkSingular") = true)));

    NEKPY_SHPTR_FIX(ContField2D, ExpList);
    NEKPY_SHPTR_FIX(ContField2D, ExpList2D);
    NEKPY_SHPTR_FIX(ContField2D, DisContField2D);

    py::class_<ContField3D, py::bases<ExpList3D>, std::shared_ptr<ContField3D>>(
        "ContField3D", py::no_init)
        .def("__init__", py::make_constructor(
                 &CreateContField3D,
                 py::default_call_policies(),
                 (py::arg("session"), py::arg("graph"), py::arg("var"),
                  py::arg("checkSingular") = true)));

    NEKPY_SHPTR_FIX(ContField3D, ExpList);
    NEKPY_SHPTR_FIX(ContField3D, ExpList3D);
    NEKPY_SHPTR_FIX(ContField3D, DisContField3D);
}
