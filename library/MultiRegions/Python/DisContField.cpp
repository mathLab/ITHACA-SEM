///////////////////////////////////////////////////////////////////////////////
//
// File: DisContField.cpp
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
// Description: Python wrapper for DisContField.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/DisContField3D.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::MultiRegions;

template<class T>
std::shared_ptr<T> CreateDisContField(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,
    const std::string &var,
    const bool setupDG)
{
    return std::make_shared<T>(session, graph, var, setupDG);
}

template<class T, class S>
void export_DisContField_Helper(const char *name)
{
    py::class_<T, py::bases<S>, std::shared_ptr<T> >(name, py::no_init)
        .def("__init__", py::make_constructor(
                 &CreateDisContField<T>,
                 py::default_call_policies(),
                 (py::arg("session"), py::arg("graph"), py::arg("var"),
                  py::arg("setupDG") = true)));

    NEKPY_SHPTR_FIX(T, ExpList);
    NEKPY_SHPTR_FIX(T, S);
}

void export_DisContField()
{
    export_DisContField_Helper<DisContField1D, ExpList1D>("DisContField1D");
    export_DisContField_Helper<DisContField2D, ExpList2D>("DisContField2D");
    export_DisContField_Helper<DisContField3D, ExpList3D>("DisContField3D");
}
