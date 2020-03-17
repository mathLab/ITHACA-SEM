///////////////////////////////////////////////////////////////////////////////
//
// File: ExpList2D.cpp
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
// Description: Python wrapper for ExpList2D.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::MultiRegions;

void export_ExpList2D()
{
    py::class_<ExpList1D, py::bases<ExpList>,
               std::shared_ptr<ExpList1D> >(
                   "ExpList1D", py::init<
                   const LibUtilities::SessionReaderSharedPtr &,
                   const SpatialDomains::MeshGraphSharedPtr &>());
    py::class_<ExpList2D, py::bases<ExpList>,
               std::shared_ptr<ExpList2D> >(
                   "ExpList2D", py::init<
                   const LibUtilities::SessionReaderSharedPtr &,
                   const SpatialDomains::MeshGraphSharedPtr &>());
    py::class_<ExpList3D, py::bases<ExpList>,
               std::shared_ptr<ExpList3D> >(
                   "ExpList3D", py::init<
                   const LibUtilities::SessionReaderSharedPtr &,
                   const SpatialDomains::MeshGraphSharedPtr &>());

    NEKPY_SHPTR_FIX(ExpList1D, ExpList);
    NEKPY_SHPTR_FIX(ExpList2D, ExpList);
    NEKPY_SHPTR_FIX(ExpList3D, ExpList);
}
