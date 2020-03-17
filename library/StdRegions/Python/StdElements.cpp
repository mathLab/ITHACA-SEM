///////////////////////////////////////////////////////////////////////////////
//
// File: StdElements.cpp
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
// Description: Python wrapper for all Nektar++ elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdHexExp.h>

using namespace Nektar;
using namespace Nektar::StdRegions;

void export_StdElements()
{
    py::class_<StdPointExp, py::bases<StdExpansion>,
               std::shared_ptr<StdPointExp> >(
                   "StdPointExp", py::init<const LibUtilities::BasisKey&>());
    py::class_<StdSegExp, py::bases<StdExpansion>,
               std::shared_ptr<StdSegExp> >(
                   "StdSegExp", py::init<const LibUtilities::BasisKey&>());
    py::class_<StdQuadExp, py::bases<StdExpansion>,
               std::shared_ptr<StdQuadExp> >(
                   "StdQuadExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
    py::class_<StdTriExp, py::bases<StdExpansion>,
               std::shared_ptr<StdTriExp> >(
                   "StdTriExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
    py::class_<StdTetExp, py::bases<StdExpansion>,
               std::shared_ptr<StdTetExp> >(
                   "StdTetExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
    py::class_<StdPrismExp, py::bases<StdExpansion>,
               std::shared_ptr<StdPrismExp> >(
                   "StdPrismExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
    py::class_<StdPyrExp, py::bases<StdExpansion>,
               std::shared_ptr<StdPyrExp> >(
                   "StdPyrExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
    py::class_<StdHexExp, py::bases<StdExpansion>,
               std::shared_ptr<StdHexExp> >(
                   "StdHexExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());

    NEKPY_SHPTR_FIX(StdSegExp, StdExpansion);
    NEKPY_SHPTR_FIX(StdQuadExp, StdExpansion);
    NEKPY_SHPTR_FIX(StdTriExp, StdExpansion);
    NEKPY_SHPTR_FIX(StdTetExp, StdExpansion);
    NEKPY_SHPTR_FIX(StdPrismExp, StdExpansion);
    NEKPY_SHPTR_FIX(StdPyrExp, StdExpansion);
    NEKPY_SHPTR_FIX(StdHexExp, StdExpansion);
}

