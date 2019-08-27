///////////////////////////////////////////////////////////////////////////////
//
// File: LocalElements.cpp
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
// Description: Python wrapper for Expansion elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/HexExp.h>

using namespace Nektar;
using namespace Nektar::LocalRegions;

void export_LocalElements()
{
    py::class_<SegExp, py::bases<Expansion, StdRegions::StdSegExp>,
               std::shared_ptr<SegExp> >(
                   "SegExp", py::init<const LibUtilities::BasisKey&,
                   const SpatialDomains::SegGeomSharedPtr &>());
    py::class_<TriExp, py::bases<Expansion, StdRegions::StdTriExp>,
               std::shared_ptr<TriExp> >(
                   "TriExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::TriGeomSharedPtr &>());
    py::class_<QuadExp, py::bases<Expansion, StdRegions::StdQuadExp>,
               std::shared_ptr<QuadExp> >(
                   "QuadExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::QuadGeomSharedPtr &>());
    py::class_<TetExp, py::bases<Expansion, StdRegions::StdTetExp>,
               std::shared_ptr<TetExp> >(
                   "TetExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::TetGeomSharedPtr &>());
    py::class_<PrismExp, py::bases<Expansion, StdRegions::StdPrismExp>,
               std::shared_ptr<PrismExp> >(
                   "PrismExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::PrismGeomSharedPtr &>());
    py::class_<PyrExp, py::bases<Expansion, StdRegions::StdPyrExp>,
               std::shared_ptr<PyrExp> >(
                   "PyrExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::PyrGeomSharedPtr &>());
    py::class_<HexExp, py::bases<Expansion, StdRegions::StdHexExp>,
               std::shared_ptr<HexExp> >(
                   "HexExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::HexGeomSharedPtr &>());

    NEKPY_SHPTR_FIX(SegExp, Expansion);
    NEKPY_SHPTR_FIX(SegExp, StdRegions::StdSegExp);
    NEKPY_SHPTR_FIX(SegExp, StdRegions::StdExpansion);
    NEKPY_SHPTR_FIX(TriExp, Expansion);
    NEKPY_SHPTR_FIX(TriExp, StdRegions::StdTriExp);
    NEKPY_SHPTR_FIX(TriExp, StdRegions::StdExpansion);
    NEKPY_SHPTR_FIX(QuadExp, Expansion);
    NEKPY_SHPTR_FIX(QuadExp, StdRegions::StdQuadExp);
    NEKPY_SHPTR_FIX(QuadExp, StdRegions::StdExpansion);
    NEKPY_SHPTR_FIX(TetExp, Expansion);
    NEKPY_SHPTR_FIX(TetExp, StdRegions::StdTetExp);
    NEKPY_SHPTR_FIX(TetExp, StdRegions::StdExpansion);
    NEKPY_SHPTR_FIX(PrismExp, Expansion);
    NEKPY_SHPTR_FIX(PrismExp, StdRegions::StdPrismExp);
    NEKPY_SHPTR_FIX(PrismExp, StdRegions::StdExpansion);
    NEKPY_SHPTR_FIX(PyrExp, Expansion);
    NEKPY_SHPTR_FIX(PyrExp, StdRegions::StdPyrExp);
    NEKPY_SHPTR_FIX(PyrExp, StdRegions::StdExpansion);
    NEKPY_SHPTR_FIX(HexExp, Expansion);
    NEKPY_SHPTR_FIX(HexExp, StdRegions::StdHexExp);
    NEKPY_SHPTR_FIX(HexExp, StdRegions::StdExpansion);
}
