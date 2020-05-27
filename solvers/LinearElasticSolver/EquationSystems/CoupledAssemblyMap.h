///////////////////////////////////////////////////////////////////////////////
//
// File CoupledAssemblyMap.h
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
// Description: Wrapper class around the library
// LocalToGlobalC0ContMap class for use in the Couplied Linearised NS
// solver.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COUPLEDASSEMBLYMAP_H
#define NEKTAR_SOLVERS_COUPLEDASSEMBLYMAP_H

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <SpatialDomains/MeshGraph.h>

namespace Nektar
{

/**
 * @brief Modified version of MultiRegions::AssemblyMapCG that allows for
 * coupled fields [u,v,w] instead of individual scalar fields u, v and w.
 */
class CoupledAssemblyMap : public MultiRegions::AssemblyMapCG
{
    typedef SpatialDomains::BoundaryConditionShPtr BoundaryCondShPtr;

public:
    CoupledAssemblyMap(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const SpatialDomains::MeshGraphSharedPtr          &graph,
        const MultiRegions::AssemblyMapCGSharedPtr        &cgMap,
        const Array<OneD, const BoundaryCondShPtr>        &boundaryConditions,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields);
};

typedef std::shared_ptr<CoupledAssemblyMap> CoupledAssemblyMapSharedPtr;

}

#endif
