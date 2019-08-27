///////////////////////////////////////////////////////////////////////////////
//
// File MeshPartitionPtScotch.h
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
// Description: PtScotch partitioner interface
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHPARTITIONPTSCOTCH_H
#define NEKTAR_SPATIALDOMAINS_MESHPARTITIONPTSCOTCH_H

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SpatialDomains/MeshPartition.h>

namespace Nektar
{
namespace SpatialDomains
{

class MeshPartitionPtScotch : public MeshPartition
{
public:
    /// Creates an instance of this class
    static MeshPartitionSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr session,
        int                                        meshDim,
        std::map<int, MeshEntity>                  element,
        CompositeDescriptor                        compMap)
    {
        return MemoryManager<MeshPartitionPtScotch>::AllocateSharedPtr(
            session, meshDim, element, compMap);
    }

    /// Name of class
    static std::string className;
    static std::string cmdSwitch;

    MeshPartitionPtScotch(const LibUtilities::SessionReaderSharedPtr session,
                          int                                        meshDim,
                          std::map<int, MeshEntity>                  element,
                          CompositeDescriptor                        compMap);
    virtual ~MeshPartitionPtScotch();

private:
    virtual void PartitionGraphImpl(int &nVerts, int &nVertConds,
                                    Nektar::Array<Nektar::OneD, int> &xadj,
                                    Nektar::Array<Nektar::OneD, int> &adjcy,
                                    Nektar::Array<Nektar::OneD, int> &vertWgt,
                                    Nektar::Array<Nektar::OneD, int> &vertSize,
                                    Nektar::Array<Nektar::OneD, int> &edgeWgt,
                                    int &nparts, int &volume,
                                    Nektar::Array<Nektar::OneD, int> &part);
};
}
}

#endif
