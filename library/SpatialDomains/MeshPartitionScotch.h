///////////////////////////////////////////////////////////////////////////////
//
// File MeshPartitionScotch.h
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
// Description: Scotch partitioner interface
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_MESHPARTITIONSCOTCH_H
#define NEKTAR_SPATIALDOMAINS_MESHPARTITIONSCOTCH_H

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SpatialDomains/MeshPartition.h>

#ifdef NEKTAR_USE_MPI
#include <ptscotch.h>
#else
#include <scotch.h>
#endif

namespace Nektar
{
namespace SpatialDomains
{
    class MeshPartitionScotch : public MeshPartition
    {
        public:
            /// Creates an instance of this class
            static MeshPartitionSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr session,
                int                                        meshDim,
                std::map<int, MeshEntity>                  element,
                CompositeDescriptor                        compMap)
            {
                return MemoryManager<MeshPartitionScotch>
                    ::AllocateSharedPtr(session, meshDim, element, compMap);
            }

            /// Name of class
            static std::string className;
            static std::string cmdSwitch;

            MeshPartitionScotch(
                const LibUtilities::SessionReaderSharedPtr session,
                int                                        meshDim,
                std::map<int, MeshEntity>                  element,
                CompositeDescriptor                        compMap);
            virtual ~MeshPartitionScotch();

        private:
            virtual void PartitionGraphImpl(
                    int&                              nVerts,
                    int&                              nVertConds,
                    Nektar::Array<Nektar::OneD, int>& xadj,
                    Nektar::Array<Nektar::OneD, int>& adjcy,
                    Nektar::Array<Nektar::OneD, int>& vertWgt,
                    Nektar::Array<Nektar::OneD, int>& vertSize,
                    Nektar::Array<Nektar::OneD, int>& edgeWgt,
                    int&                              nparts,
                    int&                              volume,
                    Nektar::Array<Nektar::OneD, int>& part);

            void PartGraphVKway(
                    const SCOTCH_Num * const    n,
                    const SCOTCH_Num * const    xadj,
                    const SCOTCH_Num * const    adjncy,
                    const SCOTCH_Num * const    vwgt,
                    const SCOTCH_Num * const    vsize,
                    const SCOTCH_Num * const    wgtflag,
                    const SCOTCH_Num * const    numflag,
                    const SCOTCH_Num * const    nparts,
                    SCOTCH_Num * const          volume,
                    SCOTCH_Num * const          part);

            int PartGraph2 (
                    const SCOTCH_Num * const    n,
                    const SCOTCH_Num * const    xadj,
                    const SCOTCH_Num * const    adjncy,
                    const SCOTCH_Num * const    vwgt,
                    const SCOTCH_Num * const    adjwgt,
                    const SCOTCH_Num * const    numflag,
                    const SCOTCH_Num * const    nparts,
                    SCOTCH_Num * const          part,
                    SCOTCH_Num                  flagval,
                    double                      kbalval);

    };

}
}

#endif
