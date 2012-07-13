///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG3D.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: C0-continuous assembly mappings specific to 3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLYMAPCG3D_H
#define MULTIREGIONS_ASSEMBLYMAPCG3D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/ExpList.h>
//#include <MultiRegions/ExpList1D.h>
//#include <MultiRegions/ExpList2D.h>


namespace Nektar
{
    namespace MultiRegions
    {
        const static map<int,pair<int, StdRegions::Orientation> > NullIntPairMap;

        class AssemblyMapCG3D;
        typedef boost::shared_ptr<AssemblyMapCG3D>  AssemblyMapCG3DSharedPtr;

        /// Constructs mappings for the C0 scalar continuous Galerkin formulation.
        class AssemblyMapCG3D: public AssemblyMapCG
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMapCG3D(
                const LibUtilities::SessionReaderSharedPtr &pSession);

            /// Constructor for the 3D expansion mappings with boundary
            /// conditions.
            MULTI_REGIONS_EXPORT AssemblyMapCG3D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                const map<int,int>& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const map<int,pair<int, StdRegions::Orientation> >& periodicFacesId);
            

            /// General constructor for expansions of all dimensions without
            /// boundary conditions.
            MULTI_REGIONS_EXPORT AssemblyMapCG3D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp);
            
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMapCG3D();
            
        private:
            /// Construct mappings for a three-dimensional scalar expansion.
            void SetUp3DExpansionC0ContMap(
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp = NullExpListSharedPtrArray,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                SpatialDomains::NullBoundaryConditionShPtrArray,
                const map<int,int>& periodicVerticesId = NullIntIntMap,
                const map<int,int>& periodicEdgesId = NullIntIntMap,
                const map<int,pair<int, StdRegions::Orientation> >& periodicFacesId = NullIntPairMap);
        }; // class
    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLYMAPCG3D_H


