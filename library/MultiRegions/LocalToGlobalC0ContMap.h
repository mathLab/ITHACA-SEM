///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalC0ContMap.h
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
// Description: C0-continuous Local to Global mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCALTOGLOBALC0CONTMAP_H
#define MULTIREGIONS_LOCALTOGLOBALC0CONTMAP_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalBaseMap.h>
#include <MultiRegions/ExpList0D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>

#include <LocalRegions/PointExp.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {


        static map<int,int> NullIntIntMap;
        const static vector<map<int,int> > NullVecIntIntMap;

        /// Constructs mappings for the C0 scalar continuous Galerkin formulation.
        class LocalToGlobalC0ContMap: public LocalToGlobalBaseMap
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT LocalToGlobalC0ContMap(
                                   const LibUtilities::CommSharedPtr &pComm);

            /// General constructor for expansions of all dimensions without
            /// boundary conditions.
            MULTI_REGIONS_EXPORT LocalToGlobalC0ContMap(
                                   const LibUtilities::CommSharedPtr &pComm,
                                   const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType);

            /// Constructor for the 1D expansion mappings with boundary
            /// conditions.
            MULTI_REGIONS_EXPORT LocalToGlobalC0ContMap(
                                   const LibUtilities::CommSharedPtr &pComm,
                                   const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType,
                                   const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId);

            /// Constructor for the 2D expansion mappings with boundary
            /// conditions.
            MULTI_REGIONS_EXPORT LocalToGlobalC0ContMap(
                                   const LibUtilities::CommSharedPtr &pComm,
                                   const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType,
                                   const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const vector<map<int,int> >& periodicVerticesId,
                                   const map<int,int>& periodicEdgesId,
                                   const bool checkIfSystemSingular);



            /// Constructor for the 3D expansion mappings with boundary
            /// conditions.
            MULTI_REGIONS_EXPORT LocalToGlobalC0ContMap(
                                   const LibUtilities::CommSharedPtr &pComm,
                                   const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType,
                                   const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId,
                                   const map<int,int>& periodicEdgesId,
                                   const map<int,int>& periodicFacesId);

            /// Destructor.
            MULTI_REGIONS_EXPORT ~LocalToGlobalC0ContMap();


            /** Construct optimal ordering a two-dimensional expansion
            /*  given a vector of boundary condition information
            */
            MULTI_REGIONS_EXPORT int SetUp2DGraphC0ContMap(
                                       const ExpList  &locExp,
                                       const GlobalSysSolnType solnType,
                                       const Array<OneD, const MultiRegions::ExpListSharedPtr>  &bndCondExp,
                                       const Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >
                                       &bndConditions,
                                       const vector<map<int,int> >& periodicVerticesId,
                                       const map<int,int>& periodicEdgesId,
                                       Array<OneD, map<int,int> > &Dofs,
                                       Array<OneD, map<int,int> > &ReorderedGraphVertId,
                                       int          &firstNonDirGraphVertID,
                                       int          &nExtraDirichlet,
                                       BottomUpSubStructuredGraphSharedPtr &bottomUpGraph,
                                       const bool checkIfSystemSingular = false,
                                       int mdswitch = 1,
                                       bool doInteriorMap = false);

        protected:
            /// Integer map of local coeffs to global space
            Array<OneD,int> m_localToGlobalMap;
            /// Integer sign of local coeffs to global space
            Array<OneD,NekDouble> m_localToGlobalSign;
            /// Bandwith of the full matrix system (no static condensation).
            int m_fullSystemBandWidth;
            /// Integer map of process coeffs to universal space
            Array<OneD,int> m_globalToUniversalMap;
            /// Integer map of unique process coeffs to universal space (signed)
            Array<OneD,int> m_globalToUniversalMapUnique;

        private:
            /// Construct mappings for a one-dimensional scalar expansion.
            void SetUp1DExpansionC0ContMap(const int numLocalCoeffs,
                                           const ExpList &locExp,
                                           const GlobalSysSolnType solnType,
                                           const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp =
                                           NullExpListSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                                           SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap);

            /// Construct mappings for a two-dimensional scalar expansion.
            void SetUp2DExpansionC0ContMap(const int numLocalCoeffs,
                                           const ExpList &locExp,
                                           const GlobalSysSolnType solnType,
                                           const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp =
                                               NullExpListSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const vector<map<int,int> >& periodicVerticesId = NullVecIntIntMap,
                                           const map<int,int>& periodicEdgesId = NullIntIntMap,
                                           const bool checkIfSystemSingular = false);

            /// Construct mappings for a three-dimensional scalar expansion.
            void SetUp3DExpansionC0ContMap(const int numLocalCoeffs,
                                           const ExpList &locExp,
                                           const GlobalSysSolnType solnType,
                                           const Array<OneD, const ExpListSharedPtr> &bndCondExp =
                                               NullExpListSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap,
                                           const map<int,int>& periodicEdgesId = NullIntIntMap,
                                           const map<int,int>& periodicFacesId = NullIntIntMap);

            void SetUpUniversalC0ContMap(const ExpList &locExp);

            /// Calculate the bandwith of the full matrix system.
            void CalculateFullSystemBandWidth();

            virtual int v_GetLocalToGlobalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

            virtual const Array<OneD, NekDouble>& v_GetLocalToGlobalSign() const;

            virtual const void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const;

            virtual const void v_LocalToGlobal(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual const void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            virtual const void v_GlobalToLocal(
                    const NekVector<const NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            virtual const void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            virtual const void v_Assemble(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual const void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            virtual const void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            virtual const int v_GetFullSystemBandWidth() const;

        };
        typedef boost::shared_ptr<LocalToGlobalC0ContMap>  LocalToGlobalC0ContMapSharedPtr;


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_LOCALTOGLOBALC0CONTMAP_H

/**
* $Log: LocalToGlobalC0ContMap.h,v $
* Revision 1.9  2009/10/30 14:02:55  pvos
* Multi-level static condensation updates
*
* Revision 1.8  2009/05/10 23:17:12  sherwin
* Updated mainly to handle doubly periodic meshes which required modification to vertex handling from a numbering perspective
*
* Revision 1.7  2009/04/02 13:06:42  sherwin
* Modified to take symmetric banded system for HDH solver
*
* Revision 1.6  2009/01/06 21:04:42  sherwin
* Constified GlobalToLocal, LocalToGlobal and Assemble calls that take input and output argument
*
* Revision 1.5  2008/12/17 17:08:53  pvos
* Performance updates
*
* Revision 1.4  2008/11/05 16:15:24  pvos
* Added bandwith calculation routines
*
* Revision 1.3  2008/10/04 19:53:04  sherwin
* Added check that input and output arrays are different
*
* Revision 1.2  2008/09/17 13:46:40  pvos
* Added LocalToGlobalC0ContMap for 3D expansions
*
* Revision 1.1  2008/09/16 13:36:06  pvos
* Restructured the LocalToGlobalMap classes
*
**/

