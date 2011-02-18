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

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalBaseMap.h>
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
            LocalToGlobalC0ContMap();

            /// General constructor for expansions of all dimensions without
            /// boundary conditions.
            LocalToGlobalC0ContMap(const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType);

            /// Constructor for the 1D expansion mappings with boundary
            /// conditions.
            LocalToGlobalC0ContMap(const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType,
                                   const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId);

            /// Constructor for the 2D expansion mappings with boundary
            /// conditions.
            LocalToGlobalC0ContMap(const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType,
                                   const Array<OneD, const ExpList1DSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const vector<map<int,int> >& periodicVerticesId,
                                   const map<int,int>& periodicEdgesId);



            /// Constructor for the 3D expansion mappings with boundary
            /// conditions.
            LocalToGlobalC0ContMap(const int numLocalCoeffs,
                                   const ExpList &locExp,
                                   const GlobalSysSolnType solnType,
                                   const Array<OneD, const ExpList2DSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId,
                                   const map<int,int>& periodicEdgesId,
                                   const map<int,int>& periodicFacesId);

            /// Destructor.
            ~LocalToGlobalC0ContMap();


            /** Construct optimal ordering a two-dimensional expansion
            /*  given a vector of boundary condition information
            */
            void SetUp2DGraphC0ContMap(
                                       const ExpList  &locExp,
                                       const GlobalSysSolnType solnType,
                                       const Array<OneD, const MultiRegions::ExpList1DSharedPtr>  &bndCondExp,
                                       const Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >
                                       &bndConditions,
                                       const vector<map<int,int> >& periodicVerticesId,
                                       const map<int,int>& periodicEdgesId,
                                       map<int,int> &vertReorderedGraphVertId,
                                       map<int,int> &edgeReorderedGraphVertId,
                                       int          &firstNonDirGraphVertID,
                                       BottomUpSubStructuredGraphSharedPtr &bottomUpGraph,
                                       map<int,int> &interiorReorderedGraphVertId = NullIntIntMap);


            inline int GetLocalToGlobalMap(const int i) const;

            inline const Array<OneD,const int>&  GetLocalToGlobalMap();

            inline NekDouble GetLocalToGlobalSign(const int i) const;

            inline const void LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const;

            inline const void LocalToGlobal(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            inline const void GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            inline const void GlobalToLocal(
                    const NekVector<const NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            inline const void Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            inline const void Assemble(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            inline int GetFullSystemBandWidth() const;

        protected:
            /// Integer map of local coeffs to global space
            Array<OneD,int> m_localToGlobalMap;
            /// Integer sign of local coeffs to global space
            Array<OneD,NekDouble> m_localToGlobalSign;
            /// Bandwith of the full matrix system (no static condensation).
            int m_fullSystemBandWidth;

        private:
            /// Construct mappings for a one-dimensional scalar expansion.
            void SetUp1DExpansionC0ContMap(const int numLocalCoeffs,
                                           const ExpList &locExp,
                                           const GlobalSysSolnType solnType,
                                           const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp =
                                               LocalRegions::NullPointExpSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap);

            /// Construct mappings for a two-dimensional scalar expansion.
            void SetUp2DExpansionC0ContMap(const int numLocalCoeffs,
                                           const ExpList &locExp,
                                           const GlobalSysSolnType solnType,
                                           const Array<OneD, const MultiRegions::ExpList1DSharedPtr> &bndCondExp =
                                               NullExpList1DSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const vector<map<int,int> >& periodicVerticesId = NullVecIntIntMap,
                                           const map<int,int>& periodicEdgesId = NullIntIntMap);

            /// Construct mappings for a three-dimensional scalar expansion.
            void SetUp3DExpansionC0ContMap(const int numLocalCoeffs,
                                           const ExpList &locExp,
                                           const GlobalSysSolnType solnType,
                                           const Array<OneD, const ExpList2DSharedPtr> &bndCondExp =
                                               NullExpList2DSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions =
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap,
                                           const map<int,int>& periodicEdgesId = NullIntIntMap,
                                           const map<int,int>& periodicFacesId = NullIntIntMap);

            /// Calculate the bandwith of the full matrix system.
            void CalculateFullSystemBandWidth();
        };
        typedef boost::shared_ptr<LocalToGlobalC0ContMap>  LocalToGlobalC0ContMapSharedPtr;

        int LocalToGlobalC0ContMap::GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalMap[i];
        }

        inline const Array<OneD,const int>&
                    LocalToGlobalC0ContMap::GetLocalToGlobalMap(void)
        {
            return m_localToGlobalMap;
        }

        inline NekDouble LocalToGlobalC0ContMap::GetLocalToGlobalSign(
                    const int i) const
        {
            if(m_signChange)
            {
                return m_localToGlobalSign[i];
            }
            else
            {
                return 1.0;
            }
        }

        inline const void LocalToGlobalC0ContMap::LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const
        {
            if(m_signChange)
            {
                Vmath::Scatr(m_numLocalCoeffs, m_localToGlobalSign.get(), loc.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Scatr(m_numLocalCoeffs, loc.get(), m_localToGlobalMap.get(), global.get());
            }
        }

        inline const void LocalToGlobalC0ContMap::LocalToGlobal(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            LocalToGlobal(loc.GetPtr(),global.GetPtr());
        }

        inline const void LocalToGlobalC0ContMap::GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            if(m_signChange)
            {
                Vmath::Gathr(m_numLocalCoeffs, m_localToGlobalSign.get(), global.get(), m_localToGlobalMap.get(), loc.get());
            }
            else
            {
                Vmath::Gathr(m_numLocalCoeffs, global.get(), m_localToGlobalMap.get(), loc.get());
            }
        }

        inline const void LocalToGlobalC0ContMap::GlobalToLocal(
                    const NekVector<const NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocal(global.GetPtr(),loc.GetPtr());
        }

        inline const void LocalToGlobalC0ContMap::Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            ASSERTL1(loc.get() != global.get(),"Local and Global Arrays cannot be the same");

            Vmath::Zero(m_numGlobalCoeffs, global.get(), 1);

            if(m_signChange)
            {
                Vmath::Assmb(m_numLocalCoeffs, m_localToGlobalSign.get(), loc.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Assmb(m_numLocalCoeffs, loc.get(), m_localToGlobalMap.get(), global.get());
            }
        }

        inline const void LocalToGlobalC0ContMap::Assemble(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            Assemble(loc.GetPtr(),global.GetPtr());
        }

        inline int LocalToGlobalC0ContMap::GetFullSystemBandWidth() const
        {
            return m_fullSystemBandWidth;
        }

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

