///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalDGMap.h
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
// Description: Local to Global DG mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////
#ifndef MULTIREGIONS_LOCALTOGLOBALDGMAP_H
#define MULTIREGIONS_LOCALTOGLOBALDGMAP_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalBaseMap.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList0D.h>

#include <LocalRegions/PointExp.h>

namespace Nektar
{
    namespace MultiRegions
    {
        ///
        class LocalToGlobalDGMap: public LocalToGlobalBaseMap
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT LocalToGlobalDGMap();
            /// Constructor for trace map for one-dimensional expansion.
            MULTI_REGIONS_EXPORT LocalToGlobalDGMap( 
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraph1D &graph1D,
                const ExpList &locExp,
                const Array<OneD, const MultiRegions::ExpListSharedPtr>
                                                                &bndConstraint,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndCond);
            /// Constructor for trace map for two-dimensional expansion.
            MULTI_REGIONS_EXPORT LocalToGlobalDGMap(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                      SpatialDomains::MeshGraph2D &graph2D,
                const ExpList1DSharedPtr &trace,
                const ExpList &locExp,
                const Array<OneD, MultiRegions::ExpListSharedPtr>
                                                                &bndContraint,
                const Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                                                                &bndCond,
                const map<int,int> &periodicEdges);

            /// Destructor.
            MULTI_REGIONS_EXPORT ~LocalToGlobalDGMap();

            /// Return the number of boundary segments on which Dirichlet
            /// boundary conditions are imposed.
            inline int GetNumDirichletBndPhys();

            inline Array<OneD, StdRegions::StdExpansion1DSharedPtr>
                    GetElmtToTrace(const int i);

            inline Array<OneD,Array<OneD,StdRegions::StdExpansion1DSharedPtr> >
                    GetElmtToTrace();

            inline AdjacentTraceOrientation
                    GetBndExpAdjacentOrient(const int i);

        protected:

        private:
            /// Number of physical dirichlet boundary values in trace
            int m_numDirichletBndPhys;
            /// list of edge expansions for a given element
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> > m_elmtToTrace;

            Array<OneD, AdjacentTraceOrientation > m_bndExpAdjacentOrient;

            void SetUpUniversalDGMap(const ExpList &locExp);

            virtual int v_GetLocalToGlobalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

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

        typedef boost::shared_ptr<LocalToGlobalDGMap>  LocalToGlobalDGMapSharedPtr;

        inline int LocalToGlobalDGMap::GetNumDirichletBndPhys()
        {
            return m_numDirichletBndPhys;
        }

        inline Array<OneD, StdRegions::StdExpansion1DSharedPtr>
                    LocalToGlobalDGMap::GetElmtToTrace(const int i)
        {
            ASSERTL1(i >= 0 && i < m_elmtToTrace.num_elements(),
                     "i is out of range");
            return m_elmtToTrace[i];
        }

        inline Array<OneD, Array< OneD, StdRegions::StdExpansion1DSharedPtr> >
                    LocalToGlobalDGMap::GetElmtToTrace()
        {
            return m_elmtToTrace;
        }

        inline AdjacentTraceOrientation
                    LocalToGlobalDGMap::GetBndExpAdjacentOrient(const int i)
        {
            return m_bndExpAdjacentOrient[i];
        }
    } // end of namespace
} // end of namespace

#endif //LOCALTOGLOBALDGMAP_H


/** $Log: LocalToGlobalDGMap.h,v $
/** Revision 1.6  2009/11/20 10:47:54  cbiotto
/** Update for creating boundary to global trace map
/**
/** Revision 1.5  2009/11/02 19:15:43  cantwell
/** Moved ContField1D to inherit from DisContField1D.
/** Moved ContField3D to inherit from DisContField3D.
/** Incorporated GenExpList1D functionality into ExpList1D.
/** Tidied up and added documentation to various classes.
/** Moved Namespace documentation and introductions to separate files along with
/** doxygen configuration.
/** Added option to use system ZLIB library instead of libboost_zlib on UNIX.
/** Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
/** Updated Linux compiling instructions.
/** Updated regDemo to use Helmholtz2D-g when built as debug.
/**
/** Revision 1.4  2009/10/30 14:02:55  pvos
/** Multi-level static condensation updates
/**
/** Revision 1.3  2009/04/02 13:06:42  sherwin
/** Modified to take symmetric banded system for HDH solver
/**
/** Revision 1.2  2008/09/16 13:36:06  pvos
/** Restructured the LocalToGlobalMap classes
/**
/** Revision 1.1  2008/08/18 08:16:23  sherwin
/** First version of this new class container for mappings
/**
 */

