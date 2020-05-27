///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapDG.h
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
// Description: Local to Global DG mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLY_MAP_DG_H
#define MULTIREGIONS_ASSEMBLY_MAP_DG_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <MultiRegions/AssemblyMap/AssemblyCommDG.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList0D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class AssemblyMapDG;
        typedef std::shared_ptr<AssemblyMapDG>  AssemblyMapDGSharedPtr;

        ///
        class AssemblyMapDG: public AssemblyMap
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMapDG();

            /// Constructor for trace map for one-dimensional expansion.
            MULTI_REGIONS_EXPORT AssemblyMapDG( 
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const ExpListSharedPtr &trace,
                const ExpList &locExp,
                const Array<OneD, const MultiRegions::ExpListSharedPtr>
                                                                &bndConstraint,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndCond,
                const PeriodicMap &periodicTrace,
                const std::string variable = "DefaultVar");

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMapDG();

            /// Return the number of boundary segments on which Dirichlet
            /// boundary conditions are imposed.
            MULTI_REGIONS_EXPORT int GetNumDirichletBndPhys();

            MULTI_REGIONS_EXPORT Array<OneD, LocalRegions::ExpansionSharedPtr>
                &GetElmtToTrace(const int i);

            MULTI_REGIONS_EXPORT 
                Array<OneD,Array<OneD,LocalRegions::ExpansionSharedPtr> >
                &GetElmtToTrace();

            /**
             * Changes toAlign quadrature point order, where the realignment is
             * given by orient, which defines the mapping needed to go between
             * the original ordering and the new desired ordering.
             *
             * @param[in,out] toAlign Data to reorder
             * @param[in] orient The transformation to perform
             * @param[in] nquad1 Quadrature points in direction 1
             * @param[in] nquad2 Quadrature points in direction 2
             */
            MULTI_REGIONS_EXPORT static void RealignTraceElement(
                    Array<OneD, int>        &toAlign,
                    StdRegions::Orientation  orient,
                    int                      nquad1,
                    int                      nquad2 = 0);

            MULTI_REGIONS_EXPORT AssemblyCommDGSharedPtr GetAssemblyCommDG();

        protected:
            /// Number of physical dirichlet boundary values in trace
            int m_numDirichletBndPhys;

            AssemblyCommDGSharedPtr m_assemblyComm;

            /// list of edge expansions for a given element
            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> > m_elmtToTrace;

            void SetUpUniversalDGMap(const ExpList &locExp);

            virtual int v_GetLocalToGlobalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

            virtual void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm = false) const;

            virtual void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            virtual void v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            virtual void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            virtual void v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            virtual void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            virtual int v_GetFullSystemBandWidth() const;
        }; // class


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLY_MAP_DG_H
