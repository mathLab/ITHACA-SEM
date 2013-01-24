///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG.h
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
// Description: C0-continuous Local to Global mapping routines, base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLYMAPCG_H
#define MULTIREGIONS_ASSEMBLYMAPCG_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>

namespace Nektar
{
    namespace MultiRegions
    {

        static map<int,int> NullIntIntMap;
        const static vector<map<int,int> > NullVecIntIntMap;

        class ExpList;
        class AssemblyMapCG;
        typedef boost::shared_ptr<AssemblyMapCG>  AssemblyMapCGSharedPtr;



        /// Constructs mappings for the C0 scalar continuous Galerkin formulation.
        class AssemblyMapCG: public AssemblyMap
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMapCG(
                                    const LibUtilities::SessionReaderSharedPtr &pSession);


            /// General constructor for expansions of all dimensions without
            /// boundary conditions.
            MULTI_REGIONS_EXPORT AssemblyMapCG(
                                    const LibUtilities::SessionReaderSharedPtr &pSession,
                                    const int numLocalCoeffs,
                                    const ExpList &locExp);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMapCG();

            MULTI_REGIONS_EXPORT map<int, vector<pair<int, int> > > 
                &GetExtraDirDofs()
            {
                return m_extraDirDofs;
            }

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
            /// Number of non Dirichlet vertex modes
            int m_numNonDirVertexModes;
            /// Number of non Dirichlet edge modes
            int m_numNonDirEdgeModes;
            /// Number of non Dirichlet face modes
            int m_numNonDirFaceModes;
            /// Maximum static condensation level.
            int m_maxStaticCondLevel;
            map<int, vector<pair<int, int> > > m_extraDirDofs;
            
            void SetUpUniversalC0ContMap(const ExpList &locExp);

            /// Calculate the bandwith of the full matrix system.
            void CalculateFullSystemBandWidth();

            MULTI_REGIONS_EXPORT virtual int v_GetLocalToGlobalMap(const int i) const;

            MULTI_REGIONS_EXPORT virtual int v_GetGlobalToUniversalMap(const int i) const;

            MULTI_REGIONS_EXPORT virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            MULTI_REGIONS_EXPORT virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            MULTI_REGIONS_EXPORT virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            MULTI_REGIONS_EXPORT virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            MULTI_REGIONS_EXPORT virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

            MULTI_REGIONS_EXPORT virtual const Array<OneD, NekDouble>& v_GetLocalToGlobalSign() const;

            MULTI_REGIONS_EXPORT virtual const void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const;

            MULTI_REGIONS_EXPORT virtual const void v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            MULTI_REGIONS_EXPORT virtual const void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT virtual const void v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT virtual const void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            MULTI_REGIONS_EXPORT virtual const void v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            MULTI_REGIONS_EXPORT virtual const void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT virtual const void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT virtual const void v_UniversalAssemble(
                Array<OneD,     NekDouble>& pGlobal,
                int offset) const;

            MULTI_REGIONS_EXPORT virtual const int v_GetFullSystemBandWidth() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirVertexModes() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirEdgeModes() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirFaceModes() const;

        };


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLYMAPCG_H

