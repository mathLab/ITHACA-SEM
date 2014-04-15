///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLcoalToGlobalC0ContMap.cpp
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
// Description: Wrapper class around the library
// LocalToGlobalC0ContMap class for use in the Couplied Linearised NS
// solver.
///////////////////////////////////////////////////////////////////////////////

#include <LinearElasticSolver/EquationSystems/CoupledAssemblyMap.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{    
    /** 
     * Take an existing assembly map and create a coupled version suitable for
     * use in the linear elasticity solver.
     */
    CoupledAssemblyMap::CoupledAssemblyMap(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const SpatialDomains::MeshGraphSharedPtr          &graph,
        const MultiRegions::AssemblyMapCGSharedPtr        &cgMap,
        const SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields) :
        AssemblyMapCG(pSession)
    {
        int nVel = fields[0]->GetCoordim(0);

        // Multi-level static condensation doesn't work yet.
        ASSERTL0(m_solnType != MultiRegions::eDirectMultiLevelStaticCond    &&
                 m_solnType != MultiRegions::eIterativeMultiLevelStaticCond &&
                 m_solnType != MultiRegions::eXxtMultiLevelStaticCond,
                 "Multi-level static condensation not supported.");

        // Copy various number of coefficient counts.
        m_numLocalDirBndCoeffs      = cgMap->GetNumLocalDirBndCoeffs()  * nVel;
        m_numLocalBndCoeffs         = cgMap->GetNumLocalBndCoeffs()     * nVel;
        m_numLocalCoeffs            = cgMap->GetNumLocalCoeffs()        * nVel;
        m_numGlobalDirBndCoeffs     = cgMap->GetNumGlobalDirBndCoeffs() * nVel;
        m_signChange                = cgMap->GetSignChange();
        m_systemSingular            = cgMap->GetSingularSystem();

        // Copy static condensation information
        m_staticCondLevel           = cgMap->GetStaticCondLevel();
        m_numPatches                = cgMap->GetNumPatches();
        m_numLocalBndCoeffsPerPatch = cgMap->GetNumLocalBndCoeffsPerPatch();
        m_numLocalIntCoeffsPerPatch = cgMap->GetNumLocalIntCoeffsPerPatch();

        // Set up local to global and boundary condition maps.
        const int nLocBndCondDofs = cgMap->
            GetBndCondCoeffsToGlobalCoeffsMap().num_elements() * nVel;

        ASSERTL0(nLocBndCondDofs == m_numLocalDirBndCoeffs,
                 "Only Dirichlet boundary conditions are supported");
        
        m_localToGlobalMap               =
            Array<OneD, int>(m_numLocalCoeffs,-1);
        m_localToGlobalBndMap            =
            Array<OneD, int>(m_numLocalBndCoeffs,-1);
        m_bndCondCoeffsToGlobalCoeffsMap =
            Array<OneD, int>(nLocBndCondDofs,-1);

        if(m_signChange)
        {
            m_localToGlobalSign               =
                Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
            m_localToGlobalBndSign            =
                Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
            m_bndCondCoeffsToGlobalCoeffsSign =
                Array<OneD, NekDouble>(nLocBndCondDofs,1.0);
        }
        else
        {
            m_localToGlobalSign               = NullNekDouble1DArray;
            m_localToGlobalBndSign            = NullNekDouble1DArray;
            m_bndCondCoeffsToGlobalCoeffsSign = NullNekDouble1DArray;
        }

        const int nGlobBndCoeffs = m_numGlobalBndCoeffs    / nVel;
        const int nGlobDirCoeffs = m_numGlobalDirBndCoeffs / nVel;

        // Set up local to global boundary mapping.
        const LocalRegions::ExpansionVector &locExpVector = *(fields[0]->GetExp());
        int i, j, n, cnt1, cnt2;

        cnt1 = 0;
        for (n = 0; n < nVel; ++n)
        {
            cnt2 = 0;
            for (i = 0; i < locExpVector.size(); ++i)
            {
                const int nBndCoeffs = locExpVector[i]->NumBndryCoeffs();
                for (j = 0; j < nBndCoeffs; ++j, ++cnt1, ++cnt2)
                {
                    const int l2g = cgMap->GetLocalToGlobalBndMap()[cnt2];

                    if (l2g < nGlobDirCoeffs)
                    {
                        m_localToGlobalMap[cnt1] = n * nGlobDirCoeffs + l2g;
                    }
                    else
                    {
                        m_localToGlobalMap[cnt1] =
                            nVel * nGlobDirCoeffs + n * nGlobBndCoeffs + l2g;
                    }

                    if (m_signChange)
                    {
                        m_localToGlobalSign[cnt1] =
                            cgMap->GetLocalToGlobalBndSign()[cnt2];
                    }
                }
            }
        }

        // Set up local to global mapping
        const int nLocalCoeffs    = m_numLocalCoeffs    / nVel;
        const int nLocalBndCoeffs = m_numLocalBndCoeffs / nVel;
        int globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;

        for (n = 0; n < nVel; ++n)
        {
            const int off1 = n * nLocalCoeffs;
            const int off2 = n * nLocalBndCoeffs;
            cnt1 = 0;

            for (i = 0; i < nLocalCoeffs; ++i)
            {
                if (m_localToGlobalMap[off1+i] == -1)
                {
                    m_localToGlobalMap[off1+i] = globalId++;
                }
                else
                {
                    if (m_signChange)
                    {
                        m_localToGlobalBndSign[off2+cnt1] =
                            m_localToGlobalSign[off1+i];
                    }
                    m_localToGlobalBndMap[cnt1++] = m_localToGlobalMap[off1+i];
                }
            }
        }

        // Set up boundary condition mapping: this is straightforward since we
        // only consider Dirichlet boundary conditions.
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp
            = fields[0]->GetBndCondExpansions();

        const int nLocalDirBndCoeffs = m_numLocalDirBndCoeffs / nVel;

        for (n = 0; n < nVel; ++n)
        {
            const int offset = n * nLocalDirBndCoeffs;

            for (i = 0; i < nLocalDirBndCoeffs; ++i)
            {
                m_bndCondCoeffsToGlobalCoeffsMap[offset+i] =
                    cgMap->GetBndCondCoeffsToGlobalCoeffsMap()[i];

                if (m_signChange)
                {
                    m_bndCondCoeffsToGlobalCoeffsSign[offset+i] =
                        cgMap->GetBndCondCoeffsToGlobalCoeffsSign(i);
                }
            }
        }
        
        m_hash = boost::hash_range(
            m_localToGlobalMap.begin(), m_localToGlobalMap.end());
    }
}
