///////////////////////////////////////////////////////////////////////////////
//
// File: CoupledAssemblyMap.cpp
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
// Description: Coupled assembly map for linear elasticity solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <LinearElasticSolver/EquationSystems/CoupledAssemblyMap.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

#include <iomanip>

namespace Nektar
{

/**
 * @brief Take an existing assembly map and create a coupled version suitable
 * for use in the linear elasticity solver.
 *
 * The linear elasticity solver requires a slight reordering of local and global
 * coefficients to support problems of the form
 *
 * [ A B C ] [ u ] = [ f_u ]
 * [ D E F ] [ v ]   [ f_v ]
 * [ G H I ] [ w ]   [ f_w ]
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

    // Copy various coefficient counts, and multiply by the dimension of the
    // problem to obtain our new values.
    m_numLocalDirBndCoeffs      = cgMap->GetNumLocalDirBndCoeffs()  * nVel;
    m_numLocalBndCoeffs         = cgMap->GetNumLocalBndCoeffs()     * nVel;
    m_numLocalCoeffs            = cgMap->GetNumLocalCoeffs()        * nVel;
    m_numGlobalBndCoeffs        = cgMap->GetNumGlobalBndCoeffs()    * nVel;
    m_numGlobalDirBndCoeffs     = cgMap->GetNumGlobalDirBndCoeffs() * nVel;
    m_numGlobalCoeffs           = cgMap->GetNumGlobalCoeffs()       * nVel;
    m_signChange                = cgMap->GetSignChange();
    m_systemSingular            = cgMap->GetSingularSystem();

    // Copy static condensation information. TODO: boundary and interior patches
    // need to be re-ordered in order to allow for multi-level static
    // condensation support.
    m_staticCondLevel           = cgMap->GetStaticCondLevel();
    m_lowestStaticCondLevel     = cgMap->GetLowestStaticCondLevel();
    m_numPatches                = cgMap->GetNumPatches();
    m_numLocalBndCoeffsPerPatch = cgMap->GetNumLocalBndCoeffsPerPatch();
    m_numLocalIntCoeffsPerPatch = cgMap->GetNumLocalIntCoeffsPerPatch();

    // Set up local to global and boundary condition maps.
    const int nLocBndCondDofs = cgMap->
        GetBndCondCoeffsToGlobalCoeffsMap().num_elements() * nVel;

    // Allocate storage for local to global maps.
    m_localToGlobalMap               =
        Array<OneD, int>(m_numLocalCoeffs,-1);
    m_localToGlobalBndMap            =
        Array<OneD, int>(m_numLocalBndCoeffs,-1);
    m_bndCondCoeffsToGlobalCoeffsMap =
        Array<OneD, int>(nLocBndCondDofs,-1);

    // Only require a sign map if we are using modal polynomials in the
    // expansion and the order is >= 3.
    if (m_signChange)
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

    const LocalRegions::ExpansionVector &locExpVector
        = *(fields[0]->GetExp());

    map<int, int> newGlobalIds;
    int i, j, n, cnt1, cnt2;

    // Order local boundary degrees of freedom. These are basically fine; we
    // reorder storage so that we loop over each element and then each component
    // of velocity, by applying a mapping l2g -> nVel*l2g + n, for 0 <= n <
    // nVel. Note that Dirichlet ordering is preserved under this
    // transformation.
    cnt1 = cnt2 = 0;
    for (i = 0; i < locExpVector.size(); ++i)
    {
        const int nBndCoeffs = locExpVector[i]->NumBndryCoeffs();

        for (n = 0; n < nVel; ++n)
        {
            for (j = 0; j < nBndCoeffs; ++j, ++cnt1)
            {
                const int l2g = cgMap->GetLocalToGlobalBndMap()[cnt2+j];
                m_localToGlobalBndMap[cnt1] = nVel * l2g + n;

                if (m_signChange)
                {
                    m_localToGlobalBndSign[cnt1] =
                        cgMap->GetLocalToGlobalBndSign()[cnt2+j];
                }

                if (n == 0)
                {
                    const int l2gnew = m_localToGlobalBndMap[cnt1];
                    if (newGlobalIds.count(l2g))
                    {
                        ASSERTL1(newGlobalIds[l2g] == l2gnew,
                                 "Consistency error");
                    }
                    newGlobalIds[l2g] = l2gnew;
                }
            }
        }

        cnt2 += nBndCoeffs;
    }

    // Grab map of extra Dirichlet degrees of freedom for parallel runs
    // (particularly in 3D).
    m_extraDirDofs = cgMap->GetExtraDirDofs();

    // Counter for remaining interior degrees of freedom.
    int globalId = m_numGlobalBndCoeffs;

    // Interior degrees of freedom are a bit more tricky -- global linear system
    // solve relies on them being in the same order as the BinvD, C and invD
    // matrices.
    cnt1 = cnt2 = 0;
    for (i = 0; i < locExpVector.size(); ++i)
    {
        const int nCoeffs    = locExpVector[i]->GetNcoeffs();
        const int nBndCoeffs = locExpVector[i]->NumBndryCoeffs();

        for (n = 0; n < nVel; ++n)
        {
            for (j = 0; j < nBndCoeffs; ++j, ++cnt1, ++cnt2)
            {
                const int l2g = m_localToGlobalBndMap[cnt2];
                m_localToGlobalMap[cnt1] = l2g;

                if (m_signChange)
                {
                    m_localToGlobalSign[cnt1] = m_localToGlobalBndSign[cnt2];
                }
            }
        }

        for (n = 0; n < nVel; ++n)
        {
            for (j = 0; j < nCoeffs - nBndCoeffs; ++j, ++cnt1)
            {
                m_localToGlobalMap[cnt1] = globalId++;
            }
        }
    }

    for (i = 0; i < m_localToGlobalMap.num_elements(); ++i)
    {
        ASSERTL1(m_localToGlobalMap[i] != -1, "Consistency error");
    }

    ASSERTL1(globalId == m_numGlobalCoeffs, "Consistency error");

    const int nLocalDirBndCoeffs =
        cgMap->GetBndCondCoeffsToGlobalCoeffsMap().num_elements();

    cnt1 = 0;
    for (n = 0; n < nVel; ++n)
    {
        for (i = 0; i < nLocalDirBndCoeffs; ++i, ++cnt1)
        {
            const int l2g = cgMap->GetBndCondCoeffsToGlobalCoeffsMap()[i];
            int newId = newGlobalIds[l2g];
            m_bndCondCoeffsToGlobalCoeffsMap[cnt1] = newId + n;

            if (m_signChange)
            {
                m_bndCondCoeffsToGlobalCoeffsSign[cnt1] =
                    cgMap->GetBndCondCoeffsToGlobalCoeffsSign(i);
            }
        }
    }

    // Finally, set up global to universal maps.
    m_globalToUniversalMap          = Array<OneD, int>(m_numGlobalCoeffs);
    m_globalToUniversalMapUnique    = Array<OneD, int>(m_numGlobalCoeffs);
    m_globalToUniversalBndMap       = Array<OneD, int>(m_numGlobalBndCoeffs);
    m_globalToUniversalBndMapUnique = Array<OneD, int>(m_numGlobalBndCoeffs);

    for (i = 0; i < cgMap->GetNumGlobalBndCoeffs(); ++i)
    {
        for (n = 0; n < nVel; ++n)
        {
            m_globalToUniversalBndMap[i*nVel + n] =
                cgMap->GetGlobalToUniversalBndMap()[i]*nVel + n;
            m_globalToUniversalMap[i*nVel + n] =
                cgMap->GetGlobalToUniversalBndMap()[i]*nVel + n;
        }
    }

    Array<OneD, long> tmp(m_numGlobalCoeffs);
    Vmath::Zero(m_numGlobalCoeffs, tmp, 1);
    Array<OneD, long> tmp2(m_numGlobalBndCoeffs, tmp);
    for (unsigned int i = 0; i < m_numGlobalBndCoeffs; ++i)
    {
        tmp[i] = m_globalToUniversalBndMap[i];
    }

    LibUtilities::CommSharedPtr vCommRow = m_comm->GetRowComm();
    m_gsh    = Gs::Init(tmp,  vCommRow);
    m_bndGsh = Gs::Init(tmp2, vCommRow);
    Gs::Unique(tmp, vCommRow);
    for (unsigned int i = 0; i < m_numGlobalCoeffs; ++i)
    {
        m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
    }
    for (unsigned int i = 0; i < m_numGlobalBndCoeffs; ++i)
    {
        m_globalToUniversalBndMapUnique[i] = (tmp2[i] >= 0 ? 1 : 0);
    }

    m_hash = boost::hash_range(
        m_localToGlobalMap.begin(), m_localToGlobalMap.end());
}

}
