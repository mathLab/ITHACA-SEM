///////////////////////////////////////////////////////////////////////////////
//
// File LinearElasticSystem.cpp
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
// Description: LinearElasticSystem solve routines 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <LinearElasticSolver/EquationSystems/LinearElasticSystem.h>

namespace Nektar
{
    string LinearElasticSystem::className = GetEquationSystemFactory().
        RegisterCreatorFunction("LinearElasticSystem",
                                LinearElasticSystem::create);

    LinearElasticSystem::LinearElasticSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession)
    {
    }

    void LinearElasticSystem::v_InitObject()
    {
        EquationSystem::v_InitObject();

        const int nVel = m_fields[0]->GetCoordim(0);
        int n;

        ASSERTL0(nVel == 2, "Linear elastic solver not set up for"
                            " this dimension (only 2D supported).");

        MultiRegions::ContField2DSharedPtr u = boost::dynamic_pointer_cast<
            MultiRegions::ContField2D>(m_fields[0]);
        m_assemblyMap = MemoryManager<CoupledAssemblyMap>
            ::AllocateSharedPtr(m_session,
                                m_graph,
                                u->GetLocalToGlobalMap(),
                                m_boundaryConditions,
                                m_fields);

        // Figure out size of matrices by looping over expansion.
        const int nEl = m_fields[0]->GetExpSize();
        LocalRegions::ExpansionSharedPtr exp;

        Array<OneD, unsigned int> sizeBnd(nEl);
        Array<OneD, unsigned int> sizeInt(nEl);

        for (n = 0; n < nEl; ++n)
        {
            exp = m_fields[0]->GetExp(m_fields[0]->GetOffset_Elmt_Id(n));
            sizeBnd[n] = nVel * exp->NumBndryCoeffs();
            sizeInt[n] = nVel * exp->GetNcoeffs() - sizeBnd[n];
        }

        // Create block matrices.
        MatrixStorage blkmatStorage = eDIAGONAL;
        m_schurCompl = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeBnd, sizeBnd, blkmatStorage);
        m_BinvD      = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeBnd, sizeInt, blkmatStorage);
        m_C          = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeInt, sizeBnd, blkmatStorage);
        m_Dinv       = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeInt, sizeInt, blkmatStorage);

        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 1.0;

        // Loop over each element and construct matrices.
        for (n = 0; n < nEl; ++n)
        {
            exp = m_fields[0]->GetExp(m_fields[0]->GetOffset_Elmt_Id(n));
            int nCoeffs = exp->GetNcoeffs();

            LocalRegions::MatrixKey matkey(StdRegions::eHelmholtz,
                                           exp->DetShapeType(),
                                           *exp, factors);

            //LocalRegions::MatrixKey matkey2(StdRegions::eLaplacian01e,
            //                                exp->DetShapeType(), *exp);

            // As a test, set up a Helmholtz matrix for each element.
            DNekMatSharedPtr loc_mat = exp->GenMatrix(matkey);
            DNekMatSharedPtr zero = MemoryManager<DNekMat>::AllocateSharedPtr(
                nCoeffs, nCoeffs, 0.0, eFULL);
            //DNekMatSharedPtr loc_mat2 = exp->GenMatrix(matkey2);

            Array<TwoD, DNekMatSharedPtr> mat(2,2);
            mat[0][0] = loc_mat;
            mat[0][1] = zero;
            mat[1][0] = zero;
            mat[1][1] = loc_mat;

            SetStaticCondBlock(n, exp, mat);
        }
    }

    void LinearElasticSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        EquationSystem::SessionSummary(s);
    }

    void LinearElasticSystem::v_DoSolve()
    {
        int i, j, n, nv;
        const int nVel = m_fields[0]->GetCoordim(0);

        // Now we've got the matrix system set up, create a GlobalLinSys object.
        MultiRegions::GlobalLinSysKey key(
            StdRegions::eLinearAdvectionReaction, m_assemblyMap);
        MultiRegions::GlobalLinSysSharedPtr linSys = MemoryManager<
            MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(
                key, m_fields[0], m_schurCompl, m_BinvD, m_C, m_Dinv, m_assemblyMap);

        const int nCoeffs = m_fields[0]->GetNcoeffs();
        const int nGlobDofs = boost::dynamic_pointer_cast<
            MultiRegions::ContField2D>(m_fields[0])->GetLocalToGlobalMap()
                                                   ->GetNumGlobalCoeffs();

        Array<OneD, Array<OneD, NekDouble> > forcing(nVel);
        EvaluateFunction(forcing, "Forcing");

        Array<OneD, NekDouble> forCoeffs(nVel * nCoeffs, 0.0);
        Array<OneD, NekDouble> inout    (nVel * nGlobDofs, 0.0);
        Array<OneD, NekDouble> rhs      (nVel * nGlobDofs, 0.0);

        int bndcnt = 0;

        for (nv = 0; nv < nVel; ++nv)
        {
            // Inner product of forcing
            Array<OneD, NekDouble> tmp(nCoeffs);
            m_fields[nv]->IProductWRTBase_IterPerExp(forcing[nv], tmp);

            // Scatter forcing into RHS vector
            for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
            {
                Array<OneD, unsigned int> bmap;
                Array<OneD, unsigned int> imap;
                m_fields[nv]->GetExp(i)->GetBoundaryMap(bmap);
                m_fields[nv]->GetExp(i)->GetInteriorMap(imap);
                int nBnd = bmap.num_elements();
                int nInt = imap.num_elements();

                int offset     = m_fields[nv]->GetCoeff_Offset(i);
                int nLocCoeffs = m_fields[nv]->GetExp(i)->GetNcoeffs();

                for (j = 0; j < nBnd; ++j)
                {
                    forCoeffs[nVel*offset + nv*nBnd + j] = tmp[offset+bmap[j]];
                }
                for (j = 0; j < nInt; ++j)
                {
                    forCoeffs[nVel*(offset + nBnd) + nv*nInt + j] = tmp[offset+imap[j]];
                }
            }

            // Impose Dirichlet boundary conditions
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp =
                m_fields[nv]->GetBndCondExpansions();
            const Array<OneD, const int> &bndMap =
                m_assemblyMap->GetBndCondCoeffsToGlobalCoeffsMap();

            for (i = 0; i < bndCondExp.num_elements(); ++i)
            {
                const Array<OneD,const NekDouble> &bndCoeffs = 
                    bndCondExp[i]->GetCoeffs();
                
                for (j = 0; j < bndCondExp[i]->GetNcoeffs(); ++j)
                {
                    NekDouble sign =
                        m_assemblyMap->GetBndCondCoeffsToGlobalCoeffsSign(
                            bndcnt);
                    inout[bndMap[bndcnt++]] = sign * bndCoeffs[j];
                }
            }
        }

        m_assemblyMap->Assemble(forCoeffs, rhs);

        // Negate RHS to be consistent with matrix definition
        Vmath::Neg(rhs.num_elements(), rhs, 1);

        // Solve
        linSys->Solve(rhs, inout, m_assemblyMap);

        Array<OneD, NekDouble> tmp(nVel * nCoeffs);

        // Backward transform
        m_assemblyMap->GlobalToLocal(inout, tmp);

        for (nv = 0; nv < nVel; ++nv)
        {
            // Scatter back to field degrees of freedom
            for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
            {
                Array<OneD, unsigned int> bmap;
                Array<OneD, unsigned int> imap;
                m_fields[nv]->GetExp(i)->GetBoundaryMap(bmap);
                m_fields[nv]->GetExp(i)->GetInteriorMap(imap);
                int nBnd = bmap.num_elements();
                int nInt = imap.num_elements();

                int offset     = m_fields[nv]->GetCoeff_Offset(i);
                int nLocCoeffs = m_fields[nv]->GetExp(i)->GetNcoeffs();

                for (j = 0; j < nBnd; ++j)
                {
                    m_fields[nv]->UpdateCoeffs()[offset+bmap[j]] = tmp[nVel*offset + nv*nBnd + j];
                }
                for (j = 0; j < nInt; ++j)
                {
                    m_fields[nv]->UpdateCoeffs()[offset+imap[j]] = tmp[nVel*(offset + nBnd) + nv*nInt + j];
                }
            }
            m_fields[nv]->BwdTrans(m_fields[nv]->GetCoeffs(),
                                   m_fields[nv]->UpdatePhys());
        }
    }

    void LinearElasticSystem::SetStaticCondBlock(
        const int                              n,
        const LocalRegions::ExpansionSharedPtr exp,
        Array<TwoD, DNekMatSharedPtr>         &mat)
    {
        int i, j, k, l;
        const int nVel = mat.GetRows();
        const int nB   = exp->NumBndryCoeffs();
        const int nI   = exp->GetNcoeffs() - nB;
        const int nBnd = exp->NumBndryCoeffs() * nVel;
        const int nInt = exp->GetNcoeffs() * nVel - nBnd;
        const MatrixStorage s = eFULL;

        // Get boundary and interior maps.
        Array<OneD, unsigned int> bmap, imap;
        exp->GetBoundaryMap(bmap);
        exp->GetInteriorMap(imap);

        DNekMatSharedPtr A =
            MemoryManager<DNekMat>::AllocateSharedPtr(nBnd, nBnd, 0.0, s);
        DNekMatSharedPtr B =
            MemoryManager<DNekMat>::AllocateSharedPtr(nBnd, nInt, 0.0, s);
        DNekMatSharedPtr C =
            MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nBnd, 0.0, s);
        DNekMatSharedPtr D =
            MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, 0.0, s);

        for (i = 0; i < nVel; ++i)
        {
            for (j = 0; j < nVel; ++j)
            {
                // Boundary-boundary and boundary-interior
                for (k = 0; k < nB; ++k)
                {
                    for (l = 0; l < nB; ++l)
                    {
                        (*A)(k + i*nB, l + j*nB) = (*mat[i][j])(bmap[k], bmap[l]);
                    }

                    for (l = 0; l < nI; ++l)
                    {
                        (*B)(k + i*nB, l + j*nI) = (*mat[i][j])(bmap[k], imap[l]);
                    }
                }

                // Interior-boundary / interior-interior
                for (k = 0; k < nI; ++k)
                {
                    for (l = 0; l < nB; ++l)
                    {
                        (*C)(k + i*nI, l + j*nB) = (*mat[i][j])(imap[k], bmap[l]);
                    }

                    for (l = 0; l < nI; ++l)
                    {
                        (*D)(k + i*nI, l + j*nI) = (*mat[i][j])(imap[k], imap[l]);
                    }
                }
            }
        }

        D->Invert();
        (*B) = (*B)*(*D);
        (*A) = (*A) - (*B)*(*C);

        DNekScalMatSharedPtr tmp_mat;
        m_schurCompl->SetBlock(
            n, n, tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                1.0, A));
        m_BinvD     ->SetBlock(
            n, n, tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                1.0, B));
        m_C         ->SetBlock(
            n, n, tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                1.0, C));
        m_Dinv      ->SetBlock(
            n, n, tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                1.0, D));
    }
}
