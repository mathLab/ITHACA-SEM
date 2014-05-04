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

        // For now only two dimensions are supported. The code below and in the
        // assembly map should be readily extendible to three dimensions
        // however.
        ASSERTL0(nVel == 2, "Linear elastic solver not set up for"
                            " this dimension (only 2D supported).");

        // Make sure that we have Young's modulus and Poisson ratio set.
        m_session->LoadParameter("E", m_E, 1.0);
        m_session->LoadParameter("nu", m_nu, 0.25);

        // Create a coupled assembly map which allows us to tie u and v fields
        // together.
        MultiRegions::ContField2DSharedPtr u = boost::dynamic_pointer_cast<
            MultiRegions::ContField2D>(m_fields[0]);
        m_assemblyMap = MemoryManager<CoupledAssemblyMap>
            ::AllocateSharedPtr(m_session,
                                m_graph,
                                u->GetLocalToGlobalMap(),
                                m_boundaryConditions,
                                m_fields);

        // Figure out size of our new matrix systems by looping over all
        // expansions and multiply number of coefficients by velocity
        // components.
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

        // Create block matrix storage for the statically condensed system.
        MatrixStorage blkmatStorage = eDIAGONAL;
        m_schurCompl = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeBnd, sizeBnd, blkmatStorage);
        m_BinvD      = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeBnd, sizeInt, blkmatStorage);
        m_C          = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeInt, sizeBnd, blkmatStorage);
        m_Dinv       = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(
            sizeInt, sizeInt, blkmatStorage);

        // Build initial matrix system.
        BuildMatrixSystem();
    }

    void LinearElasticSystem::BuildMatrixSystem()
    {
        const int nEl = m_fields[0]->GetExpSize();

        LocalRegions::ExpansionSharedPtr exp;
        int n;

        // Factors map for matrix keys.
        StdRegions::ConstFactorMap factors;

        // Calculate various constants
        NekDouble a = m_E * (1.0 - m_nu) / (1.0 + m_nu) / (1.0 - 2.0*m_nu);
        NekDouble b = m_E * 0.5 / (1.0 + m_nu);
        NekDouble c = m_E * m_nu / (1.0 + m_nu) / (1.0 - 2.0*m_nu);

        // Loop over each element and construct matrices.
        for (n = 0; n < nEl; ++n)
        {
            exp = m_fields[0]->GetExp(m_fields[0]->GetOffset_Elmt_Id(n));
            const int nPhys = exp->GetTotPoints();

            StdRegions::VarCoeffMap varcoeffA, varcoeffD;
            varcoeffA[StdRegions::eVarCoeffD00] =
                Array<OneD, NekDouble>(nPhys, a);
            varcoeffA[StdRegions::eVarCoeffD11] =
                Array<OneD, NekDouble>(nPhys, b);
            varcoeffD[StdRegions::eVarCoeffD00] =
                Array<OneD, NekDouble>(nPhys, b);
            varcoeffD[StdRegions::eVarCoeffD11] =
                Array<OneD, NekDouble>(nPhys, a);

            LocalRegions::MatrixKey matkeyA(StdRegions::eLaplacian,
                                            exp->DetShapeType(),
                                            *exp, factors, varcoeffA);
            LocalRegions::MatrixKey matkeyD(StdRegions::eLaplacian,
                                            exp->DetShapeType(),
                                            *exp, factors, varcoeffD);

            /*
             * mat holds the linear operator [ A B ] acting on [ u ].
             *                               [ C D ]           [ v ]
             *
             * In this case it is just a diagonal matrix with each component
             * being a Helmholtz matrix, so that the u and v fields are not
             * coupled at all.
             */
            Array<TwoD, DNekMatSharedPtr> mat(2,2);
            mat[0][0] = exp->GenMatrix(matkeyA);
            mat[0][1] = BuildLaplacianIJMatrix(1, 0, c, exp);
            mat[1][0] = BuildLaplacianIJMatrix(0, 1, c, exp);
            mat[1][1] = exp->GenMatrix(matkeyD);

            // Set up the statically condensed block for this element.
            SetStaticCondBlock(n, exp, mat);
        }
    }

    void LinearElasticSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        EquationSystem::SessionSummary(s);
    }

    void LinearElasticSystem::v_DoSolve()
    {
        int i, j, nv;
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

        // Evaluate the forcing function from the XML file.
        Array<OneD, Array<OneD, NekDouble> > forcing(nVel);
        EvaluateFunction(forcing, "Forcing");

        // Set up some temporary storage.
        //
        // - forCoeffs holds the forcing coefficients in a local ordering;
        //   however note that the ordering is different and dictated by the
        //   assembly map. We loop over each element, then the boundary degrees
        //   of freedom for u, boundary for v, followed by the interior for u
        //   and then interior for v.
        // - rhs is the global assembly of forCoeffs.
        // - inout holds the Dirichlet degrees of freedom in the global
        //   ordering, which have been assembled from the boundary expansion.
        Array<OneD, NekDouble> forCoeffs(nVel * nCoeffs, 0.0);
        Array<OneD, NekDouble> inout    (nVel * nGlobDofs, 0.0);
        Array<OneD, NekDouble> rhs      (nVel * nGlobDofs, 0.0);

        // Counter for the local Dirichlet boundary to global ordering.
        int bndcnt = 0;

        for (nv = 0; nv < nVel; ++nv)
        {
            // Take the inner product of the forcing function.
            Array<OneD, NekDouble> tmp(nCoeffs);
            m_fields[nv]->IProductWRTBase_IterPerExp(forcing[nv], tmp);

            // Scatter forcing into RHS vector according to the ordering
            // dictated in the comment above.
            for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
            {
                Array<OneD, unsigned int> bmap;
                Array<OneD, unsigned int> imap;
                m_fields[nv]->GetExp(i)->GetBoundaryMap(bmap);
                m_fields[nv]->GetExp(i)->GetInteriorMap(imap);
                int nBnd = bmap.num_elements();
                int nInt = imap.num_elements();

                int offset     = m_fields[nv]->GetCoeff_Offset(i);

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
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp
                = m_fields[nv]->GetBndCondExpansions();
            const Array<OneD, const int> &bndMap
                = m_assemblyMap->GetBndCondCoeffsToGlobalCoeffsMap();

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

        // Assemble forcing into the RHS.
        m_assemblyMap->Assemble(forCoeffs, rhs);

        // Negate RHS to be consistent with matrix definition.
        Vmath::Neg(rhs.num_elements(), rhs, 1);

        // Solve.
        linSys->Solve(rhs, inout, m_assemblyMap);

        // Scatter the global ordering back to the alternate local ordering.
        Array<OneD, NekDouble> tmp(nVel * nCoeffs);
        m_assemblyMap->GlobalToLocal(inout, tmp);

        // Finally, scatter back to field degrees of freedom
        for (nv = 0; nv < nVel; ++nv)
        {
            for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
            {
                Array<OneD, unsigned int> bmap;
                Array<OneD, unsigned int> imap;
                m_fields[nv]->GetExp(i)->GetBoundaryMap(bmap);
                m_fields[nv]->GetExp(i)->GetInteriorMap(imap);
                int nBnd   = bmap.num_elements();
                int nInt   = imap.num_elements();
                int offset = m_fields[nv]->GetCoeff_Offset(i);

                for (j = 0; j < nBnd; ++j)
                {
                    m_fields[nv]->UpdateCoeffs()[offset+bmap[j]] =
                        tmp[nVel*offset + nv*nBnd + j];
                }
                for (j = 0; j < nInt; ++j)
                {
                    m_fields[nv]->UpdateCoeffs()[offset+imap[j]] =
                        tmp[nVel*(offset + nBnd) + nv*nInt + j];
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

    DNekMatSharedPtr LinearElasticSystem::BuildLaplacianIJMatrix(
        const int                        k1,
        const int                        k2,
        const NekDouble                  scale,
        LocalRegions::ExpansionSharedPtr exp)
    {
        const int nCoeffs = exp->GetNcoeffs();
        const int nPhys   = exp->GetTotPoints();
        int i;

        DNekMatSharedPtr ret = MemoryManager<DNekMat>::AllocateSharedPtr(
            nCoeffs, nCoeffs, 0.0, eFULL);

        Array<OneD, NekDouble> tmp2(nPhys);
        Array<OneD, NekDouble> tmp3(nPhys);

        for (i = 0; i < nCoeffs; ++i)
        {
            Array<OneD, NekDouble> tmp1(nCoeffs, 0.0);
            tmp1[i] = 1.0;

            exp->BwdTrans            (    tmp1, tmp2);
            exp->PhysDeriv           (k1, tmp2, tmp3);
            exp->IProductWRTDerivBase(k2, tmp3, tmp1);

            Vmath::Smul(nCoeffs, scale, &tmp1[0], 1, &(ret->GetPtr())[0]+i*nCoeffs,1);
        }

        return ret;
    }
}
