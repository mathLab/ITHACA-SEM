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

#include <algorithm>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/Preconditioner.h>
#include <LinearElasticSolver/EquationSystems/LinearElasticSystem.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdNodalTriExp.h>

namespace Nektar
{
    string LinearElasticSystem::className = GetEquationSystemFactory().
        RegisterCreatorFunction("LinearElasticSystem",
                                LinearElasticSystem::create);

    inline void idealToRef(NekDouble &x, NekDouble &y)
    {
        x = x + 0.5 * y + 0.5;
        y = sqrt(3)/2 * y - (1 - sqrt(3)/2);
    }
    
    /* Mapping that requires the coordinates of the three vertices (x,y) of the triangular element and outputs the function coefficients that defines the mapping to the reference element for each coordinate (xfunc and yfunc)
    */
    inline void MappingIdealToRef(Array<OneD, NekDouble> x,
                                  Array<OneD, NekDouble> y,
                                  Array<OneD, NekDouble> xf,
                                  Array<OneD, NekDouble> yf)
    {
        char jobvl = 'N', jobvr = 'V';
        
        int n = x.num_elements();
        int worklen = 8*n, info;
        Array<OneD, NekDouble> tmp ((n+1)*(n+1), 0.0);
        
        DNekMat map(n, n, 0.0, eFULL);
        DNekMat mapinv(n, n, 0.0, eFULL);
        
        map(0,0) = x[0];map(0,1) = y[0];map(0,2) = 1.0;
        map(1,0) = x[1];map(1,1) = y[1];map(1,2) = 1.0;
        map(2,0) = x[2];map(2,1) = y[2];map(2,2) = 1.0;
        
        mapinv = map;
        mapinv.Invert();
        
        Array<OneD, NekDouble> xref (n,0.0), yref (n,0.0);
        //Array<OneD, NekDouble> xf   (n,0.0),   yf (n,0.0);
        
        xref[0] = -1.0;xref[1] =  1.0;xref[2] = -1.0;
        yref[0] = -1.0;yref[1] = -1.0;yref[2] =  1.0;
        
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                xf[i] += mapinv(i,j)*xref[j];
                yf[i] += mapinv(i,j)*yref[j];
            }
        }
    }

    LinearElasticSystem::LinearElasticSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession)
    {
    }

    /**
     * @brief Set up the linear elasticity system.
     *
     * This routine loads the E and nu variables from the session file, creates
     * a coupled assembly map and creates containers for the statically
     * condensed block matrix system.
     */
    void LinearElasticSystem::v_InitObject()
    {
        EquationSystem::v_InitObject();
        
        const int nVel = m_fields[0]->GetCoordim(0);
        int n;

        ASSERTL0(nVel > 1, "Linear elastic solver not set up for"
                           " this dimension (only 2D/3D supported).");
       
        // Make sure that we have Young's modulus and Poisson ratio set.
        m_session->LoadParameter("E", m_E, 1.0);
        m_session->LoadParameter("nu", m_nu, 0.25);
        m_session->LoadParameter("Beta", m_beta, 1.0);
        
        // Create a coupled assembly map which allows us to tie u, v and w
        // fields together.
        if (nVel == 2)
        {
            MultiRegions::ContField2DSharedPtr u = boost::dynamic_pointer_cast<
                MultiRegions::ContField2D>(m_fields[0]);
            m_assemblyMap = MemoryManager<CoupledAssemblyMap>
                ::AllocateSharedPtr(m_session,
                                    m_graph,
                                    u->GetLocalToGlobalMap(),
                                    m_boundaryConditions,
                                    m_fields);
        }
        
        if (nVel == 3)
        {
            MultiRegions::ContField3DSharedPtr u = boost::dynamic_pointer_cast<
                MultiRegions::ContField3D>(m_fields[0]);
            m_assemblyMap = MemoryManager<CoupledAssemblyMap>
                ::AllocateSharedPtr(m_session,
                                    m_graph,
                                    u->GetLocalToGlobalMap(),
                                    m_boundaryConditions,
                                    m_fields);
        }
        
        // Figure out size of our new matrix systems by looping over all
        // expansions and multiply number of coefficients by velocity
        // components.
        const int nEl = m_fields[0]->GetExpSize();
        LocalRegions::ExpansionSharedPtr exp;

        Array<OneD, unsigned int> sizeBnd(nEl);
        Array<OneD, unsigned int> sizeInt(nEl);

        // Allocate storage for boundary and interior map caches.
        m_bmap = Array<OneD, Array<OneD, unsigned int> >(nEl);
        m_imap = Array<OneD, Array<OneD, unsigned int> >(nEl);

        // Cache interior and boundary maps to avoid recalculating
        // constantly. Should really be handled by manager in StdRegions but not
        // really working yet.
        for (n = 0; n < nEl; ++n)
        {
            exp = m_fields[0]->GetExp(m_fields[0]->GetOffset_Elmt_Id(n));
            sizeBnd[n] = nVel * exp->NumBndryCoeffs();
            sizeInt[n] = nVel * exp->GetNcoeffs() - sizeBnd[n];
            exp->GetBoundaryMap(m_bmap[n]);
            exp->GetInteriorMap(m_imap[n]);
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
    }

    /**
     * @brief Build matrix system for linear elasticity equations.
     *
     * This routine constructs the matrix discretisation arising from the weak
     * formulation of the linear elasticity equations. The resulting matrix
     * system is then passed to LinearElasticSystem::SetStaticCondBlock in order
     * to construct the statically condensed system.
     *
     * All of the matrices involved in the construction of the divergence of the
     * stress tensor are Laplacian-like. We use the variable coefficient
     * functionality within the library to multiply the Laplacian by the
     * appropriate constants for all diagonal terms, and off-diagonal terms use
     * LinearElasticSystem::BuildLaplacianIJMatrix to construct matrices
     * containing
     *
     * \f[ \int \partial_{x_i} \phi_k \partial_{x_j} \phi_l \f]
     *
     * where mixed derivative terms are present. Symmetry (in terms of k,l) is
     * exploited to avoid constructing this matrix repeatedly.
     *
     * @todo Make static condensation optional and construct full system
     *       instead.
     */
    void LinearElasticSystem::BuildMatrixSystem()
    {
        const int nEl = m_fields[0]->GetExpSize();
        const int nVel = m_fields[0]->GetCoordim(0);
        
        LocalRegions::ExpansionSharedPtr exp;
        int n;

        // Factors map for matrix keys.
        StdRegions::ConstFactorMap factors;

        // Calculate various constants

        // mu
        NekDouble a = m_E * (1.0 - m_nu) / (1.0 + m_nu) / (1.0 - 2.0*m_nu);
        // lambda
        NekDouble b = m_E * 0.5 / (1.0 + m_nu);
        // lambda + mu?
        NekDouble c = m_E * m_nu / (1.0 + m_nu) / (1.0 - 2.0*m_nu);

        bool verbose = m_session->DefinesCmdLineArgument("verbose");
        bool root    = m_session->GetComm()->GetRank() == 0;

        // Loop over each element and construct matrices.
        if (nVel == 2)
        {
            for (n = 0; n < nEl; ++n)
            {
                exp = m_fields[0]->GetExp(m_fields[0]->GetOffset_Elmt_Id(n));
                const int nPhys = exp->GetTotPoints();
                Array<OneD, NekDouble> aArr(nPhys, a);
                Array<OneD, NekDouble> bArr(nPhys, b);
                
                StdRegions::VarCoeffMap varcoeffA, varcoeffD;
                varcoeffA[StdRegions::eVarCoeffD00] = aArr;
                varcoeffA[StdRegions::eVarCoeffD11] = bArr;
                varcoeffD[StdRegions::eVarCoeffD00] = bArr;
                varcoeffD[StdRegions::eVarCoeffD11] = aArr;

                LocalRegions::MatrixKey matkeyA(StdRegions::eLaplacian,
                                                exp->DetShapeType(),
                                                *exp, factors, varcoeffA);
                LocalRegions::MatrixKey matkeyD(StdRegions::eLaplacian,
                                                exp->DetShapeType(),
                                                *exp, factors, varcoeffD);

                /*
                 * mat holds the linear operator [ A B ] acting on [ u ].
                 *                               [ C D ]           [ v ]
                 */
                Array<TwoD, DNekMatSharedPtr> mat(2,2);
                mat[0][0] = exp->GenMatrix(matkeyA);
                mat[0][1] = BuildLaplacianIJMatrix(1, 0, c, exp);
                mat[1][0] = mat[0][1];
                mat[1][1] = exp->GenMatrix(matkeyD);

                // Set up the statically condensed block for this element.
                SetStaticCondBlock(n, exp, mat);

                if (verbose && root)
                {
                    cout << "\rBuilding matrix system: "
                         << (int)(100.0 * n / nEl) << "%" << flush;
                }
            }
        }
        else if (nVel == 3)
        {
            for (n = 0; n < nEl; ++n)
            {
                exp = m_fields[0]->GetExp(m_fields[0]->GetOffset_Elmt_Id(n));
                const int nPhys = exp->GetTotPoints();
                Array<OneD, NekDouble> aArr(nPhys, a);
                Array<OneD, NekDouble> bArr(nPhys, b);

                StdRegions::VarCoeffMap varcoeffA, varcoeffE, varcoeffI;
                varcoeffA[StdRegions::eVarCoeffD00] = aArr;
                varcoeffA[StdRegions::eVarCoeffD11] = bArr;
                varcoeffA[StdRegions::eVarCoeffD22] = bArr;
                varcoeffE[StdRegions::eVarCoeffD00] = bArr;
                varcoeffE[StdRegions::eVarCoeffD11] = aArr;
                varcoeffE[StdRegions::eVarCoeffD22] = bArr;
                varcoeffI[StdRegions::eVarCoeffD00] = bArr;
                varcoeffI[StdRegions::eVarCoeffD11] = bArr;
                varcoeffI[StdRegions::eVarCoeffD22] = aArr;

                LocalRegions::MatrixKey matkeyA(StdRegions::eLaplacian,
                                                exp->DetShapeType(),
                                                *exp, factors, varcoeffA);
                LocalRegions::MatrixKey matkeyE(StdRegions::eLaplacian,
                                                exp->DetShapeType(),
                                                *exp, factors, varcoeffE);
                LocalRegions::MatrixKey matkeyI(StdRegions::eLaplacian,
                                                exp->DetShapeType(),
                                                *exp, factors, varcoeffI);

                /*
                 * mat holds the linear operator [ A B C ] acting on [ u ].
                 *                               [ D E F ]           [ v ]
                 *                               [ G H I ]           [ w ]
                 */
                Array<TwoD, DNekMatSharedPtr> mat(3,3);
                mat[0][0] = exp->GenMatrix(matkeyA);
                mat[0][1] = BuildLaplacianIJMatrix(1, 0, c, exp);
                mat[0][2] = BuildLaplacianIJMatrix(2, 0, c, exp);
                
                mat[1][0] = mat[0][1];
                mat[1][1] = exp->GenMatrix(matkeyE);
                mat[1][2] = BuildLaplacianIJMatrix(2, 1, c, exp);
                
                mat[2][0] = mat[0][2];
                mat[2][1] = mat[1][2];
                mat[2][2] = exp->GenMatrix(matkeyI);
                
                // Set up the statically condensed block for this element.
                SetStaticCondBlock(n, exp, mat);

                if (verbose && root)
                {
                    cout << "\rBuilding matrix system: "
                         << (int)(100.0 * n / nEl) << "%" << flush;
                }
            }
        }

        if (verbose && root)
        {
            cout << "\rBuilding matrix system: done." << endl;
        }
    }

    /**
     * @brief Generate summary at runtime.
     */
    void LinearElasticSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        EquationSystem::SessionSummary(s);
    }

    /**
     * @brief Solve elliptic linear elastic system.
     *
     * The solve proceeds as follows:
     *
     * - Create a MultiRegions::GlobalLinSys object.
     * - Evaluate a forcing function from the session file.
     * - If the temperature term is enabled, evaluate this and add to forcing.
     * - Apply Dirichlet boundary conditions.
     * - Scatter forcing into the correct ordering according to the coupled
     *   assembly map.
     * - Do the solve.
     * - Scatter solution back to fields and backwards transform to physical
     *   space.
     */
    void LinearElasticSystem::v_DoSolve()
    {
        
        int i, j, k, l, m, nv;
        const int nVel = m_fields[0]->GetCoordim(0);

        // Build initial matrix system.
        BuildMatrixSystem();

        // Now we've got the matrix system set up, create a GlobalLinSys object.
        MultiRegions::GlobalLinSysKey key(
            StdRegions::eLinearAdvectionReaction, m_assemblyMap);        
        MultiRegions::GlobalLinSysSharedPtr linSys;

        // Currently either direct or iterative static condensation is
        // supported.
        if (m_assemblyMap->GetGlobalSysSolnType() == MultiRegions::eDirectStaticCond)
        {
            linSys = MemoryManager<
                MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(
                    key, m_fields[0], m_schurCompl, m_BinvD, m_C, m_Dinv,
                    m_assemblyMap);
        }
        else if (m_assemblyMap->GetGlobalSysSolnType() ==
                     MultiRegions::eIterativeStaticCond)
        {
            linSys = MemoryManager<
                MultiRegions::GlobalLinSysIterativeStaticCond>::AllocateSharedPtr(
                    key, m_fields[0], m_schurCompl, m_BinvD, m_C, m_Dinv,
                    m_assemblyMap, MultiRegions::NullPreconditionerSharedPtr);
        }

        const int nCoeffs = m_fields[0]->GetNcoeffs();
        const int nGlobDofs = m_assemblyMap->GetNumGlobalCoeffs();

        //
        // -- Evaluate forcing functions
        //

        // Evaluate the forcing function from the XML file.
        Array<OneD, Array<OneD, NekDouble> > forcing(nVel);
        EvaluateFunction(forcing, "Forcing");

        // Add temperature term
        string tempEval;
        m_session->LoadSolverInfo("Temperature", tempEval, "None");
        string refElement;
        m_session->LoadSolverInfo("RefElement", refElement, "None");
        
        if (tempEval == "Jacobian")
        {
            // Allocate storage
            m_temperature = Array<OneD, Array<OneD, NekDouble> >(nVel);

            for (nv = 0; nv < nVel; ++nv)
            {
                Array<OneD, NekDouble> tmp;
                m_temperature[nv] = Array<OneD, NekDouble>(
                    m_fields[nv]->GetNpoints());
                
                for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
                {
                    // Calculate element area
                    LocalRegions::ExpansionSharedPtr exp =
                        m_fields[0]->GetExp(i);
                    LibUtilities::PointsKeyVector pkey =
                        exp->GetPointsKeys();
                    Array<OneD, NekDouble> jac =
                        exp->GetMetricInfo()->GetJac(pkey);

                    int offset = m_fields[0]->GetPhys_Offset(i);

                    if (exp->GetMetricInfo()->GetGtype() ==
                        SpatialDomains::eDeformed)
                    {
                        Vmath::Smul(exp->GetTotPoints(), m_beta, jac, 1,
                                    tmp = m_temperature[nv] + offset, 1);
                    }
                    else
                    {
                        Vmath::Fill(exp->GetTotPoints(), m_beta*jac[0],
                                    tmp = m_temperature[nv] + offset, 1);
                    }
                }
                m_fields[nv]->PhysDeriv(nv, m_temperature[nv], forcing[nv]);
                
            }
        }
        if (tempEval == "Metric")
        {
            m_temperature = Array<OneD, Array<OneD, NekDouble> >(nVel);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmpstress(nVel);

            for (nv = 0; nv < nVel; ++nv)
            {
                m_temperature[nv] = Array<OneD, NekDouble>(
                    m_fields[nv]->GetNpoints());
                tmpstress[nv] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                for (i = 0; i < nVel; ++i)
                {
                    tmpstress[nv][i] = Array<OneD, NekDouble>(
                        m_fields[nv]->GetNpoints());
                }
            }

            for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
            {
                // Calculate element area
                LocalRegions::ExpansionSharedPtr exp =
                    m_fields[0]->GetExp(i);
                LibUtilities::PointsKeyVector pkey = exp->GetPointsKeys();
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > deriv
                    = exp->GetMetricInfo()->GetDeriv(pkey);
                int offset = m_fields[0]->GetPhys_Offset(i);
                
                // Compute metric tensor
                Array<OneD, NekDouble> tmp(nVel*nVel, 0.0);
                for (j = 0; j < deriv[0][0].num_elements(); ++j)
                {
                    // Setting up the metric tensor
                    for (k = 0; k < nVel; ++k)
                    {
                        for (l = 0; l < nVel; ++l)
                        {
                            tmp[k*nVel+l] = 0.0;
                            
                            for (int m = 0; m < nVel; ++m)
                            {
                                tmp[k*nVel+l] += deriv[k][m][j] * deriv[l][m][j];
                            }
                        }
                    }
                
                    // Compute eigenvalues/eigenvectors.
                    char jobvl = 'N', jobvr = 'V';
                    int worklen = 8*nVel, info;
                    
                    DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
                    DNekMat evec   (nVel, nVel, 0.0, eFULL);
                    DNekMat evecinv(nVel, nVel, 0.0, eFULL);
                    Array<OneD, NekDouble> vl  (nVel*nVel);
                    Array<OneD, NekDouble> work(worklen);
                    Array<OneD, NekDouble> wi  (nVel);
                    
                    Lapack::Dgeev(jobvl, jobvr, nVel, &tmp[0], nVel,
                                  &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                                  &(evec.GetPtr())[0], nVel,
                                  &work[0], worklen, info);
                    
                    evecinv = evec;
                    evecinv.Invert();

                    // rescaling of the eigenvalues
                    for (nv = 0; nv < nVel; ++nv)
                    {
                        eval(nv,nv) = m_beta * (sqrt(eval(nv,nv))-1.0);
                    }

                    DNekMat beta = evec * eval * evecinv;
                    
                    NekDouble term = 0.0;

                    tmpstress[0][0][offset+j] = -beta(0,0);
                    tmpstress[1][0][offset+j] = -beta(1,0);
                    tmpstress[0][1][offset+j] = -beta(0,1);
                    tmpstress[1][1][offset+j] = -beta(1,1);
                }

                if (deriv[0][0].num_elements() != exp->GetTotPoints())
                {
                    Array<OneD, NekDouble> tmp;
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[0][0][offset],
                                tmp = tmpstress[0][0] + offset, 1);
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[1][0][offset],
                                tmp = tmpstress[1][0] + offset, 1);
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[0][1][offset],
                                tmp = tmpstress[0][1] + offset, 1);
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[1][1][offset],
                                tmp = tmpstress[1][1] + offset, 1);
                }
            }

            Array<OneD, NekDouble> tmpderiv(m_fields[0]->GetNpoints());
            for (nv = 0; nv < nVel; ++nv)
            {
                m_fields[nv]->PhysDeriv(0, tmpstress[nv][0], tmpderiv);
                m_fields[nv]->PhysDeriv(1, tmpstress[nv][1], m_temperature[nv]);
                Vmath::Vadd(m_fields[nv]->GetNpoints(), tmpderiv, 1, m_temperature[nv], 1, m_temperature[nv], 1);
                Vmath::Vcopy(m_fields[nv]->GetNpoints(), m_temperature[nv], 1, forcing[nv], 1);
            }
        }
        if (tempEval == "LinearIdealMetric")
        {
            m_temperature = Array<OneD, Array<OneD, NekDouble> >(nVel);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmpstress(nVel);

            for (nv = 0; nv < nVel; ++nv)
            {
                m_temperature[nv] = Array<OneD, NekDouble>(
                    m_fields[nv]->GetNpoints());
                tmpstress[nv] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                for (i = 0; i < nVel; ++i)
                {
                    tmpstress[nv][i] = Array<OneD, NekDouble>(
                        m_fields[nv]->GetNpoints());
                }
            }

            // Grab existing basis keys for first element.
            LibUtilities::BasisKey bkey0 =
                m_fields[0]->GetExp(0)->GetBasis(0)->GetBasisKey();
            LibUtilities::BasisKey bkey1 =
                m_fields[0]->GetExp(0)->GetBasis(1)->GetBasisKey();
            
            for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
            {
                // Calculate element area
                LocalRegions::ExpansionSharedPtr exp =
                    m_fields[0]->GetExp(i);
                
                int offset = m_fields[0]->GetPhys_Offset(i);
                int nvert = exp->GetGeom()->GetNumVerts();
                int nqel = exp->GetTotPoints();
                
                // Create metric for ref to element
                // hence, g = dxi/dx
                
                LibUtilities::PointsKeyVector pkey = exp->GetPointsKeys();
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > deriv
                = exp->GetMetricInfo()->GetDeriv(pkey);
                
                Array<OneD, NekDouble> tmp_met(nVel*nVel, 0.0);
                for (int p = 0; p < deriv[0][0].num_elements(); ++p)
                {
                    // Setting up the metric tensor
                    for (int k = 0; k < nVel; ++k)
                    {
                        for (int l = 0; l < nVel; ++l)
                        {
                            tmp_met[k*nVel+l] = 0.0;
                            
                            for (int m = 0; m < nVel; ++m)
                            {
                                tmp_met[k*nVel+l] += deriv[k][m][p] * deriv[l][m][p];
                            }
                        }
                    }
                }
                
                // Ignore linear elements
                if (deriv[0][0].num_elements() == 1)
                {
                    Array<OneD, NekDouble> tmp;
                    for (nv = 0; nv < nVel; ++nv)
                    {
                        Vmath::Zero(exp->GetTotPoints(),
                                    tmp = m_temperature[nv] + offset, 1);
                    }
                    continue;
                }

                // Get distribution of nodal points
                
                Array<OneD, NekDouble> tmp(nVel*nVel, 0.0);
                
                
                Array<OneD, NekDouble> xVertLIN(nvert, 0.0);
                Array<OneD, NekDouble> yVertLIN(nvert, 0.0);
                Array<OneD, NekDouble> zVertLIN(nvert, 0.0);
                
                // Determine mapping of the linear (ideal) element to the reference element
                
                for (j = 0; j < nvert; ++j)
                {
                    m_fields[0]->GetExp(i)->GetGeom()->
                           GetVertex(j)->GetCoords(xVertLIN[j],
                                                   yVertLIN[j],
                                                   zVertLIN[j]);
                    
                }
                
                // I THINK THIS IS NOT NECESSARY IN THE END
                
                // Create linear element ONLY IMPLEMENTED FOR TRIANGLES STILL NEEDS TO BE GENERALISED.
                /*
                const int zero = 0;
                const int one=1;
                const int two=2;
                const double dZero = 0.0;
                
                StdRegions::Orientation edgeDir = StdRegions::eForwards;
                SpatialDomains::PointGeomSharedPtr verts[3];
                
                verts[0] = MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(two,zero,xVertLIN[0],yVertLIN[0],dZero);
                verts[1] = MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(two,one,xVertLIN[1],yVertLIN[1],dZero);
                verts[2] = MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(two,two,xVertLIN[2],yVertLIN[2],dZero);
                
                // Set up Edges
                SpatialDomains::SegGeomSharedPtr edges[3];
                edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,verts[0],verts[1]);
                edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one,verts[1],verts[2]);
                edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two,verts[2],verts[0]);
            
                StdRegions::Orientation eorient[3];
                
                eorient[0] = edgeDir;
                eorient[1] = edgeDir;
                eorient[2] = edgeDir;
                
                SpatialDomains::TriGeomSharedPtr geom = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(zero,verts,edges,eorient);
                
                geom->SetOwnData();
                
                LocalRegions::NodalTriExpSharedPtr triNodal =
                MemoryManager<LocalRegions::NodalTriExp>
                ::AllocateSharedPtr(bkey0, bkey1, LibUtilities::eNodalTriElec, geom);
                */
                
                Array<OneD, NekDouble> xf(3.0, 0.0);
                Array<OneD, NekDouble> yf(3.0, 0.0);
                
                MappingIdealToRef(xVertLIN,
                                  yVertLIN,
                                  xf,
                                  yf);
                
                DNekMat Jac   (nVel, nVel, 0.0, eFULL);
                
                Jac(0,0) = xf[0];Jac(0,1) = xf[1];
                Jac(1,0) = yf[0];Jac(1,1) = yf[1];
                
                for (j = 0; j < nqel; ++j)
                {
                    // Using the ideal to ref mapping functions, determine the coordinates in the reference element
                
                    Array<OneD, NekDouble> Jnew(nVel*nVel,0.0);
                    Array<OneD, NekDouble> Gnew(nVel*nVel,0.0);
    
                    // Apply chain rule so,
                    // Jnew = dgi_detaj = dgi_dxi1*dxi1_detaj + dgi_dxi2*dxi2_detaj
                    Jnew[0] = Jac(0,0)*tmp_met[0] + Jac(0,1)*tmp_met[1];
                    Jnew[1] = Jac(1,0)*tmp_met[0] + Jac(1,1)*tmp_met[1];
                    Jnew[2] = Jac(0,0)*tmp_met[2] + Jac(0,1)*tmp_met[3];
                    Jnew[3] = Jac(1,0)*tmp_met[2] + Jac(1,1)*tmp_met[3];
                    
                    // Determine the new metric tensor as Gnew = Jnew*Jnew^t
                    Gnew[0] = Jnew[0]*Jnew[0] + Jnew[1]*Jnew[1];
                    Gnew[1] = Jnew[0]*Jnew[2] + Jnew[1]*Jnew[3];
                    Gnew[2] = Jnew[0]*Jnew[2] + Jnew[1]*Jnew[3];
                    Gnew[3] = Jnew[2]*Jnew[2] + Jnew[3]*Jnew[3];
                
                    // Compute eigenvalues/eigenvectors.
                    char jobvl = 'N', jobvr = 'V';
                    
                    int worklen = 8*nVel, info;
                    
                    DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
                    DNekMat evec   (nVel, nVel, 0.0, eFULL);
                    DNekMat evecinv(nVel, nVel, 0.0, eFULL);
                    
                    Array<OneD, NekDouble> vl  (nVel*nVel);
                    Array<OneD, NekDouble> work(worklen);
                    Array<OneD, NekDouble> wi  (nVel);
                    
                    Lapack::Dgeev(jobvl, jobvr, nVel, &Gnew[0], nVel,
                                  &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                                  &(evec.GetPtr())[0], nVel,
                                  &work[0], worklen, info);
                    
                    evecinv = evec;
                    evecinv.Invert();
                    
                    // rescaling of the eigenvalues
                    for (nv = 0; nv < nVel; ++nv)
                    {
                        eval(nv,nv) = sqrt(eval(nv,nv))-1;
                    }

                    
                    DNekMat beta = evec * eval * evecinv;
                    
                    NekDouble term = 0.0;
                
                    tmpstress[0][0][offset+j] = -beta(0,0);
                    tmpstress[1][0][offset+j] = -beta(0,1);
                    tmpstress[0][1][offset+j] = -beta(1,0);
                    tmpstress[1][1][offset+j] = -beta(1,1);
                }
                
                cout << endl;
                Array<OneD, NekDouble> tmp2;
                
                if (deriv[0][0].num_elements() != exp->GetTotPoints())
                {
                    Array<OneD, NekDouble> tmp3;
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[0][0][offset],
                                tmp3 = tmpstress[0][0] + offset, 1);
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[1][0][offset],
                                tmp3 = tmpstress[1][0] + offset, 1);
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[0][1][offset],
                                tmp3 = tmpstress[0][1] + offset, 1);
                    Vmath::Fill(exp->GetTotPoints(), tmpstress[1][1][offset],
                                tmp3 = tmpstress[1][1] + offset, 1);
                }
            }
            
            Array<OneD, NekDouble> tmpderiv(m_fields[0]->GetNpoints());
            
            for (nv = 0; nv < nVel; ++nv)
            {

                m_fields[nv]->PhysDeriv(0,
                                        tmpstress[nv][0],
                                        tmpderiv);
                
                m_fields[nv]->PhysDeriv(1,
                                        tmpstress[nv][1],
                                        m_temperature[nv]);
                
                Vmath::Vadd(m_fields[nv]->GetNpoints(),
                            tmpderiv, 1,
                            m_temperature[nv], 1,
                            m_temperature[nv], 1);
                
                Vmath::Vcopy(m_fields[nv]->GetNpoints(),
                             m_temperature[nv], 1,
                             forcing[nv], 1);
            }
        }
        else if (tempEval == "None")
        {
            ASSERTL0(false, "Unknown temperature form: " + tempEval);
        }


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
        Array<OneD, NekDouble> inout    (nGlobDofs, 0.0);
        Array<OneD, NekDouble> rhs      (nGlobDofs, 0.0);

        for (nv = 0; nv < nVel; ++nv)
        {
            // Take the inner product of the forcing function.
            Array<OneD, NekDouble> tmp(nCoeffs);
            m_fields[nv]->IProductWRTBase_IterPerExp(forcing[nv], tmp);

            // Scatter forcing into RHS vector according to the ordering
            // dictated in the comment above.
            for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
            {
                int nBnd   = m_bmap[i].num_elements();
                int nInt   = m_imap[i].num_elements();
                int offset = m_fields[nv]->GetCoeff_Offset(i);

                for (j = 0; j < nBnd; ++j)
                {
                    forCoeffs[nVel*offset + nv*nBnd + j] = tmp[offset+m_bmap[i][j]];
                }
                for (j = 0; j < nInt; ++j)
                {
                    forCoeffs[nVel*(offset + nBnd) + nv*nInt + j] = tmp[offset+m_imap[i][j]];
                }
            }
        }

        // -- Impose Dirichlet boundary conditions.

        // First try to do parallel assembly: the intention here is that
        // Dirichlet values at some edges/vertices need to be communicated to
        // processes which don't contain the entire boundary region. See
        // ContField2D::v_ImposeDirichletConditions for more detail.
        map<int, vector<MultiRegions::ExtraDirDof> > &extraDirDofs =
            m_assemblyMap->GetExtraDirDofs();
        map<int, vector<MultiRegions::ExtraDirDof> >::iterator it;

        for (nv = 0; nv < nVel; ++nv)
        {
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp
                = m_fields[nv]->GetBndCondExpansions();

            // First try to do parallel stuff
            for (it = extraDirDofs.begin(); it != extraDirDofs.end(); ++it)
            {
                for (i = 0; i < it->second.size(); ++i)
                {
                    inout[it->second.at(i).get<1>()*nVel + nv] =
                        bndCondExp[it->first]->GetCoeffs()[
                            it->second.at(i).get<0>()]*it->second.at(i).get<2>();
                }
            }
        }

        m_assemblyMap->UniversalAssemble(inout);

        // Counter for the local Dirichlet boundary to global ordering.
        int bndcnt = 0;

        // Now assemble local boundary contributions.
        for (nv = 0; nv < nVel; ++nv)
        {
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

        //
        // -- Perform solve
        //

        // Assemble forcing into the RHS.
        m_assemblyMap->Assemble(forCoeffs, rhs);

        // Negate RHS to be consistent with matrix definition.
        Vmath::Neg(rhs.num_elements(), rhs, 1);

        // Solve.
        linSys->Solve(rhs, inout, m_assemblyMap);

        // Scatter the global ordering back to the alternate local ordering.
        Array<OneD, NekDouble> tmp(nVel * nCoeffs);
        m_assemblyMap->GlobalToLocal(inout, tmp);

        //
        // -- Postprocess
        //

        // Scatter back to field degrees of freedom
        for (nv = 0; nv < nVel; ++nv)
        {
            for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
            {
                int nBnd   = m_bmap[i].num_elements();
                int nInt   = m_imap[i].num_elements();
                int offset = m_fields[nv]->GetCoeff_Offset(i);

                for (j = 0; j < nBnd; ++j)
                {
                    m_fields[nv]->UpdateCoeffs()[offset+m_bmap[i][j]] =
                        tmp[nVel*offset + nv*nBnd + j];
                }
                for (j = 0; j < nInt; ++j)
                {
                    m_fields[nv]->UpdateCoeffs()[offset+m_imap[i][j]] =
                        tmp[nVel*(offset + nBnd) + nv*nInt + j];
                }
            }
            m_fields[nv]->BwdTrans(m_fields[nv]->GetCoeffs(),
                                   m_fields[nv]->UpdatePhys());
        }
    }

    /**
     * @brief Given a block matrix for an element, construct its static
     * condensation matrices.
     *
     * This routine essentially duplicates the logic present in many of the
     * LocalRegions matrix generation routines to construct the statically
     * condensed equivalent of mat to pass to the GlobalLinSys solver.
     *
     * @param n    Element/block number.
     * @param exp  Pointer to expansion.
     * @param mat  Block matrix containing matrix operator.
     */
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
        const MatrixStorage s = eFULL; // Maybe look into doing symmetric
                                       // version of this?

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
                        (*A)(k + i*nB, l + j*nB) =
                            (*mat[i][j])(m_bmap[n][k], m_bmap[n][l]);
                    }

                    for (l = 0; l < nI; ++l)
                    {
                        (*B)(k + i*nB, l + j*nI) =
                            (*mat[i][j])(m_bmap[n][k], m_imap[n][l]);
                    }
                }

                // Interior-boundary / interior-interior
                for (k = 0; k < nI; ++k)
                {
                    for (l = 0; l < nB; ++l)
                    {
                        (*C)(k + i*nI, l + j*nB) =
                            (*mat[i][j])(m_imap[n][k], m_bmap[n][l]);
                    }

                    for (l = 0; l < nI; ++l)
                    {
                        (*D)(k + i*nI, l + j*nI) =
                            (*mat[i][j])(m_imap[n][k], m_imap[n][l]);
                    }
                }
            }
        }

        // Construct static condensation matrices.
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

    /**
     * @brief Construct a LaplacianIJ matrix for a given expansion.
     *
     * This routine constructs matrices whose entries contain the evaluation of
     *
     * \f[ \partial_{\tt k1} \phi_i \partial_{\tt k2} \phi_j \,dx  \f]
     */
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

    void LinearElasticSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        const int nVel    = m_fields[0]->GetCoordim(0);
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        static char *dimStr[] = { "X", "Y", "Z" };
        
        if (m_temperature.num_elements() == 0)
        {
            return;
        }

        for (int i = 0; i < nVel; ++i)
        {
            Array<OneD, NekDouble> tFwd(nCoeffs);
            m_fields[i]->FwdTrans_IterPerExp(m_temperature[i], tFwd);
            fieldcoeffs.push_back(tFwd);
            variables.push_back("Temp" + boost::lexical_cast<std::string>(dimStr[i]));
        }
    }
}
