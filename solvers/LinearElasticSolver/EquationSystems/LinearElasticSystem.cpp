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

    /* 
     * @brief Generate mapping between a triangle and the reference element.
     *
     * Mapping that requires the coordinates of the three vertices (x,y) of the
     * triangular element and outputs the function coefficients that defines the
     * mapping to the reference element for each coordinate (xfunc and yfunc)
     *
     * @param x  x co-ordinates of triangle
     * @param y  y co-ordinates of triangle
     * @param xf  x co-ordinates of mapping
     * @param yf  y co-ordinates of mapping
     */
    /*inline void MappingIdealToRef(Array<OneD, NekDouble> x,
                                  Array<OneD, NekDouble> y,
                                  Array<OneD, NekDouble> xf,
                                  Array<OneD, NekDouble> yf)
    {
        int n = x.num_elements();
        
        DNekMat map(n, n, 0.0, eFULL);
        DNekMat mapinv(n, n, 0.0, eFULL);
        cout << "NEW SET OF POITS" << endl;
        cout << endl;
        cout << " x(1)=" << x[0] << "; x(2)=" << x[1] << "; x(3)=" << x[2] << endl;
        cout << " y(1)=" << y[0] << "; y(2)=" << y[1] << "; y(3)=" << y[2] << endl;
        
        map(0,0) = x[0];map(0,1) = y[0];map(0,2) = 1.0;
        map(1,0) = x[1];map(1,1) = y[1];map(1,2) = 1.0;
        map(2,0) = x[2];map(2,1) = y[2];map(2,2) = 1.0;
        
        cout << endl;
        
        mapinv = map;
        mapinv.Invert();
        
        Array<OneD, NekDouble> xref (n,0.0), yref (n,0.0);
        
        xref[0] = -1.0;xref[1] =  1.0;xref[2] = -1.0;
        yref[0] = -1.0;yref[1] = -1.0;yref[2] =  1.0;
        
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                xf[i] += mapinv(i,j)*xref[j];
                yf[i] += mapinv(i,j)*yref[j];
                
            }
            cout << "xf(" << i+1 << ")=" << xf[i] << ";"
               << " " << "yf(" << i+1 << ")=" << yf[i] << ";" << endl;
        }
    }*/
    
    inline void MappingIdealToRef(Array<OneD, NekDouble> x,
                                  Array<OneD, NekDouble> y,
                                  Array<OneD, NekDouble> xf,
                                  Array<OneD, NekDouble> yf)
    {
        int n = x.num_elements();
        
        DNekMat map(n, n, 0.0, eFULL);
        DNekMat mapref(n, n, 0.0, eFULL);
        DNekMat mapinv(n, n, 0.0, eFULL);
        DNekMat newmap(n, n, 0.0, eFULL);
        
        map(0,0) = x[0];map(1,0) = y[0];map(2,0) = 1.0;
        map(0,1) = x[1];map(1,1) = y[1];map(2,1) = 1.0;
        map(0,2) = x[2];map(1,2) = y[2];map(2,2) = 1.0;
        
        mapref(0,0) = -1.0;mapref(1,0) = -1.0;mapref(2,0) = 1.0;
        mapref(0,1) =  1.0;mapref(1,1) = -1.0;mapref(2,1) = 1.0;
        mapref(0,2) = -1.0;mapref(1,2) =  1.0;mapref(2,2) = 1.0;
        
        mapinv = map;
        mapinv.Invert();
        
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                NekDouble tmp = 0;
                
                for (int k = 0; k < n; ++k)
                {
                    tmp += mapref(i,k)*mapinv(k,j);
                }
                
                newmap(i,j) = tmp;
            }
        }
        
        xf[0] = newmap(0,0);xf[1] = newmap(0,1);xf[2] = newmap(0,2);
        yf[0] = newmap(1,0);yf[1] = newmap(1,1);yf[2] = newmap(1,2);
    }
    

    /**
     * @brief Default constructor.
     */
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
        m_session->LoadParameter("E",    m_E,    1.00);
        m_session->LoadParameter("nu",   m_nu,   0.25);
        m_session->LoadParameter("Beta", m_beta, 1.00);
        
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

        AddSummaryItem(s, "Young's modulus", m_E);
        AddSummaryItem(s, "Poisson ratio", m_nu);
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
        
        int i, j, k, l, nv;
        const int nVel = m_fields[0]->GetCoordim(0);

        // Build initial matrix system.
        BuildMatrixSystem();

        // Now we've got the matrix system set up, create a GlobalLinSys
        // object. We mask ourselves as LinearAdvectionReaction to create a full
        // matrix instead of symmetric storage.
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

        linSys->Initialise(m_assemblyMap);

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
        else if (tempEval == "Metric")
        {
            ASSERTL0(m_fields[0]->GetCoordim(0) == 2 &&
                     m_graph->GetAllQuadGeoms().size() == 0,
                     "LinearIdealMetric temperature only implemented for "
                     "two-dimensional triangular meshes.");

            m_temperature = Array<OneD, Array<OneD, NekDouble> >(nVel);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmpstress(nVel);
            
            m_Eig = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(nVel);
            
            
            for (nv = 0; nv < nVel; ++nv)
            {
                m_temperature[nv] = Array<OneD, NekDouble>(
                    m_fields[nv]->GetNpoints());
                
                m_Eig[nv] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                tmpstress[nv] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                for (i = 0; i < nVel; ++i)
                {
                    tmpstress[nv][i] = Array<OneD, NekDouble>(
                        m_fields[nv]->GetNpoints());
                    
                    m_Eig[nv][i] = Array<OneD, NekDouble>(
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

                Array<OneD, NekDouble> xIdeal(3), yIdeal(3);
                xIdeal[0] = exp->GetGeom()->GetVertex(0)->x();
                xIdeal[1] = exp->GetGeom()->GetVertex(1)->x();
                xIdeal[2] = exp->GetGeom()->GetVertex(2)->x();
                yIdeal[0] = exp->GetGeom()->GetVertex(0)->y();
                yIdeal[1] = exp->GetGeom()->GetVertex(1)->y();
                yIdeal[2] = exp->GetGeom()->GetVertex(2)->y();

                Array<OneD, NekDouble> i2rx(3, 0.0), i2ry(3, 0.0);
                MappingIdealToRef(xIdeal, yIdeal, i2rx, i2ry);

                // Compute metric tensor
                Array<OneD, NekDouble> tmp (nVel*nVel, 0.0);
                Array<TwoD, NekDouble> tmp2(nVel, nVel);

                for (j = 0; j < deriv[0][0].num_elements(); ++j)
                {
                    // Setting up the metric tensor
                    NekDouble dxi1dx1 = deriv[0][0][j];
                    NekDouble dxi2dx1 = deriv[0][1][j];
                    NekDouble dxi1dx2 = deriv[1][0][j];
                    NekDouble dxi2dx2 = deriv[1][1][j];

                    tmp2[0][0] = dxi1dx1 * i2rx[0] + dxi1dx2 * i2ry[0];
                    tmp2[0][1] = dxi1dx1 * i2rx[1] + dxi1dx2 * i2ry[1];
                    tmp2[1][0] = dxi2dx1 * i2rx[0] + dxi2dx2 * i2ry[0];
                    tmp2[1][1] = dxi2dx1 * i2rx[1] + dxi2dx2 * i2ry[1];


                    for (k = 0; k < nVel; ++k)
                    {
                        for (l = 0; l < nVel; ++l)
                        {
                            tmp[k*nVel+l] = 0.0;
                            
                            for (int m = 0; m < nVel; ++m)
                            {
                                tmp[k*nVel+l] += tmp2[k][m] * tmp2[l][m];
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
                        eval(nv,nv) = m_beta * (sqrt(eval(nv,nv)) - 1.0);
                    }

                    DNekMat beta = evec * eval * evecinv;
                    
                    tmpstress[0][0][offset+j] = beta(0,0);
                    tmpstress[1][0][offset+j] = beta(1,0);
                    tmpstress[0][1][offset+j] = beta(0,1);
                    tmpstress[1][1][offset+j] = beta(1,1);
                    
                    /*
                    m_Eig[0][0][offset+j] = evec(0,0);
                    m_Eig[0][1][offset+j] = evec(0,1);
                    m_Eig[1][0][offset+j] = evec(1,0);
                    m_Eig[1][1][offset+j] = evec(1,1);
                    */
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
            
            Vmath::Vcopy(m_fields[0]->GetNpoints(), tmpstress[0][0], 1, m_Eig[0][0], 1);
            Vmath::Vcopy(m_fields[0]->GetNpoints(), tmpstress[0][1], 1, m_Eig[0][1], 1);
            Vmath::Vcopy(m_fields[0]->GetNpoints(), tmpstress[1][0], 1, m_Eig[1][0], 1);
            Vmath::Vcopy(m_fields[0]->GetNpoints(), tmpstress[1][1], 1, m_Eig[1][1], 1);

            Array<OneD, NekDouble> tmpderiv(m_fields[0]->GetNpoints());
            for (nv = 0; nv < nVel; ++nv)
            {
                m_fields[nv]->PhysDeriv(0, tmpstress[nv][0], tmpderiv);
                m_fields[nv]->PhysDeriv(1, tmpstress[nv][1], m_temperature[nv]);
                
                Vmath::Vadd(m_fields[nv]->GetNpoints(), tmpderiv, 1, m_temperature[nv], 1, m_temperature[nv], 1);
                Vmath::Vcopy(m_fields[nv]->GetNpoints(), m_temperature[nv], 1, forcing[nv], 1);
            }
        }
        else if (tempEval == "LinearIdealMetric")
        {
            // Allocate storage
            m_temperature = Array<OneD, Array<OneD, NekDouble> >(nVel);
            
            for (nv = 0; nv < nVel; ++nv)
            {
                Array<OneD, NekDouble> tmp;
                m_temperature[nv] = Array<OneD, NekDouble>(
                                m_fields[nv]->GetNpoints());
                
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
                    LibUtilities::PointsKeyVector pkey =
                    exp->GetPointsKeys();
                    Array<OneD, NekDouble> jac =
                    exp->GetMetricInfo()->GetJac(pkey);
                    int offset = m_fields[0]->GetPhys_Offset(i);
                    
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > deriv
                    = exp->GetMetricInfo()->GetDeriv(pkey);
                    
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                                   der_n(nVel);
                    
                    for (int o = 0; o < nVel; ++o)
                    {
                        der_n[o] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                        
                        for (int p = 0; p < nVel; ++p)
                        {
                            der_n[o][p] = Array<OneD, NekDouble>(
                                            deriv[0][0].num_elements(),0.0);
                        }
                    }
                    
                    int nvert = exp->GetGeom()->GetNumVerts();
                    
                    Array<OneD, NekDouble> xVertLIN(3.0, 0.0);
                    Array<OneD, NekDouble> yVertLIN(3.0, 0.0);
                    Array<OneD, NekDouble> zVertLIN(3.0, 0.0);
                
                    Array<OneD, NekDouble> xf(3.0, 0.0);
                    Array<OneD, NekDouble> yf(3.0, 0.0);
                    
                    for (j = 0; j < nvert; ++j)
                    {
                        m_fields[0]->GetExp(i)->GetGeom()->
                        GetVertex(j)->GetCoords(xVertLIN[j],
                                                yVertLIN[j],
                                                zVertLIN[j]);
                    }

                    MappingIdealToRef(xVertLIN,
                                      yVertLIN,
                                      xf,
                                      yf);
                    
                    DNekMat JacNew   (nVel, nVel, 0.0, eFULL);
                    
                    JacNew(0,0) = xf[0]; JacNew(0,1) = xf[1];
                    JacNew(1,0) = yf[0]; JacNew(1,1) = yf[1];
                    
                    Array<OneD, NekDouble> jac2(deriv[0][0].num_elements(), 0.0);
                    
                    for (int k = 0; k < deriv[0][0].num_elements(); ++k)
                    {
                        JacNew(0,0) = xf[0]*deriv[0][0][k]+xf[1]*deriv[0][1][k];
                        JacNew(0,1) = yf[0]*deriv[0][0][k]+yf[1]*deriv[0][1][k];
                        
                        JacNew(1,0) = xf[0]*deriv[1][0][k]+xf[1]*deriv[1][1][k];
                        JacNew(1,1) = yf[0]*deriv[1][0][k]+yf[1]*deriv[1][1][k];
                        
                        der_n[0][0][k] = xf[0]*deriv[0][0][k]+xf[1]*deriv[0][1][k];
                        der_n[0][1][k] = yf[0]*deriv[0][0][k]+yf[1]*deriv[0][1][k];
                        
                        der_n[1][0][k] = xf[0]*deriv[1][0][k]+xf[1]*deriv[1][1][k];
                        der_n[1][1][k] = yf[0]*deriv[1][0][k]+yf[1]*deriv[1][1][k];
                        
                        Array<OneD, NekDouble> GMAT(nVel*nVel, 0.0);
                        
                        GMAT[0] = der_n[0][0][k]*der_n[0][0][k]
                                      + der_n[1][0][k]*der_n[1][0][k];
                        
                        GMAT[1] = der_n[0][0][k]*der_n[0][1][k]
                                        + der_n[1][0][k]*der_n[1][1][k];
                        
                        GMAT[2] = der_n[0][0][k]*der_n[0][1][k]
                                      + der_n[1][0][k]*der_n[1][1][k];
                        
                        GMAT[3] = der_n[1][1][k]*der_n[1][1][k]
                                      + der_n[0][1][k]*der_n[0][1][k];
                        
                        jac2[k] = JacNew(0,0)*JacNew(1,1)-JacNew(0,1)*JacNew(1,0);
                        
                        // Compute eigenvalues/eigenvectors.
                        char jobvl = 'N', jobvr = 'V';
                        int worklen = 8*nVel, info;
                        
                        DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
                        DNekMat evec   (nVel, nVel, 0.0, eFULL);
                        DNekMat evecinv(nVel, nVel, 0.0, eFULL);
                        Array<OneD, NekDouble> vl  (nVel*nVel);
                        Array<OneD, NekDouble> work(worklen);
                        Array<OneD, NekDouble> wi  (nVel);
                        
                        Lapack::Dgeev(jobvl, jobvr, nVel, &GMAT[0], nVel,
                                      &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                                      &(evec.GetPtr())[0], nVel,
                                      &work[0], worklen, info);
                        
                        evecinv = evec;
                        evecinv.Invert();
                        
                        // rescaling of the eigenvalues
                        for (nv = 0; nv < nVel; ++nv)
                        {
                            eval(nv,nv) = (sqrt(eval(nv,nv))-1.0);
                        }
                        
                        DNekMat beta = evec * eval * evecinv;
                        
                        tmpstress[0][0][offset+k] = -beta(0,0);
                        tmpstress[1][0][offset+k] = -beta(1,0);
                        tmpstress[0][1][offset+k] = -beta(0,1);
                        tmpstress[1][1][offset+k] = -beta(1,1);
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
                Array<OneD, NekDouble> tmpderiv2(m_fields[0]->GetNpoints());
                
                Vmath::Vadd(m_fields[0]->GetNpoints(),
                            tmpstress[1][0], 1,
                            tmpstress[1][0], 1,
                            m_temperature[0], 1);
                
                Vmath::Vadd(m_fields[0]->GetNpoints(),
                            tmpstress[0][1], 1,
                            tmpstress[0][1], 1,
                            m_temperature[1], 1);
                
                Vmath::Smul(m_fields[0]->GetNpoints(), m_beta, m_temperature[0], 1, forcing[0], 1);
                
                Vmath::Smul(m_fields[0]->GetNpoints(), m_beta, m_temperature[1], 1, forcing[1], 1);
                
                /*for (nv = 0; nv < nVel; ++nv)
                {
                    m_fields[nv]->PhysDeriv(0, tmpstress[0][nv], tmpderiv);
                    m_fields[nv]->PhysDeriv(1, tmpstress[1][nv], tmpderiv2);
                    Vmath::Vadd(m_fields[nv]->GetNpoints(), tmpstress[0][nv], 1, tmpstress[1][nv], 1, m_temperature[nv], 1);
                    
                    //Vmath::Vadd(m_fields[nv]->GetNpoints(), tmpderiv, 1, tmpderiv2, 1, m_temperature[nv], 1);
                    
                    //cout << Vmath::Vmax(m_temperature[nv].num_elements(),m_temperature[nv], 1) << endl;
                    //Vmath::Zero(m_temperature[nv].num_elements(),m_temperature[nv], 1);
                    Vmath::Smul(m_fields[nv]->GetNpoints(), m_beta, m_temperature[nv], 1, forcing[nv], 1);
                }*/
            }
        }
        else if (tempEval != "None")
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

            Vmath::Smul(
                nCoeffs, scale, &tmp1[0], 1, &(ret->GetPtr())[0]+i*nCoeffs, 1);
        }

        return ret;
    }

    void LinearElasticSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        const int nVel    = m_fields[0]->GetCoordim(0);
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        static string dimStr[4] = { "XX", "XY", "YX", "YY" };
       
        if (m_temperature.num_elements() == 0)
        {
            return;
        }
        
        Array<OneD, NekDouble> Eig1Comp1Fwd(nCoeffs);
        Array<OneD, NekDouble> Eig1Comp2Fwd(nCoeffs);
        Array<OneD, NekDouble> Eig2Comp1Fwd(nCoeffs);
        Array<OneD, NekDouble> Eig2Comp2Fwd(nCoeffs);
        
        m_fields[0]->FwdTrans(m_Eig[0][0], Eig1Comp1Fwd);
        fieldcoeffs.push_back(Eig1Comp1Fwd);
        
        m_fields[0]->FwdTrans(m_Eig[0][1], Eig1Comp2Fwd);
        fieldcoeffs.push_back(Eig1Comp2Fwd);
        
        m_fields[0]->FwdTrans(m_Eig[1][0], Eig2Comp1Fwd);
        fieldcoeffs.push_back(Eig2Comp1Fwd);
        
        m_fields[0]->FwdTrans(m_Eig[1][1], Eig2Comp2Fwd);
        fieldcoeffs.push_back(Eig2Comp2Fwd);
        
        for (int i = 0; i < nVel*2; ++i)
        {
            //Array<OneD, NekDouble> tFwd(nCoeffs);
            //m_fields[i]->FwdTrans_IterPerExp(m_temperature[i], tFwd);
            //fieldcoeffs.push_back(tFwd);
            
            variables.push_back("Temp" + boost::lexical_cast<std::string>(
                                    dimStr[i]));
        
        }
    }
}
