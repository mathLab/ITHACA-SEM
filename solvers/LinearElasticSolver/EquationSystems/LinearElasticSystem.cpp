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
#include <iomanip>

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

#ifdef NEKTAR_USE_MPI
#include <MultiRegions/GlobalLinSysXxtStaticCond.h>
#endif

#ifdef NEKTAR_USE_PETSC
#include <MultiRegions/GlobalLinSysPETScStaticCond.h>
#endif

using namespace std;

namespace Nektar
{

string LinearElasticSystem::className = GetEquationSystemFactory().
    RegisterCreatorFunction("LinearElasticSystem",
                            LinearElasticSystem::create);

/*
 * @brief Generate mapping between a simplex and the reference element.
 *
 * Mapping that requires the coordinates of the three vertices (x,y) of the
 * triangular element and outputs the function coefficients that defines the
 * mapping to the reference element for each coordinate (xfunc and yfunc)
 *
 * @param x   co-ordinates of triangle/tetrahedron
 * @param xf  components of mapping
 */
inline DNekMat MappingIdealToRef(SpatialDomains::GeometrySharedPtr geom)
{
    int n = geom->GetNumVerts(), i, j;

    DNekMat map   (n, n, 1.0, eFULL);
    DNekMat mapref(n, n, 1.0, eFULL);

    // Extract coordinate information.
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n-1; ++j)
        {
            map(j,i) = (*geom->GetVertex(i))[j];
        }
    }

    // Set up reference triangle or tetrahedron mapping.
    if (n == 3)
    {
        mapref(0,0) = -1.0; mapref(1,0) = -1.0;
        mapref(0,1) =  1.0; mapref(1,1) = -1.0;
        mapref(0,2) = -1.0; mapref(1,2) =  1.0;
    }
    else if (n == 4)
    {
        mapref(0,0) = -1.0; mapref(1,0) = -1.0; mapref(2,0) = -1.0;
        mapref(0,1) =  1.0; mapref(1,1) = -1.0; mapref(2,1) = -1.0;
        mapref(0,2) = -1.0; mapref(1,2) =  1.0; mapref(2,2) = -1.0;
        mapref(0,3) = -1.0; mapref(1,3) = -1.0; mapref(2,3) =  1.0;
    }

    map.Invert();

    DNekMat newmap = mapref * map;
    DNekMat mapred(n-1, n-1, 1.0, eFULL);

    for (i = 0; i < n-1; ++i)
    {
        for (j = 0; j < n-1; ++j)
        {
            mapred(i,j) = newmap(i,j);
        }
    }

    return mapred;
}


/**
 * @brief Default constructor.
 */
LinearElasticSystem::LinearElasticSystem(
    const LibUtilities::SessionReaderSharedPtr& pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph)
{
}

/**
 * @brief Set up the linear elasticity system.
 *
 * This routine loads the E and nu variables from the session file, creates a
 * coupled assembly map and creates containers for the statically condensed
 * block matrix system.
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
        MultiRegions::ContField2DSharedPtr u = std::dynamic_pointer_cast<
            MultiRegions::ContField2D>(m_fields[0]);
        m_assemblyMap = MemoryManager<CoupledAssemblyMap>
            ::AllocateSharedPtr(m_session,
                                m_graph,
                                u->GetLocalToGlobalMap(),
                                m_fields[0]->GetBndConditions(),
                                m_fields);
    }

    if (nVel == 3)
    {
        MultiRegions::ContField3DSharedPtr u = std::dynamic_pointer_cast<
            MultiRegions::ContField3D>(m_fields[0]);
        m_assemblyMap = MemoryManager<CoupledAssemblyMap>
            ::AllocateSharedPtr(m_session,
                                m_graph,
                                u->GetLocalToGlobalMap(),
                                m_fields[0]->GetBndConditions(),
                                m_fields);
    }

    // Figure out size of our new matrix systems by looping over all expansions
    // and multiply number of coefficients by velocity components.
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
        exp = m_fields[0]->GetExp(n);
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
 * formulation of the linear elasticity equations. The resulting matrix system
 * is then passed to LinearElasticSystem::SetStaticCondBlock in order to
 * construct the statically condensed system.
 *
 * All of the matrices involved in the construction of the divergence of the
 * stress tensor are Laplacian-like. We use the variable coefficient
 * functionality within the library to multiply the Laplacian by the appropriate
 * constants for all diagonal terms, and off-diagonal terms use
 * LinearElasticSystem::BuildLaplacianIJMatrix to construct matrices containing
 *
 * \f[ \int \partial_{x_i} \phi_k \partial_{x_j} \phi_l \f]
 *
 * where mixed derivative terms are present. Symmetry (in terms of k,l) is
 * exploited to avoid constructing this matrix repeatedly.
 *
 * @todo Make static condensation optional and construct full system instead.
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
    NekDouble mu     = m_E * 0.5  / (1.0 + m_nu);
    NekDouble lambda = m_E * m_nu / (1.0 + m_nu) / (1.0 - 2.0*m_nu);

    bool verbose = m_session->DefinesCmdLineArgument("verbose");
    bool root    = m_session->GetComm()->GetRank() == 0;

    // Loop over each element and construct matrices.
    if (nVel == 2)
    {
        for (n = 0; n < nEl; ++n)
        {
            exp = m_fields[0]->GetExp(n);
            const int nPhys = exp->GetTotPoints();
            Array<OneD, NekDouble> aArr(nPhys, lambda + 2*mu);
            Array<OneD, NekDouble> bArr(nPhys, mu);

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
            mat[0][1] = BuildLaplacianIJMatrix(1, 0, mu+lambda, exp);
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
            exp = m_fields[0]->GetExp(n);
            const int nPhys = exp->GetTotPoints();
            Array<OneD, NekDouble> aArr(nPhys, lambda + 2*mu);
            Array<OneD, NekDouble> bArr(nPhys, mu);

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
            mat[0][1] = BuildLaplacianIJMatrix(1, 0, mu + lambda, exp);
            mat[0][2] = BuildLaplacianIJMatrix(2, 0, mu + lambda, exp);

            mat[1][0] = mat[0][1];
            mat[1][1] = exp->GenMatrix(matkeyE);
            mat[1][2] = BuildLaplacianIJMatrix(2, 1, mu + lambda, exp);

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
 * - Scatter solution back to fields and backwards transform to physical space.
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
#ifdef NEKTAR_USE_PETSC
    else if (m_assemblyMap->GetGlobalSysSolnType() ==
             MultiRegions::ePETScStaticCond)
    {
        linSys = MemoryManager<
            MultiRegions::GlobalLinSysPETScStaticCond>::AllocateSharedPtr(
                key, m_fields[0], m_schurCompl, m_BinvD, m_C, m_Dinv,
                m_assemblyMap);
    }
#endif
#ifdef NEKTAR_USE_MPI
    else if (m_assemblyMap->GetGlobalSysSolnType() ==
             MultiRegions::eXxtStaticCond)
    {
        linSys = MemoryManager<
            MultiRegions::GlobalLinSysXxtStaticCond>::AllocateSharedPtr(
                key, m_fields[0], m_schurCompl, m_BinvD, m_C, m_Dinv,
                m_assemblyMap);
    }
#endif

    linSys->Initialise(m_assemblyMap);

    const int nCoeffs = m_fields[0]->GetNcoeffs();

    //
    // -- Evaluate forcing functions
    //

    // Evaluate the forcing function from the XML file.
    Array<OneD, Array<OneD, NekDouble> > forcing(nVel);
    GetFunction("Forcing")->Evaluate(forcing);

    // Add temperature term
    string tempEval;
    m_session->LoadSolverInfo("Temperature", tempEval, "None");

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
        ASSERTL0((m_fields[0]->GetCoordim(0)         == 2 &&
                  m_graph->GetAllQuadGeoms().size() == 0) ||
                 (m_fields[0]->GetCoordim(0)         == 3 &&
                  m_graph->GetAllPrismGeoms().size() == 0 &&
                  m_graph->GetAllPyrGeoms  ().size() == 0 &&
                  m_graph->GetAllHexGeoms  ().size() == 0),
                 "LinearIdealMetric temperature only implemented for "
                 "two-dimensional triangular meshes or three-dimensional "
                 "tetrahedral meshes.");

        m_temperature = Array<OneD, Array<OneD, NekDouble> >(nVel);
        m_stress = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(nVel);

        for (nv = 0; nv < nVel; ++nv)
        {
            m_temperature[nv] = Array<OneD, NekDouble>(
                m_fields[nv]->GetNpoints(), 0.0);
            m_stress[nv] = Array<OneD, Array<OneD, NekDouble> >(nVel);

            for (i = 0; i < nVel; ++i)
            {
                m_stress[nv][i] = Array<OneD, NekDouble>(
                    m_fields[nv]->GetNpoints(), 0.0);
            }
        }

        for (i = 0; i < m_fields[0]->GetExpSize(); ++i)
        {
            // Calculate element area
            LocalRegions::ExpansionSharedPtr exp =
                m_fields[0]->GetExp(i);
            LibUtilities::PointsKeyVector pkey = exp->GetPointsKeys();
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > deriv =
                exp->GetMetricInfo()->GetDeriv(pkey);
            int offset = m_fields[0]->GetPhys_Offset(i);

            DNekMat i2rm = MappingIdealToRef(exp->GetGeom());

            // Compute metric tensor
            DNekMat jac     (nVel, nVel, 0.0, eFULL);
            DNekMat jacIdeal(nVel, nVel, 0.0, eFULL);
            DNekMat metric  (nVel, nVel, 0.0, eFULL);

            for (j = 0; j < deriv[0][0].size(); ++j)
            {
                for (k = 0; k < nVel; ++k)
                {
                    for (l = 0; l < nVel; ++l)
                    {
                        jac(l,k) = deriv[k][l][j];
                    }
                }

                jacIdeal = jac * i2rm;
                metric = Transpose(jacIdeal) * jacIdeal;

                // Compute eigenvalues/eigenvectors of metric tensor using
                // ideal mapping.
                char jobvl = 'N', jobvr = 'V';
                int worklen = 8*nVel, info;

                DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
                DNekMat evec   (nVel, nVel, 0.0, eFULL);
                DNekMat evecinv(nVel, nVel, 0.0, eFULL);
                Array<OneD, NekDouble> vl  (nVel*nVel);
                Array<OneD, NekDouble> work(worklen);
                Array<OneD, NekDouble> wi  (nVel);

                Lapack::Dgeev(jobvl, jobvr, nVel, metric.GetRawPtr(), nVel,
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

                for (k = 0; k < nVel; ++k)
                {
                    for (l = 0; l < nVel; ++l)
                    {
                        m_stress[k][l][offset+j] = beta(k,l);
                    }
                }
            }

            if (deriv[0][0].size() != exp->GetTotPoints())
            {
                Array<OneD, NekDouble> tmp;
                for (k = 0; k < nVel; ++k)
                {
                    for (l = 0; l < nVel; ++l)
                    {
                        Vmath::Fill(
                            exp->GetTotPoints(), m_stress[k][l][offset],
                            tmp = m_stress[k][l] + offset, 1);
                    }
                }
            }
        }

        // Calculate divergence of stress tensor.
        Array<OneD, NekDouble> tmpderiv(m_fields[0]->GetNpoints());
        for (i = 0; i < nVel; ++i)
        {
            for (j = 0; j < nVel; ++j)
            {
                m_fields[i]->PhysDeriv(j, m_stress[i][j], tmpderiv);
                Vmath::Vadd (m_fields[i]->GetNpoints(), tmpderiv, 1,
                             m_temperature[i], 1, m_temperature[i], 1);
            }

            Vmath::Vcopy(m_fields[i]->GetNpoints(),
                         m_temperature[i], 1, forcing[i], 1);
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
    //   assembly map. We loop over each element, then the boundary degrees of
    //   freedom for u, boundary for v, followed by the interior for u and then
    //   interior for v.
    // - inout holds the Dirichlet degrees of freedom in the local ordering,
    //   which have the boundary conditions imposed 
    Array<OneD, NekDouble> forCoeffs(nVel * nCoeffs, 0.0);
    Array<OneD, NekDouble> inout    (nVel * nCoeffs, 0.0);
    Array<OneD, NekDouble> tmp      (nCoeffs);

    for (nv = 0; nv < nVel; ++nv)
    {
        // Take the inner product of the forcing function.
        m_fields[nv]->IProductWRTBase_IterPerExp(forcing[nv], tmp);

        // Impose Dirichlet condition on field which should be initialised 
        Array<OneD, NekDouble> loc_inout = m_fields[nv]->UpdateCoeffs(); 
        m_fields[nv]->ImposeDirichletConditions(loc_inout);
        
        // Scatter forcing into RHS vector according to the ordering dictated in
        // the comment above.
        for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
        {
            int nBnd   = m_bmap[i].size();
            int nInt   = m_imap[i].size();
            int offset = m_fields[nv]->GetCoeff_Offset(i);

            for (j = 0; j < nBnd; ++j)
            {
                forCoeffs[nVel*offset + nv*nBnd + j] =
                    -1.0*tmp[offset+m_bmap[i][j]];
                inout[nVel*offset + nv*nBnd + j] =
                    loc_inout[offset+m_bmap[i][j]];
            }
            for (j = 0; j < nInt; ++j)
            {
                forCoeffs[nVel*(offset + nBnd) + nv*nInt + j] =
                    -1.0*tmp[offset+m_imap[i][j]];
            }
        }
    }

    //
    // -- Perform solve
    //
    // Solve.
    linSys->Solve(forCoeffs, inout, m_assemblyMap);

    //
    // -- Postprocess
    //

    // Scatter back to field degrees of freedom
    for (nv = 0; nv < nVel; ++nv)
    {
        for (i = 0; i < m_fields[nv]->GetExpSize(); ++i)
        {
            int nBnd   = m_bmap[i].size();
            int nInt   = m_imap[i].size();
            int offset = m_fields[nv]->GetCoeff_Offset(i);

            for (j = 0; j < nBnd; ++j)
            {
                m_fields[nv]->UpdateCoeffs()[offset+m_bmap[i][j]] =
                    inout[nVel*offset + nv*nBnd + j];
            }
            for (j = 0; j < nInt; ++j)
            {
                m_fields[nv]->UpdateCoeffs()[offset+m_imap[i][j]] =
                    inout[nVel*(offset + nBnd) + nv*nInt + j];
            }
        }
        m_fields[nv]->BwdTrans(m_fields[nv]->GetCoeffs(),
                               m_fields[nv]->UpdatePhys());
    }
}

/**
 * @brief Given a block matrix for an element, construct its static condensation
 * matrices.
 *
 * This routine essentially duplicates the logic present in many of the
 * LocalRegions matrix generation routines to construct the statically condensed
 * equivalent of mat to pass to the GlobalLinSys solver.
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

    if (m_temperature.size() == 0)
    {
        return;
    }

    for (int i = 0; i < nVel; ++i)
    {
        Array<OneD, NekDouble> tFwd(nCoeffs);
        m_fields[i]->FwdTrans(m_temperature[i], tFwd);
        fieldcoeffs.push_back(tFwd);
        variables.push_back(
            "ThermStressDiv" + boost::lexical_cast<std::string>(i));
    }

    if (m_stress.size() == 0)
    {
        return;
    }

    for (int i = 0; i < nVel; ++i)
    {
        for (int j = 0; j < nVel; ++j)
        {
            Array<OneD, NekDouble> tFwd(nCoeffs);
            m_fields[i]->FwdTrans(m_stress[i][j], tFwd);
            fieldcoeffs.push_back(tFwd);
            variables.push_back(
                "ThermStress"
                + boost::lexical_cast<std::string>(i)
                + boost::lexical_cast<std::string>(j));
        }
    }
}

}
