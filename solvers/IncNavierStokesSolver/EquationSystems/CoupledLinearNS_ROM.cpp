///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS_ROM.cpp
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
// Description: Coupled  Solver for the Linearised Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>

#include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS_ROM.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <MultiRegions/ContField.h>


using namespace std;

namespace Nektar
{

    string CoupledLinearNS_ROM::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearisedNS_ROM", CoupledLinearNS_ROM::create);


    /**
     *  @class CoupledLinearNS_ROM
     *
     * Set up expansion field for velocity and pressure, the local to
     * global mapping arrays and the basic memory definitions for
     * coupled matrix system
     */
    CoupledLinearNS_ROM::CoupledLinearNS_ROM(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          IncNavierStokes(pSession, pGraph),
          m_zeroMode(false)
    {
    }

    void CoupledLinearNS_ROM::v_InitObject()
    {
        IncNavierStokes::v_InitObject();

        int  i;
        int  expdim = m_graph->GetMeshDimension();

        // Get Expansion list for orthogonal expansion at p-2
        const SpatialDomains::ExpansionInfoMap
            &pressure_exp = GenPressureExp(m_graph->GetExpansionInfo("u"));

        m_nConvectiveFields = m_fields.size();
        if(boost::iequals(m_boundaryConditions->GetVariable(m_nConvectiveFields-1), "p"))
        {
            ASSERTL0(false,"Last field is defined as pressure but this is not suitable for this solver, please remove this field as it is implicitly defined");
        }
        // Decide how to declare explist for pressure.
        if(expdim == 2)
        {
            int nz;

            if(m_HomogeneousType == eHomogeneous1D)
            {
                ASSERTL0(m_fields.size() > 2,"Expect to have three at least three components of velocity variables");
                LibUtilities::BasisKey Homo1DKey = m_fields[0]->GetHomogeneousBasis()->GetBasisKey();

                m_pressure = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(m_session, Homo1DKey, m_LhomZ, m_useFFT,m_homogen_dealiasing, pressure_exp);

                ASSERTL1(m_npointsZ%2==0,"Non binary number of planes have been specified");
                nz = m_npointsZ/2;

            }
            else
            {
                m_pressure = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr
                    (m_session, pressure_exp);
                nz = 1;
            }

            Array<OneD, MultiRegions::ExpListSharedPtr> velocity(m_velocity.size());
            for(i =0 ; i < m_velocity.size(); ++i)
            {
                velocity[i] = m_fields[m_velocity[i]];
            }

            // Set up Array of mappings
            m_locToGloMap = Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> (nz);

            if(m_singleMode)
            {

                ASSERTL0(nz <=2 ,"For single mode calculation can only have  nz <= 2");
                if(m_session->DefinesSolverInfo("BetaZero"))
                {
                    m_zeroMode = true;
                }
                int nz_loc = 2;
                m_locToGloMap[0] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure,nz_loc,m_zeroMode);
            }
            else
            {
                // base mode
                int nz_loc = 1;
                m_locToGloMap[0] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure,nz_loc);

                if(nz > 1)
                {
                    nz_loc = 2;
                    // Assume all higher modes have the same boundary conditions and re-use mapping
                    m_locToGloMap[1] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure->GetPlane(2),nz_loc,false);
                    // Note high order modes cannot be singular
                    for(i = 2; i < nz; ++i)
                    {
                        m_locToGloMap[i] = m_locToGloMap[1];
                    }
                }
            }
        }
        else if (expdim == 3)
        {
            //m_pressure = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(pressure_exp);
            ASSERTL0(false,"Setup mapping aray");
        }
        else
        {
            ASSERTL0(false,"Exp dimension not recognised");
        }

        // creation of the extrapolation object
        if(m_equationType == eUnsteadyNavierStokes)
        {
            std::string vExtrapolation = "Standard";

            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = m_session->GetSolverInfo("Extrapolation");
            }

            m_extrapolation = GetExtrapolateFactory().CreateInstance(
                vExtrapolation,
                m_session,
                m_fields,
                m_pressure,
                m_velocity,
                m_advObject);
        }

    }

    /**
     * Set up a coupled linearised Naviers-Stokes solve. Primarily a
     * wrapper fuction around the full mode by mode version of
     * SetUpCoupledMatrix.
     *
     */
    void CoupledLinearNS_ROM::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation)
    {

        int nz;
        if(m_singleMode)
        {

            NekDouble lambda_imag;

            // load imaginary component of any potential shift
            // Probably should be called from DriverArpack but not yet
            // clear how to do this
            m_session->LoadParameter("imagShift",lambda_imag,NekConstants::kNekUnsetDouble);
            nz = 1;
            m_mat  = Array<OneD, CoupledSolverMatrices> (nz);

            ASSERTL1(m_npointsZ <=2,"Only expected a maxmimum of two planes in single mode linear NS solver");

            if(m_zeroMode)
            {
                SetUpCoupledMatrix(lambda,Advfield,IsLinearNSEquation,0,m_mat[0],m_locToGloMap[0],lambda_imag);
            }
            else
            {
                NekDouble beta =  2*M_PI/m_LhomZ;
                NekDouble lam = lambda + m_kinvis*beta*beta;

                SetUpCoupledMatrix(lam,Advfield,IsLinearNSEquation,1,m_mat[0],m_locToGloMap[0],lambda_imag);
            }
        }
        else
        {
            int n;
            if(m_npointsZ > 1)
            {
                nz = m_npointsZ/2;
            }
            else
            {
                nz =  1;
            }

            m_mat  = Array<OneD, CoupledSolverMatrices> (nz);

            // mean mode or 2D mode.
            SetUpCoupledMatrix(lambda,Advfield,IsLinearNSEquation,0,m_mat[0],m_locToGloMap[0]);

            for(n = 1; n < nz; ++n)
            {
                NekDouble beta = 2*M_PI*n/m_LhomZ;

                NekDouble lam = lambda + m_kinvis*beta*beta;

                SetUpCoupledMatrix(lam,Advfield,IsLinearNSEquation,n,m_mat[n],m_locToGloMap[n]);
            }
        }

    }


    /**
     * Set up a coupled linearised Naviers-Stokes solve in the
     * following manner:
     *
     * Consider the linearised Navier-Stokes element matrix denoted as
     *
     * \f[ \left [ \begin{array}{ccc} A
     * & -D_{bnd}^T & B\\ -D_{bnd}& 0 & -D_{int}\\ \tilde{B}^T &
     * -D_{int}^T & C \end{array} \right ] \left [ \begin{array}{c} {\bf
     * v}_{bnd} \\ p\\ {\bf v}_{int} \end{array} \right ] = \left [
     * \begin{array}{c} {\bf f}_{bnd} \\ 0\\ {\bf f}_{int} \end{array}
     * \right ] \f]
     *
     * where \f${\bf v}_{bnd}, {\bf f}_{bnd}\f$ denote the degrees of
     * freedom of the elemental velocities on the boundary of the
     * element, \f${\bf v}_{int}, {\bf f}_{int}\f$ denote the degrees
     * of freedom of the elemental velocities on the interior of the
     * element and \f$p\f$ is the piecewise continuous pressure.  The
     * matrices have the interpretation
     *
     * \f[  A[n,m]  = (\nabla \phi^b_n, \nu \nabla
     * \phi^b_m) + (\phi^b_n,{\bf U \cdot \nabla} \phi^b_m) +
     * (\phi^b_n \nabla^T {\bf U} \phi^b_m), \f]
     *
     * \f[ B[n,m] = (\nabla \phi^b_n, \nu \nabla \phi^i_m) +
     * (\phi^b_n,{\bf U \cdot \nabla} \phi^i_m) + (\phi^b_n \nabla^T
     * {\bf U} \phi^i_m),\f]
     *
     * \f[\tilde{B}^T[n,m] = (\nabla \phi^i_n, \nu \nabla \phi^b_m) +
     * (\phi^i_n,{\bf U \cdot \nabla} \phi^b_m) + (\phi^i_n \nabla^T
     * {\bf U} \phi^b_m) \f]
     *
     * \f[ C[n,m] = (\nabla \phi^i_n, \nu \nabla
     * \phi^i_m) + (\phi^i_n,{\bf U \cdot \nabla} \phi^i_m) +
     * (\phi^i_n \nabla^T {\bf U} \phi^i_m),\f]
     *
     * \f[ D_{bnd}[n,m] = (\psi_m,\nabla \phi^b),\f]
     *
     * \f[ D_{int}[n,m] = (\psi_m,\nabla \phi^i) \f]
     *
     * where \f$\psi\f$ is the space of pressures typically at order
     * \f$P-2\f$ and \f$\phi\f$ is the velocity vector space of
     * polynomial order \f$P\f$. (Note the above definitions for the
     * transpose terms shoudl be interpreted with care since we have
     * used a tensor differential for the \f$\nabla^T \f$ terms)
     *
     * Note \f$B = \tilde{B}^T\f$ if just considering the stokes
     * operator and then \f$C\f$ is also symmetric. Also note that
     * \f$A,B\f$ and \f$C\f$ are block diagonal in the Oseen equation
     * when \f$\nabla^T {\bf U}\f$ are zero.
     *
     * Since \f$C\f$ is invertible we can premultiply the governing
     * elemental equation by the following matrix:
     *
     * \f[  \left [ \begin{array}{ccc} I & 0 &
     * -BC^{-1}\\ 0& I & D_{int}C^{-1}\\ 0 & 0 & I \end{array}
     * \right ] \left \{ \left [ \begin{array}{ccc} A & -D_{bnd}^T &
     * B\\ -D_{bnd}& 0 & -D_{int}\\ \tilde{B}^T & -D_{int}^T & C
     * \end{array} \right ] \left [ \begin{array}{c} {\bf v}_{bnd} \\
     * p\\ {\bf v_{int}} \end{array} \right ] = \left [
     * \begin{array}{c} {\bf f}_{bnd} \\ 0\\ {\bf f_{int}} \end{array}
     * \right ] \right \} \f]
     *
     *  which if we multiply out the matrix equation we obtain
     *
     * \f[ \left [ \begin{array}{ccc} A - B
     * C^{-1}\tilde{B}^T & -D_{bnd}^T +B C^{-1} D_{int}^T& 0\\
     * -D_{bnd}+D_{int} C^{-1} \tilde{B}^T & -D_{int} C^{-1}
     * D_{int}^T & 0\\ \tilde{B}^T & -D_{int}^T & C \end{array} \right ]
     * \left [ \begin{array}{c} {\bf v}_{bnd} \\ p\\ {\bf v_{int}}
     * \end{array} \right ] = \left [ \begin{array}{c} {\bf f}_{bnd}
     * -B C^{-1} {\bf f}_{int}\\ f_p = D_{int} C^{-1} {\bf
     * f}_{int}\\ {\bf f_{int}} \end{array} \right ] \f]
     *
     * In the above equation the \f${\bf v}_{int}\f$ degrees of freedom
     * are decoupled and so we need to solve for the \f${\bf v}_{bnd},
     * p\f$ degrees of freedom.  The final step is to perform a second
     * level of static condensation but where we will lump the mean
     * pressure mode (or a pressure degree of freedom containing a
     * mean component) with the velocity boundary degrees of
     * freedom. To do we define \f${\bf b} = [{\bf v}_{bnd}, p_0]\f$ where
     * \f$p_0\f$ is the mean pressure mode and \f$\hat{p}\f$ to be the
     * remainder of the pressure space.  We now repartition the top \f$2
     * \times 2\f$ block of matrices of previous matrix equation  as
     *
     * \f[ \left [ \begin{array}{cc} \hat{A} & \hat{B}\\ \hat{C} &
     * \hat{D} \end{array} \right ] \left [ \begin{array}{c} {\bf b}
     * \\ \hat{p} \end{array} \right ] = \left [ \begin{array}{c}
     * \hat{\bf f}_{bnd} \\ \hat{f}_p \end{array} \right ]
     * \label{eqn.linNS_fac2} \f]
     *
     * where
     *
     * \f[ \hat{A}_{ij} = \left [ \begin{array}{cc} A - B
     * C^{-1}\tilde{B}^T & [-D_{bnd}^T +B C^{-1} D_{int}^T]_{i0}\\
     * {[}-D_{bnd}+D_{int} C^{-1} \tilde{B}^T]_{0j} & -[D_{int}
     * C^{-1} D_{int}^T ]_{00} \end{array} \right ] \f]
     *
     * \f[ \hat{B}_{ij} = \left [ \begin{array}{c} [-D_{bnd}^T +B
     * C^{-1} D_{int}^T]_{i,j+1} \\ {[} -D_{int} C^{-1} D^T_{int}
     * ]_{0j}\end{array} \right ] \f]
     *
     * \f[ \hat{C}_{ij}  =  \left [\begin{array}{cc} -D_{bnd} +
     * D_{int} C^{-1} \tilde{B}^T, & {[} -D_{int} C^{-1} D^T_{int}
     * ]_{i+1,0}\end{array} \right ] \f]
     *
     * \f[ \hat{D}_{ij}  =  \left [\begin{array}{c} {[} -D_{int}
     * C^{-1} D^T_{int} ]_{i+1,j+1}\end{array} \right ] \f]
     *
     * and
     *
     * \f[ fh\_{bnd} = \left [ \begin{array}{c} {\bf
     * f}_{bnd} -B C^{-1} {\bf f}_{int}\\ {[}D_{int} C^{-1} {\bf
     * f}_{int}]_0 \end{array}\right ] \hspace{1cm} [fh\_p_{i} =
     * \left [ \begin{array}{c} {[}D_{int} C^{-1} {\bf f}_{iint}]_{i+1}
     * \end{array}\right ] \f]
     *
     * Since the \f$\hat{D}\f$ is decoupled and invertible we can now
     * statically condense the previous matrix equationto decouple
     * \f${\bf b}\f$ from \f$\hat{p}\f$ by solving the following system
     *
     * \f[ \left [ \begin{array}{cc} \hat{A} - \hat{B} \hat{D}^{-1}
     * \hat{C} & 0 \\ \hat{C} & \hat{D} \end{array} \right ] \left
     * [ \begin{array}{c} {\bf b} \\ \hat{p} \end{array} \right ] =
     * \left [ \begin{array}{c} \hat{\bf f}_{bnd} - \hat{B}
     * \hat{D}^{-1} \hat{f}_p\\ \hat{f}_p \end{array} \right ] \f]
     *
     * The matrix \f$\hat{A} - \hat{B} \hat{D}^{-1} \hat{C}\f$ has
     * to be globally assembled and solved iteratively or
     * directly. One we obtain the solution to \f${\bf b}\f$ we can use
     * the second row of equation fourth matrix equation to solve for
     * \f$\hat{p}\f$ and finally the last row of equation second
     * matrix equation to solve for \f${\bf v}_{int}\f$.
     *
     */

    void CoupledLinearNS_ROM::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation,const int HomogeneousMode, CoupledSolverMatrices &mat, CoupledLocalToGlobalC0ContMapSharedPtr &locToGloMap, const NekDouble lambda_imag)
    {
        int  n,i,j,k;
        int  nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int  nvel   = m_velocity.size();

        // if Advfield is defined can assume it is an Oseen or LinearNS equation
        bool AddAdvectionTerms = (Advfield ==  NullNekDoubleArrayOfArray)? false: true;
        //bool AddAdvectionTerms = true; // Temporary debugging trip

        // call velocity space Static condensation and fill block
        // matrices.  Need to set this up differently for Oseen and
        // Lin NS.  Ideally should make block diagonal for Stokes and
        // Oseen problems.
        DNekScalMatSharedPtr loc_mat;
        StdRegions::StdExpansionSharedPtr locExp;
        NekDouble one = 1.0;
        int nint,nbndry;
        int rows, cols;
        NekDouble zero = 0.0;
        Array<OneD, unsigned int> bmap,imap;

        Array<OneD,unsigned int> nsize_bndry   (nel);
        Array<OneD,unsigned int> nsize_bndry_p1(nel);
        Array<OneD,unsigned int> nsize_int     (nel);
        Array<OneD,unsigned int> nsize_p       (nel);
        Array<OneD,unsigned int> nsize_p_m1    (nel);

        int nz_loc;

        if(HomogeneousMode) // Homogeneous mode flag
        {
            nz_loc = 2;
        }
        else
        {
            if(m_singleMode)
            {
                nz_loc = 2;
            }
            else
            {
                nz_loc = 1;
            }
        }

        // Set up block matrix sizes -
        for(n = 0; n < nel; ++n)
        {
            nsize_bndry[n] = nvel*m_fields[m_velocity[0]]->GetExp(n)->NumBndryCoeffs()*nz_loc;
            nsize_bndry_p1[n] = nsize_bndry[n]+nz_loc;
            nsize_int  [n] = (nvel*m_fields[m_velocity[0]]->GetExp(n)->GetNcoeffs()*nz_loc - nsize_bndry[n]);
            nsize_p[n] = m_pressure->GetExp(n)->GetNcoeffs()*nz_loc;
            nsize_p_m1[n] = nsize_p[n]-nz_loc;
        }

        MatrixStorage blkmatStorage = eDIAGONAL;
        DNekScalBlkMatSharedPtr pAh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry_p1,nsize_bndry_p1,blkmatStorage);
        mat.m_BCinv = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
        mat.m_Btilde = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
        mat.m_Cinv = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_int,nsize_int,blkmatStorage);

        mat.m_D_bnd = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p,nsize_bndry,blkmatStorage);

        mat.m_D_int = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p,nsize_int,blkmatStorage);

        // Final level static condensation matrices.
        DNekScalBlkMatSharedPtr pBh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry_p1,nsize_p_m1,blkmatStorage);
        DNekScalBlkMatSharedPtr pCh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p_m1,nsize_bndry_p1,blkmatStorage);
        DNekScalBlkMatSharedPtr pDh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p_m1,nsize_p_m1,blkmatStorage);

        LibUtilities::Timer timer;
        timer.Start();

        for(n = 0; n < nel; ++n)
        {
            nbndry = nsize_bndry[n];
            nint   = nsize_int[n];
            k = nsize_bndry_p1[n];
            DNekMatSharedPtr Ah = MemoryManager<DNekMat>::AllocateSharedPtr(k,k,zero);
            Array<OneD, NekDouble> Ah_data = Ah->GetPtr();
            int AhRows = k;
            DNekMatSharedPtr B  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            Array<OneD, NekDouble> B_data = B->GetPtr();
            DNekMatSharedPtr C  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            Array<OneD, NekDouble> C_data = C->GetPtr();
            DNekMatSharedPtr D  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);
            Array<OneD, NekDouble> D_data = D->GetPtr();

            DNekMatSharedPtr Dbnd = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p[n],nsize_bndry[n],zero);
            DNekMatSharedPtr Dint = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p[n],nsize_int[n],zero);

            locExp = m_fields[m_velocity[0]]->GetExp(n);
            locExp->GetBoundaryMap(bmap);
            locExp->GetInteriorMap(imap);
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorLambda] = lambda/m_kinvis;
            LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


            int ncoeffs = m_fields[m_velocity[0]]->GetExp(n)->GetNcoeffs();
            int nphys   = m_fields[m_velocity[0]]->GetExp(n)->GetTotPoints();
            int nbmap = bmap.size();
            int nimap = imap.size();

            Array<OneD, NekDouble> coeffs(ncoeffs);
            Array<OneD, NekDouble> phys  (nphys);
            int psize   = m_pressure->GetExp(n)->GetNcoeffs();
            int pqsize  = m_pressure->GetExp(n)->GetTotPoints();

            Array<OneD, NekDouble> deriv  (pqsize);
            Array<OneD, NekDouble> pcoeffs(psize);

            if(AddAdvectionTerms == false) // use static condensed managed matrices
            {
                // construct velocity matrices using statically
                // condensed elemental matrices and then construct
                // pressure matrix systems
                DNekScalBlkMatSharedPtr CondMat;
                CondMat = locExp->GetLocStaticCondMatrix(helmkey);

                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    DNekScalMat &SubBlock = *CondMat->GetBlock(0,0);
                    rows = SubBlock.GetRows();
                    cols = SubBlock.GetColumns();
                    for(i = 0; i < rows; ++i)
                    {
                        for(j = 0; j < cols; ++j)
                        {
                            (*Ah)(i+k*rows,j+k*cols) = m_kinvis*SubBlock(i,j);
                        }
                    }
                }

                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    DNekScalMat &SubBlock  = *CondMat->GetBlock(0,1);
                    DNekScalMat &SubBlock1 = *CondMat->GetBlock(1,0);
                    rows = SubBlock.GetRows();
                    cols = SubBlock.GetColumns();
                    for(i = 0; i < rows; ++i)
                    {
                        for(j = 0; j < cols; ++j)
                        {
                            (*B)(i+k*rows,j+k*cols) = SubBlock(i,j);
                            (*C)(i+k*rows,j+k*cols) = m_kinvis*SubBlock1(j,i);
                        }
                    }
                }

                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    DNekScalMat &SubBlock = *CondMat->GetBlock(1,1);
                    NekDouble inv_kinvis = 1.0/m_kinvis;
                    rows = SubBlock.GetRows();
                    cols = SubBlock.GetColumns();
                    for(i = 0; i < rows; ++i)
                    {
                        for(j = 0; j < cols; ++j)
                        {
                            (*D)(i+k*rows,j+k*cols) = inv_kinvis*SubBlock(i,j);
                        }
                    }
                }


                // Loop over pressure space and construct boundary block matrices.
                for(i = 0; i < bmap.size(); ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[bmap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(n)->BwdTrans(coeffs,phys);

                    // Differentiation & Inner product wrt base.
                    for(j = 0; j < nvel; ++j)
                    {
                        if( (nz_loc == 2)&&(j == 2)) // handle d/dz derivative
                        {
                            NekDouble beta =  -2*M_PI*HomogeneousMode/m_LhomZ;

                            Vmath::Smul(m_fields[m_velocity[0]]->GetExp(n)->GetTotPoints(), beta, phys,1,deriv,1);

                            m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);

                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() +
                                        ((nz_loc*j+1)*bmap.size()+i)*nsize_p[n],1);

                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() +
                                        ((nz_loc*j)*bmap.size()+i)*nsize_p[n]+psize,1);

                        }
                        else
                        {
                            if(j < 2) // required for mean mode of homogeneous expansion
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[j],phys,deriv);
                                m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);
                                // copy into column major storage.
                                for(k = 0; k < nz_loc; ++k)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dbnd->GetRawPtr() +
                                                ((nz_loc*j+k)*bmap.size()+i)*nsize_p[n]+ k*psize,1);
                                }
                            }
                        }
                    }
                }

                for(i = 0; i < imap.size(); ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[imap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(n)->BwdTrans(coeffs,phys);

                    // Differentiation & Inner product wrt base.
                    for(j = 0; j < nvel; ++j)
                    {
                        if( (nz_loc == 2)&&(j == 2)) // handle d/dz derivative
                        {
                            NekDouble beta = -2*M_PI*HomogeneousMode/m_LhomZ;

                            Vmath::Smul(m_fields[m_velocity[0]]->GetExp(n)->GetTotPoints(), beta, phys,1,deriv,1);

                            m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);

                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() +
                                        ((nz_loc*j+1)*imap.size()+i)*nsize_p[n],1);

                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() +
                                        ((nz_loc*j)*imap.size()+i)*nsize_p[n]+psize,1);

                        }
                        else
                        {
                            if(j < 2) // required for mean mode of homogeneous expansion
                            {
                                //m_fields[m_velocity[0]]->GetExp(n)->PhysDeriv(j,phys, deriv);
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[j],phys,deriv);

                                m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);

                                // copy into column major storage.
                                for(k = 0; k < nz_loc; ++k)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dint->GetRawPtr() +
                                                ((nz_loc*j+k)*imap.size()+i)*nsize_p[n]+ k*psize,1);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // construct velocity matrices and pressure systems at
                // the same time resusing differential of velocity
                // space

                DNekScalMat &HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                         ->GetLocMatrix(helmkey));
                DNekScalMatSharedPtr MassMat;

                Array<OneD, const NekDouble> HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
                NekDouble HelmMatScale = HelmMat.Scale();
                int HelmMatRows = HelmMat.GetRows();

                if((lambda_imag != NekConstants::kNekUnsetDouble)&&(nz_loc == 2))
                {
                    LocalRegions::MatrixKey masskey(StdRegions::eMass,
                                                    locExp->DetShapeType(),
                                                    *locExp);
                    MassMat = locExp->as<LocalRegions::Expansion>()
                        ->GetLocMatrix(masskey);
                }

                Array<OneD, NekDouble> Advtmp;
                Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
                // Use ExpList phys array for temporaary storage
                Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
                int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(n);
                int nv;
                int npoints = locExp->GetTotPoints();

                // Calculate derivative of base flow
                if(IsLinearNSEquation)
                {
                    int nv1;
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(nv = 0; nv < nvel; ++nv)
                    {
                        for(nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }

                            if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
                            {
                                Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
                            }
                            else
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                            }
                        }
                    }
                }


                for(i = 0; i < nbmap; ++i)
                {

                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[bmap[i]] = 1.0;
                    locExp->BwdTrans(coeffs,phys);

                    for(k = 0; k < nvel*nz_loc; ++k)
                    {
                        for(j = 0; j < nbmap; ++j)
                        {
                            //                            Ah_data[i+k*nbmap + (j+k*nbmap)*AhRows] += m_kinvis*HelmMat(bmap[i],bmap[j]);
                            Ah_data[i+k*nbmap + (j+k*nbmap)*AhRows] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i] + HelmMatRows*bmap[j]];
                        }

                        for(j = 0; j < nimap; ++j)
                        {
                            B_data[i+k*nbmap + (j+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i]+HelmMatRows*imap[j]];
                        }
                    }

                    if((lambda_imag != NekConstants::kNekUnsetDouble)&&(nz_loc == 2))
                    {
                        for(k = 0; k < nvel; ++k)
                        {
                            for(j = 0; j < nbmap; ++j)
                            {
                                Ah_data[i+2*k*nbmap + (j+(2*k+1)*nbmap)*AhRows] -= lambda_imag*(*MassMat)(bmap[i],bmap[j]);
                            }

                            for(j = 0; j < nbmap; ++j)
                            {
                                Ah_data[i+(2*k+1)*nbmap + (j+2*k*nbmap)*AhRows] += lambda_imag*(*MassMat)(bmap[i],bmap[j]);
                            }

                            for(j = 0; j < nimap; ++j)
                            {
                                B_data[i+2*k*nbmap + (j+(2*k+1)*nimap)*nbndry] -= lambda_imag*(*MassMat)(bmap[i],imap[j]);
                            }

                            for(j = 0; j < nimap; ++j)
                            {
                                B_data[i+(2*k+1)*nbmap + (j+2*k*nimap)*nbndry] += lambda_imag*(*MassMat)(bmap[i],imap[j]);
                            }

                        }
                    }



                    for(k = 0; k < nvel; ++k)
                    {
                        if((nz_loc == 2)&&(k == 2)) // handle d/dz derivative
                        {
                            NekDouble beta = -2*M_PI*HomogeneousMode/m_LhomZ;

                            // Real Component
                            Vmath::Smul(npoints,beta,phys,1,deriv,1);

                            m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() +
                                        ((nz_loc*k+1)*bmap.size()+i)*nsize_p[n],1);

                            // Imaginary Component
                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() +
                                        ((nz_loc*k)*bmap.size()+i)*nsize_p[n]+psize,1);

                            // now do advection terms
                            Vmath::Vmul(npoints,
                                        Advtmp = Advfield[k] + phys_offset,
                                        1,deriv,1,tmpphys,1);

                            locExp->IProductWRTBase(tmpphys,coeffs);


                            // real contribution
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    Ah_data[j+2*nv*nbmap + (i+(2*nv+1)*nbmap)*AhRows] +=
                                        coeffs[bmap[j]];
                                }

                                for(j = 0; j < nimap; ++j)
                                {
                                    C_data[i+(2*nv+1)*nbmap + (j+2*nv*nimap)*nbndry] += 
                                        coeffs[imap[j]];
                                }
                            }

                            Vmath::Neg(ncoeffs,coeffs,1);
                            // imaginary contribution
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    Ah_data[j+(2*nv+1)*nbmap + (i+2*nv*nbmap)*AhRows] +=
                                        coeffs[bmap[j]];
                                }

                                for(j = 0; j < nimap; ++j)
                                {
                                    C_data[i+2*nv*nbmap + (j+(2*nv+1)*nimap)*nbndry] += 
                                        coeffs[imap[j]];
                                }
                            }
                        }
                        else
                        {
                            if(k < 2)
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[k],phys, deriv);
                                Vmath::Vmul(npoints,
                                            Advtmp = Advfield[k] + phys_offset,
                                            1,deriv,1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);


                                for(nv = 0; nv < nvel*nz_loc; ++nv)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        Ah_data[j+nv*nbmap + (i+nv*nbmap)*AhRows] +=
                                            coeffs[bmap[j]];
                                    }

                                    for(j = 0; j < nimap; ++j)
                                    {
                                        C_data[i+nv*nbmap + (j+nv*nimap)*nbndry] += 
                                            coeffs[imap[j]];
                                    }
                                }

                                // copy into column major storage.
                                m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);
                                for(j = 0; j < nz_loc; ++j)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dbnd->GetRawPtr() +
                                                ((nz_loc*k+j)*bmap.size() + i)*nsize_p[n]+ j*psize,1);
                                }
                            }
                        }

                        if(IsLinearNSEquation)
                        {
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms
                                Vmath::Vmul(npoints,phys,1,
                                            AdvDeriv[k*nvel+nv],
                                            1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);

                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        Ah_data[j+(k*nz_loc+n1)*nbmap +
                                        (i+(nv*nz_loc+n1)*nbmap)*AhRows] +=
                                        coeffs[bmap[j]];
                                    }

                                    for(j = 0; j < nimap; ++j)
                                    {
                                        C_data[i+(nv*nz_loc+n1)*nbmap +
                                        (j+(k*nz_loc+n1)*nimap)*nbndry] +=
                                        coeffs[imap[j]];
                                    }
                                }
                            }
                        }
                    }
                }


                for(i = 0; i < nimap; ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[imap[i]] = 1.0;
                    locExp->BwdTrans(coeffs,phys);

                    for(k = 0; k < nvel*nz_loc; ++k)
                    {
                        for(j = 0; j < nbmap; ++j) // C set up as transpose
                        {
                            C_data[j+k*nbmap + (i+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*bmap[j]];
                        }

                        for(j = 0; j < nimap; ++j)
                        {
                            D_data[i+k*nimap + (j+k*nimap)*nint] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*imap[j]];
                        }
                    }

                    if((lambda_imag != NekConstants::kNekUnsetDouble)&&(nz_loc == 2))
                    {
                        for(k = 0; k < nvel; ++k)
                        {
                            for(j = 0; j < nbmap; ++j) // C set up as transpose
                            {
                                C_data[j+2*k*nbmap + (i+(2*k+1)*nimap)*nbndry] += lambda_imag*(*MassMat)(bmap[j],imap[i]);
                            }

                            for(j = 0; j < nbmap; ++j) // C set up as transpose
                            {
                                C_data[j+(2*k+1)*nbmap + (i+2*k*nimap)*nbndry] -= lambda_imag*(*MassMat)(bmap[j],imap[i]);
                            }

                            for(j = 0; j < nimap; ++j)
                            {
                                D_data[i+2*k*nimap + (j+(2*k+1)*nimap)*nint] -= lambda_imag*(*MassMat)(imap[i],imap[j]);
                            }

                            for(j = 0; j < nimap; ++j)
                            {
                                D_data[i+(2*k+1)*nimap + (j+2*k*nimap)*nint] += lambda_imag*(*MassMat)(imap[i],imap[j]);
                            }
                        }
                    }


                    for(k = 0; k < nvel; ++k)
                    {
                        if((nz_loc == 2)&&(k == 2)) // handle d/dz derivative
                        {
                            NekDouble beta = -2*M_PI*HomogeneousMode/m_LhomZ;

                            // Real Component
                            Vmath::Smul(npoints,beta,phys,1,deriv,1);
                            m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() +
                                        ((nz_loc*k+1)*imap.size()+i)*nsize_p[n],1);
                            // Imaginary Component
                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() +
                                        ((nz_loc*k)*imap.size()+i)*nsize_p[n]+psize,1);

                            // Advfield[k] *d/dx_k to all velocity
                            // components on diagonal
                            Vmath::Vmul(npoints, Advtmp = Advfield[k] + phys_offset,1,deriv,1,tmpphys,1);
                            locExp->IProductWRTBase(tmpphys,coeffs);

                            // Real Components
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    B_data[j+2*nv*nbmap + (i+(2*nv+1)*nimap)*nbndry] +=
                                    coeffs[bmap[j]];
                                }

                                for(j = 0; j < nimap; ++j)
                                {
                                    D_data[j+2*nv*nimap + (i+(2*nv+1)*nimap)*nint] +=
                                    coeffs[imap[j]];
                                }
                            }
                            Vmath::Neg(ncoeffs,coeffs,1);
                            // Imaginary
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    B_data[j+(2*nv+1)*nbmap + (i+2*nv*nimap)*nbndry] +=
                                    coeffs[bmap[j]];
                                }

                                for(j = 0; j < nimap; ++j)
                                {
                                    D_data[j+(2*nv+1)*nimap + (i+2*nv*nimap)*nint] +=
                                    coeffs[imap[j]];
                                }
                            }

                        }
                        else
                        {
                            if(k < 2)
                            {
                                // Differentiation & Inner product wrt base.
                                //locExp->PhysDeriv(k,phys, deriv);
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[k],phys,deriv);
                                Vmath::Vmul(npoints,
                                            Advtmp = Advfield[k] + phys_offset,
                                            1,deriv,1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);


                                for(nv = 0; nv < nvel*nz_loc; ++nv)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        B_data[j+nv*nbmap + (i+nv*nimap)*nbndry] +=
                                            coeffs[bmap[j]];
                                    }

                                    for(j = 0; j < nimap; ++j)
                                    {
                                        D_data[j+nv*nimap + (i+nv*nimap)*nint] +=
                                        coeffs[imap[j]];
                                    }
                                }
                                // copy into column major storage.
                                m_pressure->GetExp(n)->IProductWRTBase(deriv,pcoeffs);
                                for(j = 0; j < nz_loc; ++j)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dint->GetRawPtr() +
                                                ((nz_loc*k+j)*imap.size() + i)*nsize_p[n]+j*psize,1);
                                }
                            }
                        }

                        if(IsLinearNSEquation)
                        {
                            int n1;
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                // u'.Grad U terms
                                Vmath::Vmul(npoints,phys,1,
                                            AdvDeriv[k*nvel+nv],
                                            1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);

                                for(n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        B_data[j+(k*nz_loc+n1)*nbmap +
                                        (i+(nv*nz_loc+n1)*nimap)*nbndry] +=
                                        coeffs[bmap[j]];
                                    }

                                    for(j = 0; j < nimap; ++j)
                                    {
                                        D_data[j+(k*nz_loc+n1)*nimap +
                                        (i+(nv*nz_loc+n1)*nimap)*nint] +=
                                        coeffs[imap[j]];
                                    }
                                }
                            }
                        }
                    }
                }


                D->Invert();
                (*B) = (*B)*(*D);


                // perform (*Ah) = (*Ah) - (*B)*(*C) but since size of
                // Ah is larger than (*B)*(*C) easier to call blas
                // directly
                Blas::Dgemm('N','T', B->GetRows(), C->GetRows(),
                            B->GetColumns(), -1.0, B->GetRawPtr(),
                            B->GetRows(), C->GetRawPtr(),
                            C->GetRows(), 1.0,
                            Ah->GetRawPtr(), Ah->GetRows());
            }

            mat.m_BCinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
            mat.m_Btilde->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,C));
            mat.m_Cinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,D));
            mat.m_D_bnd->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dbnd));
            mat.m_D_int->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dint));

            // Do matrix manipulations and get final set of block matries
            // reset boundary to put mean mode into boundary system.

            DNekMatSharedPtr Cinv,BCinv,Btilde;
            DNekMat  DintCinvDTint, BCinvDTint_m_DTbnd, DintCinvBTtilde_m_Dbnd;

            Cinv   = D;
            BCinv  = B;
            Btilde = C;

            DintCinvDTint      = (*Dint)*(*Cinv)*Transpose(*Dint);
            BCinvDTint_m_DTbnd = (*BCinv)*Transpose(*Dint) - Transpose(*Dbnd);

            // This could be transpose of BCinvDint in some cases
            DintCinvBTtilde_m_Dbnd = (*Dint)*(*Cinv)*Transpose(*Btilde) - (*Dbnd);

            // Set up final set of matrices.
            DNekMatSharedPtr Bh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_bndry_p1[n],nsize_p_m1[n],zero);
            DNekMatSharedPtr Ch = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n],nsize_bndry_p1[n],zero);
            DNekMatSharedPtr Dh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n], nsize_p_m1[n],zero);
            Array<OneD, NekDouble> Dh_data = Dh->GetPtr();

            // Copy matrices into final structures.
            int n1,n2;
            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(i = 0; i < psize-1; ++i)
                {
                    for(n2 = 0; n2 < nz_loc; ++n2)
                    {
                        for(j = 0; j < psize-1; ++j)
                        {
                            //(*Dh)(i+n1*(psize-1),j+n2*(psize-1)) =
                            //-DintCinvDTint(i+1+n1*psize,j+1+n2*psize);
                            Dh_data[(i+n1*(psize-1)) + (j+n2*(psize-1))*Dh->GetRows()] =
                            -DintCinvDTint(i+1+n1*psize,j+1+n2*psize);
                        }
                    }
                }
            }

            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(i = 0; i < nsize_bndry_p1[n]-nz_loc; ++i)
                {
                    (*Ah)(i,nsize_bndry_p1[n]-nz_loc+n1) = BCinvDTint_m_DTbnd(i,n1*psize);
                    (*Ah)(nsize_bndry_p1[n]-nz_loc+n1,i) = DintCinvBTtilde_m_Dbnd(n1*psize,i);
                }
            }

            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(n2 = 0; n2 < nz_loc; ++n2)
                {
                    (*Ah)(nsize_bndry_p1[n]-nz_loc+n1,nsize_bndry_p1[n]-nz_loc+n2) =
                    -DintCinvDTint(n1*psize,n2*psize);
                }
            }

            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(j = 0; j < psize-1; ++j)
                {
                    for(i = 0; i < nsize_bndry_p1[n]-nz_loc; ++i)
                    {
                        (*Bh)(i,j+n1*(psize-1)) = BCinvDTint_m_DTbnd(i,j+1+n1*psize);
                        (*Ch)(j+n1*(psize-1),i) = DintCinvBTtilde_m_Dbnd(j+1+n1*psize,i);
                    }
                }
            }

            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(n2 = 0; n2 < nz_loc; ++n2)
                {
                    for(j = 0; j < psize-1; ++j)
                    {
                        (*Bh)(nsize_bndry_p1[n]-nz_loc+n1,j+n2*(psize-1)) = -DintCinvDTint(n1*psize,j+1+n2*psize);
                        (*Ch)(j+n2*(psize-1),nsize_bndry_p1[n]-nz_loc+n1) = -DintCinvDTint(j+1+n2*psize,n1*psize);
                    }
                }
            }

            // Do static condensation
            Dh->Invert();
            (*Bh) = (*Bh)*(*Dh);
            //(*Ah) = (*Ah) - (*Bh)*(*Ch);
            Blas::Dgemm('N','N', Bh->GetRows(), Ch->GetColumns(), Bh->GetColumns(), -1.0,
                        Bh->GetRawPtr(), Bh->GetRows(), Ch->GetRawPtr(), Ch->GetRows(),
                        1.0, Ah->GetRawPtr(), Ah->GetRows());

            // Set matrices for later inversion. Probably do not need to be
            // attached to class
            pAh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ah));
            pBh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Bh));
            pCh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ch));
            pDh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Dh));
        }
        timer.Stop();
        if (debug_mode)
        {
                cout << "Matrix Setup Costs: " << timer.TimePerTest(1) << endl;
	}

        timer.Start();
        // Set up global coupled boundary solver.
        // This is a key to define the solution matrix type
        // currently we are giving it a argument of eLInearAdvectionReaction
        // since this then makes the matrix storage of type eFull
        MultiRegions::GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,locToGloMap);
        mat.m_CoupledBndSys = MemoryManager<MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(key,m_fields[0],pAh,pBh,pCh,pDh,locToGloMap);
        mat.m_CoupledBndSys->Initialise(locToGloMap);
    }

    void CoupledLinearNS_ROM::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Solver Type", "Coupled Linearised NS");
    }

    void CoupledLinearNS_ROM::load_snapshots()
    {
    
    	cout << " Loading FOM snapshots from files ... " << endl;
    	
	// fill the fields snapshot_x_collection and snapshot_y_collection

	int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);  // number of phys points
        }
               
        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(m_session->GetVariable(i));
        }

	snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
	snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);

        for(int i = 0; i < number_of_snapshots; ++i)
        {
		// generate the correct string
		std::stringstream sstm;
		sstm << "TestSnap" << i+1;
		std::string result = sstm.str();
		const char* rr = result.c_str();

//	        EvaluateFunction(fieldStr, test_load_snapshot, result);
	        GetFunction(result)->Evaluate(fieldStr, test_load_snapshot);
	        
            //     now:   GetFunction( "Restart")->Evaluate(fieldStr,  Restart);
            //     PREV:    EvaluateFunction(fieldStr, Restart, "Restart");
	        
		snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		for (int j=0; j < GetNpoints(); ++j)
		{
			snapshot_x_collection[i][j] = test_load_snapshot[0][j];
			snapshot_y_collection[i][j] = test_load_snapshot[1][j];
		}
        }    	

	
    	cout << " ... finished loading FOM snapshots from files" << endl;
    	
    	if (use_fine_grid_VV_and_load_ref && use_fine_grid_VV)
    	{
    		cout << "start loading Verification & Validation snapshots ... " << endl;
		snapshot_x_collection_VV = Array<OneD, Array<OneD, NekDouble> > (fine_grid_dir0);
		snapshot_y_collection_VV = Array<OneD, Array<OneD, NekDouble> > (fine_grid_dir0);    	
	        for(int i = 0; i < fine_grid_dir0; ++i)
	        {
			// generate the correct string
			std::stringstream sstm;
			sstm << "VV" << i+1;
			std::string result = sstm.str();
			const char* rr = result.c_str();
		        GetFunction(result)->Evaluate(fieldStr, test_load_snapshot);
		        
			snapshot_x_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
			snapshot_y_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
			for (int j=0; j < GetNpoints(); ++j)
			{
				snapshot_x_collection_VV[i][j] = test_load_snapshot[0][j];
				snapshot_y_collection_VV[i][j] = test_load_snapshot[1][j];
			}
	        }    	
    		cout << " ... finished loading Verification & Validation snapshots " << endl;
    	}
    }
    
    void CoupledLinearNS_ROM::load_session_parameters()
    {
    
    	cout << " Loading ROM parameters ..." << endl;
    
    	ROM_started = 0;
    	ongoing_snapshot_computation = 0;
    	use_Newton = 0;
    
    	load_snapshot_data_from_files = m_session->GetParameter("load_snapshot_data_from_files");
	number_of_snapshots = m_session->GetParameter("number_of_snapshots");
	POD_tolerance = m_session->GetParameter("POD_tolerance");
	if (m_session->DefinesParameter("parameter_space_dimension")) 
	{
		parameter_space_dimension = m_session->GetParameter("parameter_space_dimension");	
	}
	else
	{
		parameter_space_dimension = 1;
	}
	if (m_session->DefinesParameter("debug_mode")) 
	{
		debug_mode = m_session->GetParameter("debug_mode");	
	}
	else
	{
		debug_mode = 1;
	}	
	
	if (m_session->DefinesParameter("do_trafo_check")) 
	{
		do_trafo_check = m_session->GetParameter("do_trafo_check");	
	}
	else
	{
		do_trafo_check = 1;
	}
	if (m_session->DefinesParameter("do_trafo_check_relative_error")) 
	{
		do_trafo_check_relative_error = m_session->GetParameter("do_trafo_check_relative_error");	
	}
	else
	{
		do_trafo_check_relative_error = 1e-9;
	}
	if (m_session->DefinesParameter("compute_smaller_model_errs")) 
	{
		compute_smaller_model_errs = m_session->GetParameter("compute_smaller_model_errs");	
	}
	else
	{
		compute_smaller_model_errs = 0;
	}
	if (m_session->DefinesParameter("use_fine_grid_VV")) 
	{
		use_fine_grid_VV = m_session->GetParameter("use_fine_grid_VV");	
	}
	else
	{
		use_fine_grid_VV = 0;
	}
	if (m_session->DefinesParameter("use_fine_grid_VV_and_load_ref")) 
	{
		use_fine_grid_VV_and_load_ref = m_session->GetParameter("use_fine_grid_VV_and_load_ref");	
	}
	else
	{
		use_fine_grid_VV_and_load_ref = 1;
	}	
	if (m_session->DefinesParameter("use_fine_grid_VV_random")) 
	{
		use_fine_grid_VV_random = m_session->GetParameter("use_fine_grid_VV_random");	
	}
	else
	{
		use_fine_grid_VV_random = 0;
	}
	if (m_session->DefinesParameter("use_sparse_poly")) 
	{
		use_sparse_poly = m_session->GetParameter("use_sparse_poly");
	}
	else
	{
		use_sparse_poly = 0;
	} 
	if (m_session->DefinesParameter("max_sparse_poly_approx_dimension")) 
	{
		max_sparse_poly_approx_dimension = m_session->GetParameter("max_sparse_poly_approx_dimension");
	}
	else
	{
		max_sparse_poly_approx_dimension = 1;
	} 
	if (m_session->DefinesParameter("Geo_trafo_load")) 
	{
		Geo_trafo_load = m_session->GetParameter("Geo_trafo_load");
	}
	else
	{
		Geo_trafo_load = 1;
	} 	
	if (m_session->DefinesParameter("replace_snapshot_with_transformed")) 
	{
		replace_snapshot_with_transformed = m_session->GetParameter("replace_snapshot_with_transformed");
	}
	else
	{
		replace_snapshot_with_transformed = 1;
	} 
	
	
	
	parameter_types = Array<OneD, int> (parameter_space_dimension); 
	parameter_types[0] = m_session->GetParameter("type_para1");
	Nmax = number_of_snapshots;
    	if (parameter_space_dimension == 1)
	{
		param_vector = Array<OneD, NekDouble> (Nmax);
		for(int i = 0; i < number_of_snapshots; ++i)
		{
			// generate the correct string
			std::stringstream sstm;
			sstm << "param" << i;
			std::string result = sstm.str();
			// const char* rr = result.c_str();
		        param_vector[i] = m_session->GetParameter(result);
	        }
	}
	if (parameter_space_dimension == 2)
	{
		parameter_types[1] = m_session->GetParameter("type_para2");
		general_param_vector = Array<OneD, Array<OneD, NekDouble> > (Nmax);
		number_of_snapshots_dir0 = m_session->GetParameter("number_of_snapshots_dir0");
		number_of_snapshots_dir1 = m_session->GetParameter("number_of_snapshots_dir1");
		int i_all = 0;
		Array<OneD, NekDouble> parameter_point(parameter_space_dimension, 0.0);
		for(int i = 0; i < Nmax; ++i)
		{
			parameter_point = Array<OneD, NekDouble> (parameter_space_dimension, 0.0);
			general_param_vector[i] = parameter_point;
		}
		for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
		{
			// generate the correct string
			std::stringstream sstm;
			sstm << "param" << i0 << "_dir0";
			std::string result = sstm.str();
		        if (i0 == 0)
				start_param_dir0 = m_session->GetParameter(result);
		        if (i0 == number_of_snapshots_dir0-1)
				end_param_dir0 = m_session->GetParameter(result);
			param_vector_dir0[i0] = m_session->GetParameter(result);
			for(int i1 = 0; i1 < number_of_snapshots_dir1; ++i1)
			{
				// generate the correct string
				std::stringstream sstm1;
				sstm1 << "param" << i1 << "_dir1";
				std::string result1 = sstm1.str();
			    if (i1 == 0)
					start_param_dir1 = m_session->GetParameter(result1);
			    if (i1 == number_of_snapshots_dir1-1)
					end_param_dir1 = m_session->GetParameter(result1);
				general_param_vector[i_all][0] = m_session->GetParameter(result);
			    general_param_vector[i_all][1] = m_session->GetParameter(result1);
				i_all = i_all + 1;
//				general_param_vector[i_all] = Array<OneD, NekDouble> (parameter_space_dimension);
				param_vector_dir1[i1] = m_session->GetParameter(result1);
			}
		}
		
		
		
	}
		if (use_fine_grid_VV && (parameter_space_dimension == 1))
		{
			fine_grid_dir0 = m_session->GetParameter("fine_grid_dir0");
			Nmax_VV = fine_grid_dir0;
			fine_general_param_vector = Array<OneD, Array<OneD, NekDouble> >(fine_grid_dir0);
			for(int i = 0; i < fine_grid_dir0; ++i)
			{
				Array<OneD, NekDouble> parameter_point = Array<OneD, NekDouble> (parameter_space_dimension, 0.0);
				fine_general_param_vector[i] = parameter_point;
			}

			// or alternatively overwrite with given random vector
			if (use_fine_grid_VV_random)
			{
				int current_index = 0;
				for (int i1 = 0; i1 < fine_grid_dir0; i1++)
				{
					// generate the correct string
					std::stringstream sstm;
					sstm << "VV_param" << i1 << "_dir0";
					std::string result = sstm.str();
					double param_dir0 = m_session->GetParameter(result);
	

					fine_general_param_vector[current_index][0] = param_dir0;
					current_index++;
 				}
				
			}

		}


	{ // currently just copy-pasted hard-coded from previous version	
		if (m_session->DefinesParameter("number_elem_trafo"))
		{
			number_elem_trafo = m_session->GetParameter("number_elem_trafo");
		}
		else
		{
			number_elem_trafo = 5;
		}
//		elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);

		switch(Geo_trafo_load) 
		{
			case 1:
				elements_trafo = set_elem_trafo( number_elem_trafo);
				break;
			case 2:
				elements_trafo = set_elem_trafo_no2( number_elem_trafo);
				break;
			case 3:
				elements_trafo = set_elem_trafo_no3( number_elem_trafo);
				break;
			case 4:
				elements_trafo = set_elem_trafo_no4( number_elem_trafo);
				break;
			case 5:
				elements_trafo = set_elem_trafo_no5( number_elem_trafo);
				break;
			default:
				cout << "error: geo trafo number not known " << endl;				
		}				




//	    <P> elem_1 = {32,30,11,10,9,8,17,16,15,14,25,24,23,22}   	</P>   
//          <P> elem_2 = {12, 28, 34, 13,21,20,18,19,26,27}   		</P>   
//	    <P> elem_3 = {29, 31}   	</P>   
//	    <P> elem_4 = {33, 35}   	</P>   
//	    <P> elem_5 = {0,1,2,3,4,5,6,7}  

/*		elements_trafo[0].insert(32); 
		elements_trafo[0].insert(30);
		elements_trafo[0].insert(11);
		elements_trafo[0].insert(10);
		elements_trafo[0].insert(9);
		elements_trafo[0].insert(8);
		elements_trafo[0].insert(17);
		elements_trafo[0].insert(16);
		elements_trafo[0].insert(15);
		elements_trafo[0].insert(14);
		elements_trafo[0].insert(25);
		elements_trafo[0].insert(24);
		elements_trafo[0].insert(23);
		elements_trafo[0].insert(22);

		elements_trafo[1].insert(12);
		elements_trafo[1].insert(28);
		elements_trafo[1].insert(34);
		elements_trafo[1].insert(13);
		elements_trafo[1].insert(21);
		elements_trafo[1].insert(20);
		elements_trafo[1].insert(18);
		elements_trafo[1].insert(19);
		elements_trafo[1].insert(26);
		elements_trafo[1].insert(27);

		elements_trafo[2].insert(29);
		elements_trafo[2].insert(31);

		elements_trafo[3].insert(33);
		elements_trafo[3].insert(35);

		elements_trafo[4].insert(0);
		elements_trafo[4].insert(1);
		elements_trafo[4].insert(2);
		elements_trafo[4].insert(3);
		elements_trafo[4].insert(4);
		elements_trafo[4].insert(5);
		elements_trafo[4].insert(6);
		elements_trafo[4].insert(7);  */
	}
	
    	cout << " ... finished loading ROM parameters" << endl;
	
    }
    
    Array<OneD, Array<OneD, NekDouble> > CoupledLinearNS_ROM::DoSolve_at_param(Array<OneD, NekDouble> init_snapshot_x, Array<OneD, NekDouble> init_snapshot_y, NekDouble parameter)
    {
//	DoInitialise();
//	DoSolve();
	bool snapshot_computation_plot_rel_errors = 1;
	double rel_err = 1.0;
	ongoing_snapshot_computation = 1;
	while (rel_err > 1e-10)
	{
		Set_m_kinvis( parameter );
		DoInitialiseAdv(init_snapshot_x, init_snapshot_y); // replaces .DoInitialise();
		DoSolve();
		// compare the accuracy
		Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = UpdateFields();
		m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys());
		m_fields_t[1]->BwdTrans(m_fields_t[1]->GetCoeffs(), m_fields_t[1]->UpdatePhys());
		Array<OneD, NekDouble> out_field_trafo_x(GetNpoints(), 0.0);
		Array<OneD, NekDouble> out_field_trafo_y(GetNpoints(), 0.0);

		Eigen::VectorXd csx0_trafo(GetNpoints());
		Eigen::VectorXd csy0_trafo(GetNpoints());
		Eigen::VectorXd csx0(GetNpoints());
		Eigen::VectorXd csy0(GetNpoints());

		CopyFromPhysField(0, out_field_trafo_x); 
		CopyFromPhysField(1, out_field_trafo_y);

		for( int index_conv = 0; index_conv < GetNpoints(); ++index_conv)
		{
			csx0_trafo(index_conv) = out_field_trafo_x[index_conv];
			csy0_trafo(index_conv) = out_field_trafo_y[index_conv];
			csx0(index_conv) = init_snapshot_x[index_conv];
			csy0(index_conv) = init_snapshot_y[index_conv];
		}

//		cout << "csx0.norm() " << csx0.norm() << endl;
//		cout << "csx0_trafo.norm() " << csx0_trafo.norm() << endl;
//		cout << "csy0.norm() " << csy0.norm() << endl;
//		cout << "csy0_trafo.norm() " << csy0_trafo.norm() << endl;
		rel_err = (csx0_trafo - csx0).norm() / csx0.norm() + (csy0_trafo - csy0).norm() / csy0.norm();
		if (snapshot_computation_plot_rel_errors)
		{
			cout << "rel_err euclidean norm " << rel_err << endl;
		}

		init_snapshot_x = out_field_trafo_x;
		init_snapshot_y = out_field_trafo_y;
	}



	Array<OneD, Array<OneD, NekDouble> > converged_solution( 2 );
	converged_solution[0] = Array<OneD, NekDouble>(GetNpoints(), 0.0);
	converged_solution[1] = Array<OneD, NekDouble>(GetNpoints(), 0.0);
	converged_solution[0] = init_snapshot_x;
	converged_solution[1] = init_snapshot_y;

//	cout << " curr_f_bnd.size()+curr_f_int.size() " <<  curr_f_bnd.size()+curr_f_int.size() << endl;
//	cout << " GetNcoeffs() " <<  GetNcoeffs() << endl;
	ongoing_snapshot_computation = 0;
	return converged_solution;
    }

    
    void CoupledLinearNS_ROM::compute_snapshots_kinvis(void)
    {
    	Array<OneD, NekDouble> zero_phys_init(GetNpoints(), 0.0);
	snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
	snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
        for(int i = 0; i < number_of_snapshots; ++i)
	{

		Array<OneD, Array<OneD, NekDouble> > converged_solution = DoSolve_at_param(zero_phys_init, zero_phys_init, param_vector[i]);
		snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		for (int j=0; j < GetNpoints(); ++j)
		{
			snapshot_x_collection[i][j] = converged_solution[0][j];
			snapshot_y_collection[i][j] = converged_solution[1][j];
		}
	}
    }
    
    double CoupledLinearNS_ROM::lagrange_interp(double curr_param, int curr_index, int sparse_poly_approx_dimension)
    {
	double lagrange_value = 1;
	for (int i = 0; i < sparse_poly_approx_dimension; ++i)
	{
		if (parameter_space_dimension == 1)
		{
			if (param_vector[curr_index] != param_vector[i])
				lagrange_value *= (curr_param - param_vector[i]) / (param_vector[curr_index] - param_vector[i]);
		}
	}
	return lagrange_value;
    }
    
    double CoupledLinearNS_ROM::lagrange_interp_tensorised_hierarchical(Array<OneD, double> curr_param, Array<OneD, int> curr_index)
    {
	double lagrange_value = 1;
//	cout << "curr_param[0] " << curr_param[0] << " curr_param[1] " << curr_param[1] << endl;
//	cout << "curr_index[0] " << curr_index[0] << " curr_index[1] " << curr_index[1] << endl; 
	for (int i = 0; i < curr_index[0]; ++i)
	{
		lagrange_value *= (curr_param[0] - param_vector_dir0[i]) / (param_vector_dir0[curr_index[0]] - param_vector_dir0[i]);
	}
	if (parameter_space_dimension > 1)
	{
		for (int i = 0; i < curr_index[1]; ++i)
		{
			lagrange_value *= (curr_param[1] - param_vector_dir1[i]) / (param_vector_dir1[curr_index[1]] - param_vector_dir1[i]);
		}
	}
	return lagrange_value;
    }    

    void CoupledLinearNS_ROM::compute_sparse_poly_approx()
    {
	int sparse_poly_approx_dimension = max_sparse_poly_approx_dimension;
	// L2 error works on the snapshot_x_collection and snapshot_y_collection
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double current_geo;
		if (parameter_space_dimension == 1)
		{
			if (parameter_types[0] == 0)
			{
				current_nu = param_vector[current_index];
			}
			if (parameter_types[0] == 1)
			{
				current_geo = param_vector[current_index];
				current_nu = m_kinvis;
			}
		}		
		Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].size());
		Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].size());
		for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
		{
			interpolant_x[i] = 0.0;
			interpolant_y[i] = 0.0;
		}
		for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
		{
			double lagrange_value = lagrange_interp(current_geo, index_interpol_op, sparse_poly_approx_dimension);
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)	
			{
				interpolant_x[i] += snapshot_x_collection[index_interpol_op][i] * lagrange_value;
				interpolant_y[i] += snapshot_y_collection[index_interpol_op][i] * lagrange_value;
			}
		}
		double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
		cout << "rel_L2error at parameter " << current_geo << " is " << rel_L2error << endl;
	}
	Array<OneD, NekDouble> collect_max(sparse_poly_approx_dimension);
	Array<OneD, NekDouble> collect_mean(sparse_poly_approx_dimension);
	double max, mean;
	for (int approx_dim = 1; approx_dim <= sparse_poly_approx_dimension; ++approx_dim)
	{
//		sparse_approx_VV(approx_dim, max, mean);
		cout << "max at dim " << approx_dim << " is " << max << endl;
		cout << "mean at dim " << approx_dim << " is " << mean << endl;
		collect_max[approx_dim-1] = max;
		collect_mean[approx_dim-1] = mean;
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_mean.txt";
	std::string sparse_conv_mean = sstm.str();
	const char* outname = sparse_conv_mean.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_mean[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_max.txt";
	std::string sparse_conv_max = sstm.str();
	const char* outname = sparse_conv_max.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_max[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
    }

	void CoupledLinearNS_ROM::compute_sparse_poly_approx_2D_lower_dim(int approx_dim, double &max_rel_L2, double &mean_rel_L2)
	{
		int sparse_poly_approx_dimension = approx_dim;
		Array<OneD,  Array<OneD, int> > index_set(max_sparse_poly_approx_dimension);
		for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
		{
			index_set[i] = Array<OneD, int>(parameter_space_dimension);
		}
		int added_indices = -1;
		int curr_sum = 0;
		while(added_indices < max_sparse_poly_approx_dimension)
		{
			for(int i=0; i <= curr_sum; ++i)
			{
				++added_indices;
				if (added_indices < max_sparse_poly_approx_dimension)
				{
					int index1 = i;
					int index2 = curr_sum - i;
					index_set[added_indices][0] = index1;
					index_set[added_indices][1] = index2;
				}
			}
			curr_sum++;
		}
		// compute the coefficients
		Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_x(max_sparse_poly_approx_dimension);
		Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_y(max_sparse_poly_approx_dimension);
		for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
		{
			sparse_poly_coefficients_x[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
			sparse_poly_coefficients_y[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
			for(int k=0; k<snapshot_x_collection[0].size(); ++k)
			{
				sparse_poly_coefficients_x[i][k] = 0;
				sparse_poly_coefficients_y[i][k] = 0;
			}
		}
		for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
		{
			// identify the correct snapshot index based on the index_set
			int index_all = index_set[i][1] + number_of_snapshots_dir1 * index_set[i][0]; 
			for(int k=0; k<snapshot_x_collection[0].size(); ++k)
			{
				sparse_poly_coefficients_x[i][k] =  snapshot_x_collection[index_all][k];
				sparse_poly_coefficients_y[i][k] =  snapshot_y_collection[index_all][k];
			}
			for (int j=0; j < i; ++j)
			{
				double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[index_all], index_set[j]);
				for(int k=0; k<snapshot_x_collection[0].size(); ++k)
				{
					sparse_poly_coefficients_x[i][k] -= sparse_poly_coefficients_x[j][k] * lith;
					sparse_poly_coefficients_y[i][k] -= sparse_poly_coefficients_y[j][k] * lith;
				}
			}
		}	
		Array<OneD, NekDouble> collect_L2(Nmax);
		// start sweeping 
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			int current_index = iter_index;
			double current_nu;
			double current_geo;
			Array<OneD, NekDouble> current_param(parameter_space_dimension, 0.0);
			current_param = general_param_vector[current_index];
			if (parameter_types[0] == 0)
			{
				current_nu = current_param[0];
				current_geo = current_param[1];
			}
			if (parameter_types[0] == 1)
			{
				current_nu = current_param[1];
				current_geo = current_param[0];
			}
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].size());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].size());
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
			{
				double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[current_index], index_set[index_interpol_op]);
				for (int i = 0; i < snapshot_x_collection[0].size(); ++i)	
				{
					interpolant_x[i] += sparse_poly_coefficients_x[index_interpol_op][i] * lith;
					interpolant_y[i] += sparse_poly_coefficients_y[index_interpol_op][i] * lith;
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
			collect_L2[iter_index] = rel_L2error;
		}
		mean_rel_L2 = 0;
		max_rel_L2 = 0;
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			if (collect_L2[iter_index] > max_rel_L2)
				max_rel_L2 = collect_L2[iter_index];
			mean_rel_L2 += collect_L2[iter_index] / Nmax;
		}
	}


    void CoupledLinearNS_ROM::compute_sparse_poly_approx_2D()
    {
	int sparse_poly_approx_dimension = max_sparse_poly_approx_dimension;
	// L2 error works on the snapshot_x_collection and snapshot_y_collection
	
	// need to decide the grid rule 
	Array<OneD,  Array<OneD, int> > index_set(max_sparse_poly_approx_dimension);
	
	for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
	{
		index_set[i] = Array<OneD, int>(parameter_space_dimension);
	}
	// 0,0
	// 1,0
	// 0,1
	// 1,1
	// 2,0
	// 0,2
	// 3,0 ..
	int added_indices = -1;
	int curr_sum = 0;
	while(added_indices < max_sparse_poly_approx_dimension)
	{
		for(int i=0; i <= curr_sum; ++i)
		{
			++added_indices;
			if (added_indices < max_sparse_poly_approx_dimension)
			{
				int index1 = i;
				int index2 = curr_sum - i;
//				cout << "index1 " << index1 << " index2 " << index2 << endl; 				
				index_set[added_indices][0] = index1;
				index_set[added_indices][1] = index2;
			}
		
		}
		curr_sum++;
	}
	
	
	
	// compute the coefficients
	Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_x(max_sparse_poly_approx_dimension);
	Array<OneD, Array<OneD, NekDouble> > sparse_poly_coefficients_y(max_sparse_poly_approx_dimension);
	for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
	{
		sparse_poly_coefficients_x[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
		sparse_poly_coefficients_y[i] = Array<OneD, NekDouble> (snapshot_x_collection[0].size());
		for(int k=0; k<snapshot_x_collection[0].size(); ++k)
		{
			sparse_poly_coefficients_x[i][k] = 0;
			sparse_poly_coefficients_y[i][k] = 0;
		}
	}
	for (int i=0; i < max_sparse_poly_approx_dimension; ++i)
	{
		// identify the correct snapshot index based on the index_set
		int index_all = index_set[i][1] + number_of_snapshots_dir1 * index_set[i][0]; 
		for(int k=0; k<snapshot_x_collection[0].size(); ++k)
		{
			sparse_poly_coefficients_x[i][k] =  snapshot_x_collection[index_all][k];
			sparse_poly_coefficients_y[i][k] =  snapshot_y_collection[index_all][k];
		}
		
		for (int j=0; j < i; ++j)
		{
//			cout << "current j: " << j << " current i: " << i << endl;
			double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[index_all], index_set[j]);
//			cout << "lith " << lith << endl;
			for(int k=0; k<snapshot_x_collection[0].size(); ++k)
			{
				sparse_poly_coefficients_x[i][k] -= sparse_poly_coefficients_x[j][k] * lith;
				sparse_poly_coefficients_y[i][k] -= sparse_poly_coefficients_y[j][k] * lith;
			}
		}
		



	}	

	Array<OneD, NekDouble> collect_L2(Nmax);
	
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double current_geo;
		Array<OneD, NekDouble> current_param(parameter_space_dimension, 0.0);
		current_param = general_param_vector[current_index];
		if (parameter_types[0] == 0)
		{
			current_nu = current_param[0];
			current_geo = current_param[1];
		}
		if (parameter_types[0] == 1)
		{
			current_nu = current_param[1];
			current_geo = current_param[0];
		}
		Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].size());
		Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].size());
		for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
		{
			interpolant_x[i] = 0.0;
			interpolant_y[i] = 0.0;
		}
		for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
		{
			double lith = lagrange_interp_tensorised_hierarchical(general_param_vector[current_index], index_set[index_interpol_op]);
			cout << "lith sweep " << lith << endl;
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)	
			{
				interpolant_x[i] += sparse_poly_coefficients_x[index_interpol_op][i] * lith;
				interpolant_y[i] += sparse_poly_coefficients_y[index_interpol_op][i] * lith;
//				interpolant_x[i] += sparse_poly_coefficients_x[index_interpol_op][i];
//				interpolant_y[i] += sparse_poly_coefficients_y[index_interpol_op][i];
			}
		}
		double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
		// cout << "rel_L2error at 2D parameter geo " << current_geo << " nu " << current_nu << " is " << rel_L2error << endl;
		collect_L2[iter_index] = rel_L2error;
	}
	double mean_rel_L2 = 0;
	double max_rel_L2 = 0;
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		if (collect_L2[iter_index] > max_rel_L2)
			max_rel_L2 = collect_L2[iter_index];
		mean_rel_L2 += collect_L2[iter_index] / Nmax;
	}
	cout << "mean_rel_L2 " << mean_rel_L2 << " max_rel_L2 " << max_rel_L2 << endl;
	Array<OneD, NekDouble> collect_max(sparse_poly_approx_dimension);
	Array<OneD, NekDouble> collect_mean(sparse_poly_approx_dimension);
	double max, mean;
	for (int approx_dim = 1; approx_dim <= sparse_poly_approx_dimension; ++approx_dim)
	{
		compute_sparse_poly_approx_2D_lower_dim(approx_dim, max_rel_L2, mean_rel_L2);
		cout << "max at dim " << approx_dim << " is " << max_rel_L2 << endl;
		cout << "mean at dim " << approx_dim << " is " << mean_rel_L2 << endl;
		collect_max[approx_dim-1] = max_rel_L2;
		collect_mean[approx_dim-1] = mean_rel_L2;
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_mean.txt";
	std::string sparse_conv_mean = sstm.str();
	const char* outname = sparse_conv_mean.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_mean[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
	{
	std::stringstream sstm;
	sstm << "sparse_conv_max.txt";
	std::string sparse_conv_max = sstm.str();
	const char* outname = sparse_conv_max.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
		{
			myfile << std::setprecision(17) << collect_max[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 
	}
    }

    
    

    void CoupledLinearNS_ROM::v_DoInitialise(void)
    {
    
    	
    	load_session_parameters();
    	if ( load_snapshot_data_from_files )
	{
		if (parameter_space_dimension == 1)
		{
			load_snapshots();
		}
		else if (use_sparse_poly)
		{
			load_snapshots();
		
		}
		else
		{
			cout << "ROM error, expected parameter_space_dim == 1 when loading snapshots" << endl;
			return;
		}
	}
	else
	{
		if  (parameter_space_dimension == 1)
		{
			if (parameter_types[0] == 0)
			{
				compute_snapshots_kinvis();
			}
			else if (parameter_types[0] == 1)
			{
				cout << "missing compute snapshots function for geometric parameters" << endl;
				return;
			}
		}
		else
		{
			cout << "ROM error, expected parameter_space_dimension == 1 when computing snapshots" << endl;
			return;
		}
	}
	
        switch(m_equationType)
        {
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            {
                //                LibUtilities::TimeIntegrationMethod intMethod;
                //                std::string TimeIntStr = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
                //                int i;
                //                for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
                //                {
                //                    if(boost::iequals(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr))
                //                    {
                //                        intMethod = (LibUtilities::TimeIntegrationMethod)i;
                //                        break;
                //                    }
                //                }
                //
                //                ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
                //
                //                m_integrationScheme = LibUtilities::GetTimeIntegrationWrapperFactory().CreateInstance(LibUtilities::TimeIntegrationMethodMap[intMethod]);
                
                // Could defind this from IncNavierStokes class? 
                m_ode.DefineOdeRhs(&CoupledLinearNS_ROM::EvaluateAdvection, this);

                m_ode.DefineImplicitSolve(&CoupledLinearNS_ROM::SolveUnsteadyStokesSystem,this);

                // Set initial condition using time t=0

                SetInitialConditions(0.0);

            }
        case eSteadyStokes:
            SetUpCoupledMatrix(0.0);
            break;
        case eSteadyOseen:
            {
                Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.size());
                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }

                ASSERTL0(m_session->DefinesFunction("AdvectionVelocity"),
                         "Advection Velocity section must be defined in "
                         "session file.");

                std::vector<std::string> fieldStr;
                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                }
                GetFunction("AdvectionVelocity")->Evaluate(fieldStr, AdvField);

                SetUpCoupledMatrix(0.0,AdvField,false);
            }
            break;
        case eSteadyNavierStokes:
            {
                m_session->LoadParameter("KinvisMin", m_kinvisMin);
                m_session->LoadParameter("KinvisPercentage", m_KinvisPercentage);
                m_session->LoadParameter("Tolerence", m_tol);
                m_session->LoadParameter("MaxIteration", m_maxIt);
                m_session->LoadParameter("MatrixSetUpStep", m_MatrixSetUpStep);
                m_session->LoadParameter("Restart", m_Restart);


                DefineForcingTerm();

                if (m_Restart == 1)
                {
                    ASSERTL0(m_session->DefinesFunction("Restart"),
                             "Restart section must be defined in session file.");

                    Array<OneD, Array<OneD, NekDouble> > Restart(m_velocity.size());
                    for(int i = 0; i < m_velocity.size(); ++i)
                    {
                        Restart[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                    }
                    std::vector<std::string> fieldStr;
                    for(int i = 0; i < m_velocity.size(); ++i)
                    {
                        fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                    }
                    GetFunction( "Restart")->Evaluate(fieldStr,  Restart);

                    for(int i = 0; i < m_velocity.size(); ++i)
                    {
                        m_fields[m_velocity[i]]->FwdTrans_IterPerExp(Restart[i], m_fields[m_velocity[i]]->UpdateCoeffs());
                    }
                    cout << "Saving the RESTART file for m_kinvis = "<< m_kinvis << " (<=> Re = " << 1/m_kinvis << ")" <<endl;
                }
                else //We solve the Stokes Problem
                {
                    
                    SetUpCoupledMatrix(0.0);						
                    m_initialStep = true;
                    m_counter=1;
                    //SolveLinearNS(m_ForcingTerm_Coeffs);
                    Solve();
                    m_initialStep = false;
                    cout << "Saving the Stokes Flow for m_kinvis = "<< m_kinvis << " (<=> Re = " << 1/m_kinvis << ")" <<endl;
                }
            }
            break;
        case eSteadyLinearisedNS:
            {
                SetInitialConditions(0.0);

                Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.size());
                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }

                ASSERTL0(m_session->DefinesFunction("AdvectionVelocity"),
                         "Advection Velocity section must be defined in "
                         "session file.");

                std::vector<std::string> fieldStr;
                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                }
                GetFunction("AdvectionVelocity")->Evaluate(fieldStr, AdvField);

                SetUpCoupledMatrix(m_lambda,AdvField,true);
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS_ROM");
        }
    }

    void CoupledLinearNS_ROM::EvaluateAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                            Array<OneD, Array<OneD, NekDouble> > &outarray,
                                            const NekDouble time)
    {
        // evaluate convection terms
        EvaluateAdvectionTerms(inarray,outarray,time);

        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, outarray, outarray, time);
        }
    }

    void CoupledLinearNS_ROM::SolveUnsteadyStokesSystem(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                                                    const NekDouble time,
                                                    const NekDouble aii_Dt)
    {
        int i;
        Array<OneD, Array< OneD, NekDouble> > F(m_nConvectiveFields);
        NekDouble  lambda = 1.0/aii_Dt;
        static NekDouble lambda_store;
        Array <OneD, Array<OneD, NekDouble> > forcing(m_velocity.size());
        // Matrix solution
        if(fabs(lambda_store - lambda) > 1e-10)
        {
            SetUpCoupledMatrix(lambda);
            lambda_store = lambda;
        }

        // Forcing for advection solve
        for(i = 0; i < m_velocity.size(); ++i)
        {
            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[m_velocity[i]]->SetWaveSpace(true);
            m_fields[m_velocity[i]]->IProductWRTBase(inarray[i],m_fields[m_velocity[i]]->UpdateCoeffs());
            m_fields[m_velocity[i]]->SetWaveSpace(waveSpace);
            Vmath::Smul(m_fields[m_velocity[i]]->GetNcoeffs(),lambda,m_fields[m_velocity[i]]->GetCoeffs(), 1,m_fields[m_velocity[i]]->UpdateCoeffs(),1);
            forcing[i] = m_fields[m_velocity[i]]->GetCoeffs();
        }

        SolveLinearNS(forcing);

        for(i = 0; i < m_velocity.size(); ++i)
        {
            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),outarray[i]);
        }
    }


    void CoupledLinearNS_ROM::v_TransCoeffToPhys(void)
    {
        int nfields = m_fields.size();
        for (int k=0 ; k < nfields; ++k)
        {
            //Backward Transformation in physical space for time evolution
            m_fields[k]->BwdTrans_IterPerExp(m_fields[k]->GetCoeffs(),
                                             m_fields[k]->UpdatePhys());
        }

    }

    void CoupledLinearNS_ROM::v_TransPhysToCoeff(void)
    {
        int nfields = m_fields.size();
        for (int k=0 ; k < nfields; ++k)
        {
            //Forward Transformation in physical space for time evolution
            m_fields[k]->FwdTrans_IterPerExp(m_fields[k]->GetPhys(),
                                             m_fields[k]->UpdateCoeffs());

        }
    }

    void CoupledLinearNS_ROM::v_DoSolve(void)
    {
        switch(m_equationType)
        {
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            //AdvanceInTime(m_steps);
            UnsteadySystem::v_DoSolve();
            break;
        case eSteadyStokes:
        case eSteadyOseen:
        case eSteadyLinearisedNS:
            {
                Solve();
                break;
            }
            case eSteadyNavierStokes:
            {
                LibUtilities::Timer Generaltimer;
                Generaltimer.Start();

                int Check(0);

                //Saving the init datas
                Checkpoint_Output(Check);
                Check++;

                cout<<"We execute INITIALLY SolveSteadyNavierStokes for m_kinvis = "<<m_kinvis<<" (<=> Re = "<<1/m_kinvis<<")"<<endl;
                SolveSteadyNavierStokes();

                while(m_kinvis > m_kinvisMin)
                {
                    if (Check == 1)
                    {
                        cout<<"We execute SolveSteadyNavierStokes for m_kinvis = "<<m_kinvis<<" (<=> Re = "<<1/m_kinvis<<")"<<endl;
                        SolveSteadyNavierStokes();
                        Checkpoint_Output(Check);
                        Check++;
                    }

                    Continuation();

                    if (m_kinvis > m_kinvisMin)
                    {
                        cout<<"We execute SolveSteadyNavierStokes for m_kinvis = "<<m_kinvis<<" (<=> Re = "<<1/m_kinvis<<")"<<endl;
                        SolveSteadyNavierStokes();
                        Checkpoint_Output(Check);
                        Check++;
                    }
                }


                Generaltimer.Stop();
                cout<<"\nThe total calculation time is : " << Generaltimer.TimePerTest(1)/60 << " minute(s). \n\n";

                break;
            }
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS_ROM");
        }
        
        if ((!ROM_started) && (!ongoing_snapshot_computation) && (!use_sparse_poly))
        {
        	ROM_started = 1;
        	ROM_offline_phase();
		ROM_online_phase();
        }
        if (use_sparse_poly)
        {
               	f_bnd_size = curr_f_bnd.size();
		f_p_size = curr_f_p.size();
		f_int_size = curr_f_int.size();
		collect_f_all = DoTrafo();
		if (parameter_space_dimension == 1)
		{
	        	compute_sparse_poly_approx();
	        }
	        if (parameter_space_dimension == 2)
		{
	        	compute_sparse_poly_approx_2D();
	        }
        }
    }
    
    void CoupledLinearNS_ROM::Set_m_kinvis(NekDouble input)
    {
	m_kinvis = input;
    }

    Eigen::MatrixXd CoupledLinearNS_ROM::DoTrafo()
    {
    
    	Eigen::MatrixXd transformed_snapshots;

	if (parameter_space_dimension == 1)
	{
	
		if (parameter_types[0] == 0)
		{
			// if it is a single kinvis parameter
			transformed_snapshots = DoTrafo_1p_kinvis();
		}
		if (parameter_types[0] == 1)
		{
			// if it is a single geometric parameter
			transformed_snapshots = DoTrafo_1p_geo();
		}
	}
	
	
	return transformed_snapshots;

    }
    
    Eigen::MatrixXd CoupledLinearNS_ROM::DoTrafo_1p_kinvis()
    {
    
    	cout << " transform snapshots to different format for single kinvis parameter ... " << endl;
    
	Eigen::MatrixXd collect_f_bnd( curr_f_bnd.size() , Nmax );
	Eigen::MatrixXd collect_f_p( curr_f_p.size() , Nmax );
	Eigen::MatrixXd collect_f_int( curr_f_int.size() , Nmax );
	for (int i=0; i<Nmax; i++)
	{
		Set_m_kinvis( param_vector[i] );	
		DoInitialiseAdv(snapshot_x_collection[i], snapshot_y_collection[i]); // replaces .DoInitialise();
		DoSolve();
	
	// compare the accuracy
		Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = UpdateFields();
		m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys());
		m_fields_t[1]->BwdTrans(m_fields_t[1]->GetCoeffs(), m_fields_t[1]->UpdatePhys());
		Array<OneD, NekDouble> out_field_trafo_x(GetNpoints(), 0.0);
		Array<OneD, NekDouble> out_field_trafo_y(GetNpoints(), 0.0);

		Eigen::VectorXd csx0_trafo(GetNpoints());
		Eigen::VectorXd csy0_trafo(GetNpoints());
		Eigen::VectorXd csx0(GetNpoints());
		Eigen::VectorXd csy0(GetNpoints());

		CopyFromPhysField(0, out_field_trafo_x); 
		CopyFromPhysField(1, out_field_trafo_y);
		for( int index_conv = 0; index_conv < GetNpoints(); ++index_conv)
		{
			csx0_trafo(index_conv) = out_field_trafo_x[index_conv];
			csy0_trafo(index_conv) = out_field_trafo_y[index_conv];
			csx0(index_conv) = snapshot_x_collection[i][index_conv];
			csy0(index_conv) = snapshot_y_collection[i][index_conv];
		}

		if (debug_mode)
		{
			cout << "csx0.norm() " << csx0.norm() << endl;
			cout << "csx0_trafo.norm() " << csx0_trafo.norm() << endl;
			cout << "csy0.norm() " << csy0.norm() << endl;
			cout << "csy0_trafo.norm() " << csy0_trafo.norm() << endl;
		}

		Eigen::VectorXd trafo_f_bnd = curr_f_bnd;
		Eigen::VectorXd trafo_f_p = curr_f_p;
		Eigen::VectorXd trafo_f_int = curr_f_int;

		collect_f_bnd.col(i) = trafo_f_bnd;
		collect_f_p.col(i) = trafo_f_p;
		collect_f_int.col(i) = trafo_f_int;
		
		
               // generate the correct string
               std::stringstream sstm;
               sstm << "Conv_Oseen_param" << i << ".fld";
               std::string filename = sstm.str();
               write_curr_field(filename);



		
	}    	
    	
    	Eigen::MatrixXd transformed_snapshots( curr_f_bnd.size()+curr_f_p.size()+curr_f_int.size() , Nmax );
	transformed_snapshots.block(0,0,collect_f_bnd.rows(),collect_f_bnd.cols()) = collect_f_bnd;
	transformed_snapshots.block(collect_f_bnd.rows(),0,collect_f_p.rows(),collect_f_p.cols()) = collect_f_p;
	transformed_snapshots.block(collect_f_bnd.rows()+collect_f_p.rows(),0,collect_f_int.rows(),collect_f_int.cols()) = collect_f_int;
    	
    	cout << " ... finished transform snapshots to different format " << endl;
    
    	
    	return transformed_snapshots;
    }
    
    void CoupledLinearNS_ROM::DoInitialiseAdv(Array<OneD, NekDouble> myAdvField_x, Array<OneD, NekDouble> myAdvField_y)
    {
	// only covers case eSteadyOseen

	// moved to .h	Array<OneD, Array<OneD, NekDouble> > myAdvField(2);
	myAdvField = Array<OneD, Array<OneD, NekDouble> > (2);
	myAdvField[0] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	myAdvField[1] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	Array<OneD, Array<OneD, NekDouble> > local_myAdvField(2);	
	local_myAdvField[0] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	local_myAdvField[1] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	for (int i = 0; i<m_fields[0]->GetTotPoints(); ++i)
	{
		myAdvField[0][i] = myAdvField_x[i];
		myAdvField[1][i] = myAdvField_y[i];
		local_myAdvField[0][i] = myAdvField_x[i];
		local_myAdvField[1][i] = myAdvField_y[i];

	}

        std::vector<std::string> fieldStr;
        for(int i = 0; i < m_velocity.size(); ++i)
        {
             fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
        }
//	cout << "fieldStr[0] " << fieldStr[0] << endl;
//        EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity"); // defined in EquationSystem

        SetUpCoupledMatrix(0.0, local_myAdvField, false);

    }
    
    
    
    void CoupledLinearNS_ROM::setDBC(Eigen::MatrixXd collect_f_all)
    {
	no_dbc_in_loc = 0;
	no_not_dbc_in_loc = 0;
	for ( int index_c_f_bnd = 0; index_c_f_bnd < curr_f_bnd.size(); index_c_f_bnd++ )
	{
		if (collect_f_all(index_c_f_bnd,0) == collect_f_all(index_c_f_bnd,1))
		{
			no_dbc_in_loc++;
			elem_loc_dbc.insert(index_c_f_bnd);
		}
		else
		{
			no_not_dbc_in_loc++;
			elem_not_loc_dbc.insert(index_c_f_bnd);
		}
	}
    }

    void CoupledLinearNS_ROM::set_f_bnd_dbc()
    {
	nBndDofs = m_locToGloMap[0]->GetNumGlobalBndCoeffs();  // number of global bnd dofs
        const Array<OneD,const int>& loctoglobndmap = m_locToGloMap[0]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign = m_locToGloMap[0]->GetLocalToGlobalBndSign();
//	Eigen::MatrixXd Mtrafo(RB_A.rows(), nBndDofs);
	M_truth_size = curr_f_bnd.size() + curr_f_p.size() + curr_f_int.size();  // compare_vec1.rows() corresponds to nBndDofs
	if (debug_mode)
	{
		cout << "Local dof size, also M_truth_size is " << curr_f_bnd.size() + curr_f_p.size() + curr_f_int.size() << endl;
	}
	M_truth_size_without_DBC = no_not_dbc_in_loc + curr_f_p.size() + curr_f_int.size();
/*	Mtrafo = Eigen::MatrixXd (f_bnd_size, nBndDofs);
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
	int nsize_bndry_p1 = loctoglobndmap.size() / nel;
	int nsize_bndry = nsize_bndry_p1-1;
	for (int curr_elem = 0; curr_elem < nel; curr_elem++)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		int cnt_no_pp = curr_elem*nsize_bndry;
		for ( int index_ele = 0; index_ele < nsize_bndry_p1; index_ele++ )
		{
			int gid1 = loctoglobndmap[cnt+index_ele];
			int sign1 = loctoglobndsign[cnt+index_ele];
			if ((gid1 >= 0) && (index_ele < nsize_bndry))
			{
				Mtrafo(cnt_no_pp + index_ele, gid1) = sign1;
			}
		}
	}  */
//	MtM = Mtrafo * Mtrafo.transpose();
	f_bnd_dbc = Eigen::VectorXd::Zero(no_dbc_in_loc);
	f_bnd_dbc_full_size = Eigen::VectorXd::Zero(PODmodes.rows());
	RB = Eigen::MatrixXd::Zero(PODmodes.rows() - no_dbc_in_loc, PODmodes.cols());
	int counter_all = 0;
	int counter_dbc = 0;
	for (int index=0; index < PODmodes.rows(); ++index)
	{
		if (!elem_loc_dbc.count(index))
		{
			RB.row(counter_all) = PODmodes.row(index);
			f_bnd_dbc_full_size(index) = 0;
			counter_all++;
		}
		else
		{
			f_bnd_dbc_full_size(index) = collect_f_all(index,0);
			f_bnd_dbc(counter_dbc) = collect_f_all(index,0);
			counter_dbc++;
		}
	}

    }



    void CoupledLinearNS_ROM::gen_phys_base_vecs()
    {
	int RBsize = RB.cols();
	PhysBaseVec_x = Array<OneD, Array<OneD, double> > (RBsize); 
	PhysBaseVec_y = Array<OneD, Array<OneD, double> > (RBsize);
	
	if (debug_mode)
	{
		cout << " number of local dofs per velocity direction " << GetNcoeffs() << endl;
		cout << " number of quadrature dofs per velocity direction " << GetNpoints() << endl;
	}

	for (int curr_trafo_iter=0; curr_trafo_iter < RBsize; curr_trafo_iter++)
	{
		Eigen::VectorXd f_bnd = PODmodes.block(0, curr_trafo_iter, curr_f_bnd.size(), 1);
		Eigen::VectorXd f_int = PODmodes.block(curr_f_bnd.size()+curr_f_p.size(), curr_trafo_iter, curr_f_int.size(), 1);
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
	        Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(GetNcoeffs());
		Array<OneD, double> field_1(GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
	        int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
	        int nz_loc = 1;
	        int  nplanecoeffs = fields[0]->GetNcoeffs();
		int  nel  = m_fields[0]->GetNumElmts();
	        for(int i = 0; i < nel; ++i) 
	        {
	            int eid  = i;
	            fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
	            fields[0]->GetExp(eid)->GetInteriorMap(imap);
	            int nbnd   = bmap.size();
	            int nint   = imap.size();
	            int offset = fields[0]->GetCoeff_Offset(eid);
	            
	            for(int j = 0; j < nvel; ++j)
	            {
	                for(int n = 0; n < nz_loc; ++n)
	                {
	                    for(int k = 0; k < nbnd; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
	                    }
	                    
	                    for(int k = 0; k < nint; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
	                    }
	                    cnt  += nbnd;
	                    cnt1 += nint;
	                }
	            }
	        }
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
		PhysBaseVec_x[curr_trafo_iter] = curr_PhysBaseVec_x;
		PhysBaseVec_y[curr_trafo_iter] = curr_PhysBaseVec_y;
		
	}

	eigen_phys_basis_x = Eigen::MatrixXd::Zero(GetNpoints(), RBsize);
	eigen_phys_basis_y = Eigen::MatrixXd::Zero(GetNpoints(), RBsize);
	for (int index_phys_base=0; index_phys_base<GetNpoints(); index_phys_base++)
	{
		for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
		{
			eigen_phys_basis_x(index_phys_base,index_RBsize) = PhysBaseVec_x[index_RBsize][index_phys_base];
			eigen_phys_basis_y(index_phys_base,index_RBsize) = PhysBaseVec_y[index_RBsize][index_phys_base];
		}
	}

	Eigen::VectorXd curr_col = eigen_phys_basis_x.col(0);
	double norm_curr_col = curr_col.norm();
	eigen_phys_basis_x.col(0) = curr_col / norm_curr_col;
	curr_col = eigen_phys_basis_y.col(0);
	norm_curr_col = curr_col.norm();
	eigen_phys_basis_y.col(0) = curr_col / norm_curr_col;
	Eigen::VectorXd orthogonal_complement;
	
	for (int orth_iter=1; orth_iter<RBsize; orth_iter++)
	{
		curr_col = eigen_phys_basis_x.col(orth_iter);
		Eigen::MatrixXd leftmostCols = eigen_phys_basis_x.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_x.col(orth_iter) = orthogonal_complement / norm_curr_col;
		curr_col = eigen_phys_basis_y.col(orth_iter);
		leftmostCols = eigen_phys_basis_y.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_y.col(orth_iter) = orthogonal_complement / norm_curr_col;
	}

	orth_PhysBaseVec_x = Array<OneD, Array<OneD, double> > (RBsize); 
	orth_PhysBaseVec_y = Array<OneD, Array<OneD, double> > (RBsize); 
	for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
	{
		Array<OneD, double> curr_iter_x(GetNpoints());
		Array<OneD, double> curr_iter_y(GetNpoints());
		for (int index_phys_base=0; index_phys_base<GetNpoints(); index_phys_base++)	
		{
			curr_iter_x[index_phys_base] = eigen_phys_basis_x(index_phys_base,index_RBsize);
			curr_iter_y[index_phys_base] = eigen_phys_basis_y(index_phys_base,index_RBsize);			
		}
		orth_PhysBaseVec_x[index_RBsize] = curr_iter_x;
		orth_PhysBaseVec_y[index_RBsize] = curr_iter_y;			
	}
    }
    

    void CoupledLinearNS_ROM::gen_proj_adv_terms()
    {
	RBsize = RB.cols();
	adv_mats_proj_x = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_mats_proj_y = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_x = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_y = Array<OneD, Eigen::VectorXd > (RBsize);
//	adv_vec_proj_x_newton = Array<OneD, Eigen::VectorXd > (RBsize);
//	adv_vec_proj_y_newton = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_x_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_y_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);


	Array<OneD, double> PhysBase_zero(GetNpoints(), 0.0);
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];

	//	InitObject();

//		DoInitialiseAdv(curr_PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
		// needs to be replaced with a more gen. term		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(RB_A.rows() + RB_Dbnd.rows() + RB_C.cols(), RB_A.cols() + RB_Dbnd.rows() + RB_B.cols() );
		Eigen::MatrixXd adv_matrix;
//		adv_matrix = Get_advection_matrix();      // <-- replace here with function -------------------------------------------------------------------------------------
		adv_matrix = gen_adv_mats_proj_x(curr_PhysBaseVec_x, use_Newton);
		Eigen::VectorXd add_to_rhs_adv(M_truth_size); // probably need this for adv and non-adv
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;
		Eigen::MatrixXd adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);

		adv_vec_proj_x_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);


		if (use_Newton)
		{
//			Eigen::VectorXd add_to_rhs_adv_newton(M_truth_size); 
//			add_to_rhs_adv_newton = adv_matrix * PODmodes * PODmodes.transpose() * collect_f_all.col(3);      
//			Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton, elem_loc_dbc);
			// alt: not working
//			adv_rhs_add_newton = adv_matrix_simplified * remove_rows(collect_f_all.col(3), elem_loc_dbc);
			// end alt
//			Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
//			adv_vec_proj_x_newton[trafo_iter] = adv_rhs_proj_newton;

			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
//				add_to_rhs_adv_newton_RB = adv_matrix_simplified * RB.col(RB_counter);
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;

				adv_vec_proj_x_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}


		Eigen::VectorXd adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		Eigen::MatrixXd adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		Eigen::VectorXd adv_rhs_proj = RB.transpose() * adv_rhs_add;

		adv_mats_proj_x[trafo_iter] = adv_mat_proj;
		adv_vec_proj_x[trafo_iter] = adv_rhs_proj;

		adv_vec_proj_y_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);
//		DoInitialiseAdv(PhysBase_zero , curr_PhysBaseVec_y ); // call with parameter in phys state
//		adv_matrix = Get_advection_matrix();      // <-- replace here with function -------------------------------------------------------------------------------------
		adv_matrix = gen_adv_mats_proj_y(curr_PhysBaseVec_y, use_Newton);
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   
		adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);
		adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		adv_rhs_proj = RB.transpose() * adv_rhs_add;
		adv_mats_proj_y[trafo_iter] = adv_mat_proj;
		adv_vec_proj_y[trafo_iter] = adv_rhs_proj;

		if (use_Newton)
		{
//			Eigen::VectorXd add_to_rhs_adv_newton(M_truth_size); 
//			add_to_rhs_adv_newton = adv_matrix  * PODmodes * PODmodes.transpose() *  collect_f_all.col(3);      
//			Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton, elem_loc_dbc);
//			Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
//			adv_vec_proj_y_newton[trafo_iter] = adv_rhs_proj_newton;


			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
//				add_to_rhs_adv_newton_RB = adv_matrix_simplified * RB.col(RB_counter);
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;

				adv_vec_proj_y_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}
	}
    }
    
    
    
    Eigen::MatrixXd CoupledLinearNS_ROM::gen_adv_mats_proj_x(Array<OneD, double> curr_PhysBaseVec_x, int use_Newton)
    {
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 
	Array<OneD, Eigen::MatrixXd > A_elem(m_fields[0]->GetNumElmts()); 
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int

        // Calculate derivative of base flow 
        Array<OneD, Array<OneD, NekDouble> > Advfield(m_velocity.size());
        for(int i = 0; i < m_velocity.size(); ++i)
        {
               Advfield[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0); // use here the input vector
	}
	Advfield[0] = curr_PhysBaseVec_x;

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
//		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
		int npoints = locExp->GetTotPoints();
	//	int eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(curr_elem);
		int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(curr_elem);
		Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
	        if(use_Newton) // formerly isLinearNSEquation
	        {
                    int nv1;
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(int nv = 0; nv < nvel; ++nv)
                    {
                        for(nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            
//                            if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
  //                          {
    //                            Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
      //                      }
        //                    else
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                            }
                        }
                    }
	        }


                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.size();
                int nimap = imap.size();
		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = curr_PhysBaseVec_x[curr_elem*nphys + i];
		}
		Array<OneD, double> Ah_ele_vec = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_x_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        Ah_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nbmap)*Ahrows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        C_ele_vec[i+(nv*nz_loc+n1)*nbmap + (j+(k*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[imap[j]];
                                    }
                                }
                            } 
			 }
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_x_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            int n1;
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u'.Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys, coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        B_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        D_ele_vec[j+(k*nz_loc+n1)*nimap + (i+(nv*nz_loc+n1)*nimap)*nsize_int] += coeffs[imap[j]];
                                    }
                                }
                            }
                        }

			} // for (int k = 0; k < 2; ++k)

		} // for (int i = 0; i < nimap; ++i)

		A_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_bndry );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				A_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		}
		D_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_int , nsize_int ); 
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		A_adv_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = A_elem[i];
		B_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = B_elem[i];
		C_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = C_elem[i];
		D_adv_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = D_elem[i];
	}

/*	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_adv_all * Mtrafo;
			adv_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_adv_all;
			adv_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_adv_all.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 2:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
	}  */
	adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
	adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
	adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
	adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;

	return adv_matrix;

	
    }

    Eigen::MatrixXd CoupledLinearNS_ROM::gen_adv_mats_proj_y(Array<OneD, double> curr_PhysBaseVec_y, int use_Newton)
    {
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 
	Array<OneD, Eigen::MatrixXd > A_elem(m_fields[0]->GetNumElmts()); 
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int

        // Calculate derivative of base flow 
        Array<OneD, Array<OneD, NekDouble> > Advfield(m_velocity.size());
        for(int i = 0; i < m_velocity.size(); ++i)
        {
               Advfield[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0); // use here the input vector
	}
	Advfield[1] = curr_PhysBaseVec_y;

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
//		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
		int npoints = locExp->GetTotPoints();
	//	int eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(curr_elem);
		int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(curr_elem);
		Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
	        if(use_Newton) // formerly isLinearNSEquation
	        {
                    int nv1;
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(int nv = 0; nv < nvel; ++nv)
                    {
                        for(nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            
//                            if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
  //                          {
    //                            Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
      //                      }
        //                    else
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                            }
                        }
                    }
	        }


                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.size();
                int nimap = imap.size();
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_y_part[i] = curr_PhysBaseVec_y[curr_elem*nphys + i];
		}
		Array<OneD, double> Ah_ele_vec = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        Ah_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nbmap)*Ahrows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        C_ele_vec[i+(nv*nz_loc+n1)*nbmap + (j+(k*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[imap[j]];
                                    }
                                }
                            } 
			 }
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_y_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            int n1;
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u'.Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys, coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        B_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        D_ele_vec[j+(k*nz_loc+n1)*nimap + (i+(nv*nz_loc+n1)*nimap)*nsize_int] += coeffs[imap[j]];
                                    }
                                }
                            }
                        }

			} // for (int k = 0; k < 2; ++k)

		} // for (int i = 0; i < nimap; ++i)

		A_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_bndry );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				A_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		}
		D_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_int , nsize_int ); 
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		A_adv_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = A_elem[i];
		B_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = B_elem[i];
		C_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = C_elem[i];
		D_adv_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = D_elem[i];
	}

/*	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_adv_all * Mtrafo;
			adv_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_adv_all;
			adv_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_adv_all.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 2:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
	} */
	adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
	adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
	adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
	adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;

	return adv_matrix;


    }
    

    Eigen::MatrixXd CoupledLinearNS_ROM::remove_cols_and_rows(Eigen::MatrixXd the_matrix, std::set<int> elements_to_be_removed)
    {

	// or move whole blocks in-place ??
	// need to use iterators
	Eigen::MatrixXd simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows() - elements_to_be_removed.size(), the_matrix.cols() - elements_to_be_removed.size());
	std::set<int>::iterator set_iterator_rows;
	std::set<int>::iterator set_iterator_cols; 
	int prev_iter_row = -1; // check if zero is to be removed
	int prev_iter_col = -1; // also need to move the last block??
	int passed_rows = 0;
	int passed_cols = 0;
//	Eigen::MatrixXd sm2;
	if (1)
	{
//		Eigen::MatrixXd simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows() - elements_to_be_removed.size(), the_matrix.cols() - elements_to_be_removed.size());
		for (set_iterator_rows = elements_to_be_removed.begin(); set_iterator_rows != elements_to_be_removed.end(); set_iterator_rows++)
		{
			if ((prev_iter_row != *set_iterator_rows) && (0 != *set_iterator_rows))
			{
				prev_iter_col = -1;
				passed_cols = 0;
				for (set_iterator_cols = elements_to_be_removed.begin(); set_iterator_cols != elements_to_be_removed.end(); set_iterator_cols++)
				{
//					cout << " output set_iterator_rows " << *set_iterator_rows << endl; // might contain many consecutive elements
//					cout << " output set_iterator_cols " << *set_iterator_cols << endl; 
					if ((prev_iter_col != *set_iterator_cols) && (0 != *set_iterator_cols))
					{
						// do move a block
						simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, *set_iterator_rows - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, *set_iterator_rows - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1);
					}
					prev_iter_col = *set_iterator_cols;
					passed_cols++;
				}
				// last block
				simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, *set_iterator_rows - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, *set_iterator_rows - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1);
			}
			prev_iter_row = *set_iterator_rows;
			passed_rows++;
		}
		prev_iter_col = -1;
		passed_cols = 0;
		for (set_iterator_cols = elements_to_be_removed.begin(); set_iterator_cols != elements_to_be_removed.end(); set_iterator_cols++)
		{
			if ((prev_iter_col != *set_iterator_cols) && (0 != *set_iterator_cols))
			{
				// do move a block
				simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, the_matrix.rows() - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, the_matrix.rows() - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1);
			}
			prev_iter_col = *set_iterator_cols;
			passed_cols++;
		}
		// last block
		simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, the_matrix.rows() - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, the_matrix.rows() - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1);
	}
	return simplified_matrix; 
    }

    Eigen::VectorXd CoupledLinearNS_ROM::remove_rows(Eigen::VectorXd the_vector, std::set<int> elements_to_be_removed)
    {
	Eigen::VectorXd simplified_vector = Eigen::VectorXd::Zero(the_vector.rows() - elements_to_be_removed.size());
	int counter_row_simplified = 0;
	for (int row_index=0; row_index < the_vector.rows(); ++row_index)
	{
		if (!elements_to_be_removed.count(row_index))
		{
			simplified_vector(counter_row_simplified) = the_vector(row_index);
			counter_row_simplified++;
		}		
	}
	return simplified_vector;
    }    
    
    Eigen::MatrixXd CoupledLinearNS_ROM::gen_no_advection_matrix_ABCD()
    {
//	double mKinvis = 1;
//	double w = 1;
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;

	Array<OneD, Eigen::MatrixXd > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();

                int nbmap = bmap.size();
                int nimap = imap.size();
		Array<OneD, double> Ah_ele_vec(Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec(nsize_int*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_0_0[int(bmap[j])] + coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_0_0[int(imap[j])] + coeffs_1_1[int(imap[j])];
				}
			} //for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					C_ele_vec[ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_0_0[int(bmap[j])] + coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_0_0[int(imap[j])] + coeffs_1_1[int(imap[j])];
				}
			} //for (int k = 0; k < 2; ++k)
		}

 		// chosen choice: redo the nektar++ approach
		// possible alternatives:
		// build instead the sing_* matrices directly
		// or copy for test purposes from setupcoupledmats??

		Ah_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 



	}

	Eigen::MatrixXd A_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	Eigen::MatrixXd B_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd D_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	Eigen::MatrixXd ABCD_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		A_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = Ah_elem[i].block(0,0,nsize_bndry,nsize_bndry);
		B_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = B_elem[i];
		C_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = C_elem[i];
		D_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = D_elem[i];
	}
/*	switch(globally_connected) {
		case 0:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 1:
			ABCD_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_all * Mtrafo;
			ABCD_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_all;
			ABCD_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_all.transpose() * Mtrafo;
			ABCD_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 2:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
	} */

	ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_all;
	ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_all;
	ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
	ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;


/*	Eigen::VectorXd add_to_rhs_ABCD(M_truth_size);
	add_to_rhs_ABCD = ABCD_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd ABCD_rhs_add = remove_rows(add_to_rhs_ABCD, elem_loc_dbc);
	ABCD_vec_proj = RB.transpose() * ABCD_rhs_add;
	Eigen::MatrixXd ABCD_matrix_simplified = remove_cols_and_rows(ABCD_matrix, elem_loc_dbc);
	Eigen::MatrixXd ABCD_mat_proj = RB.transpose() * ABCD_matrix_simplified * RB;
*/

	return ABCD_matrix;



    }    
    
    
    Eigen::MatrixXd CoupledLinearNS_ROM::gen_no_advection_matrix_pressure()
    {
//	double mKinvis = 1;
//	double w = 1;
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 

	// verify and transform to bnd / p / int the snapshot data

        int nz_loc;
        nz_loc = 1;

        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;

	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
	

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();

//		cout << "ncoeffs " << ncoeffs << endl;
//		cout << "nphys " << nphys << endl;
//		cout << "pqsize " << pqsize << endl;   // pqsize == nphys and ncoeffs == nphys / 2 when?

                int nbmap = bmap.size();
                int nimap = imap.size();
//		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
//		Array<OneD, double> curr_snap_y_part(nphys, 0.0);

		Array<OneD, double> Dbnd_ele_vec(nsize_p*nsize_bndry, 0.0);
		Array<OneD, double> Dint_ele_vec(nsize_p*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
			            	int psize = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		}

		Dbnd_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_bndry);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				Dbnd_elem[curr_elem](i,j) = Dbnd_ele_vec[ i + j*nsize_p];
			}
		} 
		Dint_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_int);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				Dint_elem[curr_elem](i,j) = Dint_ele_vec[ i + j*nsize_p];
			}
		} 
	}
//	cout << " nsize_p " << nsize_p << endl;
//	cout << " nsize_bndry " << nsize_bndry << endl;
//	cout << "Dbnd_elem[0] 1..4 " << Dbnd_elem[0].block(0,0,3,3) << endl;
	Eigen::MatrixXd Dbnd_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_bndry*nel );
	Eigen::MatrixXd Dint_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_int*nel );
	Eigen::MatrixXd press_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		Dbnd_all.block(i*nsize_p, i*nsize_bndry, nsize_p, nsize_bndry) = Dbnd_elem[i];
		Dint_all.block(i*nsize_p, i*nsize_int, nsize_p, nsize_int) = Dint_elem[i];
	}

/*	switch(globally_connected) {
		case 0:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -MtM * Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 1:
			press_matrix.block(0, nBndDofs, nBndDofs, f_p_size) = -Mtrafo.transpose() * Dbnd_all.transpose();
			press_matrix.block(nBndDofs, 0, f_p_size, nBndDofs) = -Dbnd_all * Mtrafo;
			press_matrix.block(nBndDofs, nBndDofs + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(nBndDofs + f_p_size, nBndDofs, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 2:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
	}	*/

	press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
	press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
	press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
	press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();


	// also need to create appropriate vector right-hand-side contribution

/*	Eigen::VectorXd add_to_rhs_press(M_truth_size); // probably need this for adv and non-adv
	add_to_rhs_press = press_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd press_rhs_add = remove_rows(add_to_rhs_press, elem_loc_dbc);
	press_vec_proj = RB.transpose() * press_rhs_add;
	Eigen::MatrixXd press_matrix_simplified = remove_cols_and_rows(press_matrix, elem_loc_dbc);
	Eigen::MatrixXd press_mat_proj = RB.transpose() * press_matrix_simplified * RB;
*/
	return press_matrix;

    }


    
    void CoupledLinearNS_ROM::gen_reference_matrices()
    {
//	double current_nu = ref_param_nu;
//	int current_index = ref_param_index;
//	Set_m_kinvis( current_nu );
	Set_m_kinvis( 1.0 );
//	DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]); // why is this necessary? 
//	the_const_one = Get_no_advection_matrix_pressure();
	the_const_one = gen_no_advection_matrix_pressure();
//	cout << "co norm " << the_const_one.block(0, f_bnd_size, 10, 10) << endl;
//	cout << "co2 norm " << the_const_one2.block(0, f_bnd_size, 10, 10) << endl;
//	Eigen::MatrixXd diff = the_const_one - the_const_one2;
//	cout << "diff norm " << diff.block(0, f_bnd_size, 10, 10) << endl;
//	the_ABCD_one = Get_no_advection_matrix_ABCD();
	the_ABCD_one = gen_no_advection_matrix_ABCD();
	the_const_one_simplified = remove_cols_and_rows(the_const_one, elem_loc_dbc);
	the_ABCD_one_simplified = remove_cols_and_rows(the_ABCD_one, elem_loc_dbc);
	the_const_one_proj = RB.transpose() * the_const_one_simplified * RB;
	the_ABCD_one_proj = RB.transpose() * the_ABCD_one_simplified * RB;
	the_ABCD_one_rhs = the_ABCD_one * f_bnd_dbc_full_size;
	the_const_one_rhs = the_const_one * f_bnd_dbc_full_size;
	the_ABCD_one_rhs_simplified = remove_rows(the_ABCD_one_rhs, elem_loc_dbc);
	the_const_one_rhs_simplified = remove_rows(the_const_one_rhs, elem_loc_dbc);
	the_ABCD_one_rhs_proj = RB.transpose() * the_ABCD_one_rhs_simplified;
	the_const_one_rhs_proj = RB.transpose() * the_const_one_rhs_simplified;
//	cout << "the_const_one_rhs_proj " << the_const_one_rhs_proj << endl;
    }


    void CoupledLinearNS_ROM::ROM_offline_phase()
    {
       	f_bnd_size = curr_f_bnd.size();
	f_p_size = curr_f_p.size();
	f_int_size = curr_f_int.size();

	collect_f_all = DoTrafo();


	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_all(collect_f_all, Eigen::ComputeThinU);
//	cout << "svd_collect_f_all.singularValues() " << svd_collect_f_all.singularValues() << endl << endl;
	Eigen::VectorXd singular_values = svd_collect_f_all.singularValues();
	if (debug_mode)
	{
		cout << "sum singular values " << singular_values.sum() << endl << endl;
	}
	Eigen::VectorXd rel_singular_values = singular_values / singular_values.sum();
//	cout << "relative singular value percents: " << rel_singular_values << endl;
	// determine RBsize corresponding to the chosen POD_tolerance
	RBsize = 1; 
	Eigen::VectorXd cum_rel_singular_values = Eigen::VectorXd::Zero(singular_values.rows());
	for (int i = 0; i < singular_values.rows(); ++i)
	{
		cum_rel_singular_values(i) = singular_values.head(i+1).sum() / singular_values.sum();
		if (cum_rel_singular_values(i) < POD_tolerance)
		{
			RBsize = i+2;
		}		
	}
	if (debug_mode)
	{
		cout << "cumulative relative singular value percentages: " << cum_rel_singular_values << endl;
		cout << "RBsize: " << RBsize << endl;
	}
	
	Eigen::MatrixXd collect_f_all_PODmodes = svd_collect_f_all.matrixU(); // this is a local variable...
	setDBC(collect_f_all); // agnostic to RBsize
	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), RBsize);  
	PODmodes = collect_f_all_PODmodes.leftCols(RBsize);
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
	set_f_bnd_dbc();
	if (debug_mode)
	{
//		cout << "M_no_dbc_in_loc " << M_no_dbc_in_loc << endl;
		cout << "no_dbc_in_loc " << no_dbc_in_loc << endl;
//		cout << "M_no_not_dbc_in_loc " << M_no_not_dbc_in_loc << endl;
		cout << "no_not_dbc_in_loc " <<	no_not_dbc_in_loc << endl;
	}	
	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc = RB;
	if (debug_mode)
	{
		cout << "c_f_all_PODmodes_wo_dbc.rows() " << c_f_all_PODmodes_wo_dbc.rows() << endl;
		cout << "c_f_all_PODmodes_wo_dbc.cols() " << c_f_all_PODmodes_wo_dbc.cols() << endl;
	}
	gen_phys_base_vecs();
	cout << "finished gen_phys_base_vecs " << endl;
	if (parameter_space_dimension == 1)
	{
		if (parameter_types[0] == 0)
		{
			gen_proj_adv_terms();
			cout << "finished gen_proj_adv_terms " << endl;
		}
		if (parameter_types[0] == 1)
		{
			gen_proj_adv_terms_geo();
			cout << "finished gen_proj_adv_terms " << endl;
		}		
	}
	if (parameter_space_dimension == 1)
	{
		if (parameter_types[0] == 0)
		{
			gen_reference_matrices();
			cout << "finished gen_reference_matrices " << endl;
		}
		if (parameter_types[0] == 1)
		{
			gen_reference_matrices_geo();
			cout << "finished gen_reference_matrices " << endl;
		}		
	}
    }



    void CoupledLinearNS_ROM::ROM_online_phase()
    {
    	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);  // is of size M_truth_size
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double current_geo;
		if (parameter_space_dimension == 1)
		{
			
			if (parameter_types[0] == 0)
			{
				current_nu = param_vector[current_index];
			}
			if (parameter_types[0] == 1)
			{
				current_geo = param_vector[current_index];
				current_nu = m_kinvis;
			}				
			
		}
		if (debug_mode)
		{
			cout << " online phase current nu " << current_nu << endl;
			cout << " online phase current geo " << current_geo << endl;
		}
		Set_m_kinvis( current_nu );
		if (use_Newton)
		{
			DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		}
		Eigen::MatrixXd curr_xy_proj = project_onto_basis(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		Eigen::MatrixXd affine_mat_proj;
		Eigen::VectorXd affine_vec_proj;
		if (parameter_space_dimension == 1)
		{
			if (parameter_types[0] == 0)
			{
				affine_mat_proj = gen_affine_mat_proj(current_nu);
				affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
			}
			if (parameter_types[0] == 1)
			{
				affine_mat_proj = gen_affine_mat_proj_geo(current_nu, current_geo);
				affine_vec_proj = gen_affine_vec_proj_geo(current_nu, current_geo, current_index);
			}

		}
		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
		Eigen::VectorXd repro_solve_affine = RB * solve_affine;
		Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
		mat_compare.col(0) = collect_f_all.col(current_index);
		if (debug_mode)
		{
			Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
			current_f_all = collect_f_all.col(current_index);
			Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
			Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
			// cout << "proj_current_f_all_wo_dbc " << proj_current_f_all_wo_dbc << endl;
			Eigen::VectorXd correctRHS = affine_mat_proj * proj_current_f_all_wo_dbc;
			Eigen::VectorXd correction_RHS = correctRHS - affine_vec_proj;
		}
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) - mat_compare.col(0);
		Eigen::VectorXd FOM_solution = mat_compare.col(0);
		Eigen::VectorXd FOM_solution_wo_dbc = remove_rows(FOM_solution, elem_loc_dbc);
		Eigen::VectorXd proj_FOM_solution_wo_dbc = RB.transpose() * FOM_solution_wo_dbc;
		Eigen::VectorXd reproj_FOM_solution_wo_dbc = RB * proj_FOM_solution_wo_dbc;
		Eigen::VectorXd reconstruct_FOM_solution = reconstruct_solution_w_dbc(reproj_FOM_solution_wo_dbc);
		Eigen::VectorXd diff_projection = reconstruct_FOM_solution - mat_compare.col(0);
		// now only in RB:
		Eigen::VectorXd diff_projection_RB = reproj_FOM_solution_wo_dbc - remove_rows(mat_compare.col(0), elem_loc_dbc);
		Eigen::VectorXd diff_RB = repro_solve_affine - remove_rows(mat_compare.col(0), elem_loc_dbc);
		if (debug_mode)
		{
			cout << "snapshot_x_collection.size() " << snapshot_x_collection.size() << " snapshot_x_collection[0].size() " << snapshot_x_collection[0].size() << endl;
		}

//		if (write_ROM_field || (qoi_dof >= 0))
//		{
//			recover_snapshot_data(reconstruct_solution, current_index); // this is setting the fields and fieldcoeffs
//		}
	  	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
		Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
		Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(GetNcoeffs());
		Array<OneD, double> field_1(GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
		int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
		int nz_loc = 1;
		int  nplanecoeffs = fields[0]->GetNcoeffs();
		int  nel  = m_fields[0]->GetNumElmts();
		for(int i = 0; i < nel; ++i) 
		{
		      int eid  = i;
		      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
		      fields[0]->GetExp(eid)->GetInteriorMap(imap);
		      int nbnd   = bmap.size();
		      int nint   = imap.size();
		      int offset = fields[0]->GetCoeff_Offset(eid);
		      for(int j = 0; j < nvel; ++j)
		      {
		           for(int n = 0; n < nz_loc; ++n)
		           {
		                    for(int k = 0; k < nbnd; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
		                    }
		                    
		                    for(int k = 0; k < nint; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
		                    }
		                    cnt  += nbnd;
		                    cnt1 += nint;
		           }
		      }
		}
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
	        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.size()+1);
	        int i;
	        for(i = 0; i < m_fields.size(); ++i)
	        {
	            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
		}
		Eigen::VectorXd diff_x_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		Eigen::VectorXd snap_x = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		Eigen::VectorXd diff_y_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		Eigen::VectorXd snap_y = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		for (int index_recr = 0; index_recr < curr_PhysBaseVec_x.size(); ++index_recr)
		{
			snap_x(index_recr) = snapshot_x_collection[current_index][index_recr];
			snap_y(index_recr) = snapshot_y_collection[current_index][index_recr];
			diff_x_RB_solve(index_recr) = curr_PhysBaseVec_x[index_recr] - snapshot_x_collection[current_index][index_recr];
			diff_y_RB_solve(index_recr) = curr_PhysBaseVec_y[index_recr] - snapshot_y_collection[current_index][index_recr];
		}

		Eigen::MatrixXd curr_xy_reproj = reproject_from_basis(curr_xy_proj);
		
		double relative_L2_error = L2norm_abs_error_ITHACA(curr_PhysBaseVec_x, curr_PhysBaseVec_y, snapshot_x_collection[current_index], snapshot_y_collection[current_index]) / L2norm_ITHACA(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		double relative_Linf_error = Linfnorm_abs_error_ITHACA(curr_PhysBaseVec_x, curr_PhysBaseVec_y, snapshot_x_collection[current_index], snapshot_y_collection[current_index]) / Linfnorm_ITHACA(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		cout << "relative L2 error: " << relative_L2_error << endl;
		cout << "relative Linf error: " << relative_Linf_error << endl;
		cout << "relative euclidean error norm in x coords: " << diff_x_RB_solve.norm() / snap_x.norm() << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean error norm in y coords: " << diff_y_RB_solve.norm() / snap_y.norm() << " of snapshot number " << iter_index << endl;
		Eigen::VectorXd diff_x_proj = curr_xy_reproj.col(0) - snap_x;
		Eigen::VectorXd diff_y_proj = curr_xy_reproj.col(1) - snap_y;
		cout << "relative euclidean projection error norm in x coords: " << diff_x_proj.norm() / snap_x.norm() << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean projection error norm in y coords: " << diff_y_proj.norm() / snap_y.norm() << " of snapshot number " << iter_index << endl;

		cout << "relative euclidean error norm total " << sqrt( diff_x_RB_solve.norm() / snap_x.norm() * diff_x_RB_solve.norm() / snap_x.norm() + diff_y_RB_solve.norm() / snap_y.norm() * diff_y_RB_solve.norm() / snap_y.norm() ) << endl;

	}

	if (compute_smaller_model_errs)
	{
		// repeat the parameter sweep with decreasing RB sizes up to 1, but in a separate function for readability
		for (int i=0; i < RBsize; ++i)
		{
			online_snapshot_check_with_smaller_basis(i);
		}

	}

	
	if (use_fine_grid_VV)
	{
		// repeat the evaluation without the accuracy check
		// fine_general_param_vector is available already
		// start sweeping 

		// Question: how to init?
		// could use all-zero or the cluster-mean
		Array<OneD, NekDouble> cluster_mean_x(snapshot_x_collection[0].size(), 0.0);
		Array<OneD, NekDouble> cluster_mean_y(snapshot_y_collection[0].size(), 0.0);
/*		for (std::set<int>::iterator it=current_cluster.begin(); it!=current_cluster.end(); ++it)
		{
			for (int i = 0; i < snapshot_x_collection[0].size(); ++i)
			{
				cluster_mean_x[i] += (1.0 / current_cluster.size()) * snapshot_x_collection[*it][i];
				cluster_mean_y[i] += (1.0 / current_cluster.size()) * snapshot_y_collection[*it][i];
			}
		} */

		Eigen::VectorXd collected_qoi = Eigen::VectorXd::Zero(Nmax_VV);
		Eigen::VectorXd collected_relative_L2errors = Eigen::VectorXd::Zero(Nmax_VV);
		Eigen::VectorXd collected_relative_Linferrors = Eigen::VectorXd::Zero(Nmax_VV);
		// Eigen::MatrixXd collected_relative_L2errors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		// Eigen::MatrixXd collected_relative_Linferrors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		int fine_grid_dir0_index = 0;
		int fine_grid_dir1_index = 0;
		double locROM_qoi;
		for (int iter_index = 0; iter_index < Nmax_VV; ++iter_index)
		{
			int current_index = iter_index;
			double current_nu;
			double current_geo;
			if (parameter_space_dimension == 1)
			{
				if (parameter_types[0] == 0)
				{
					current_nu = fine_general_param_vector[current_index][0];
				}
				if (parameter_types[0] == 1)
				{
					current_geo = fine_general_param_vector[current_index][0];
					current_nu = m_kinvis;
				}				
			}
			if (debug_mode)
			{
	//			cout << " VV online phase current nu " << current_nu << endl;
	//			cout << " VV online phase current w " << w << endl;
			}
			Set_m_kinvis( current_nu );
			if (use_Newton)
			{
				DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
			}



			Eigen::MatrixXd curr_xy_proj = project_onto_basis(cluster_mean_x, cluster_mean_y);
			Eigen::MatrixXd affine_mat_proj;
			Eigen::VectorXd affine_vec_proj;
			if (parameter_space_dimension == 1)
			{
				if (parameter_types[0] == 0)
				{
					affine_mat_proj = gen_affine_mat_proj(current_nu);
					affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
				}
				if (parameter_types[0] == 1)
				{
					affine_mat_proj = gen_affine_mat_proj_geo(current_nu, current_geo);
					affine_vec_proj = gen_affine_vec_proj_geo(current_nu, current_geo, current_index);
				}
			}
			Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
			double relative_change_error;
			int no_iter=0;
			Array<OneD, double> field_x;
			Array<OneD, double> field_y;
			// now start looping
			do
			{
				// for now only Oseen // otherwise need to do the DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
				Eigen::VectorXd prev_solve_affine = solve_affine;
				Eigen::VectorXd repro_solve_affine = RB * solve_affine;
				Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);

				recover_snapshot_loop(reconstruct_solution, field_x, field_y);
				if (use_Newton)
				{
					DoInitialiseAdv(field_x, field_y);
				}
				curr_xy_proj = project_onto_basis(field_x, field_y);
				if (parameter_space_dimension == 1)
				{
					if (parameter_types[0] == 0)
					{
						affine_mat_proj = gen_affine_mat_proj(current_nu);
						affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
					}
					if (parameter_types[0] == 1)
					{
						affine_mat_proj = gen_affine_mat_proj_geo(current_nu, current_geo);
						affine_vec_proj = gen_affine_vec_proj_geo(current_nu, current_geo, current_index);
					}
				}
				else if (parameter_space_dimension == 2)
				{
					//cout << " VV online phase current nu " << current_nu << endl;
					//cout << " VV online phase current w " << w << endl;
//					affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
//					affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
				}
				solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
				relative_change_error = (solve_affine - prev_solve_affine).norm() / prev_solve_affine.norm();
				//cout << "relative_change_error " << relative_change_error << " no_iter " << no_iter << endl;
				no_iter++;
			} 
			while( ((relative_change_error > 1e-12) && (no_iter < 100)) );
//			cout << "ROM solve no iters used " << no_iter << endl;
			Eigen::VectorXd repro_solve_affine = RB * solve_affine;
			Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
			recover_snapshot_loop(reconstruct_solution, field_x, field_y);
			Eigen::VectorXd  snap_x = Eigen::VectorXd::Zero(cluster_mean_x.size());
			Eigen::VectorXd  snap_y = Eigen::VectorXd::Zero(cluster_mean_x.size());
			for (int index_recr = 0; index_recr < cluster_mean_x.size(); ++index_recr)
			{
				snap_x(index_recr) = field_x[index_recr];
				snap_y(index_recr) = field_y[index_recr];
			}
//			cout << "solved snap_x.norm() " << snap_x.norm() << endl;
//			cout << "solved snap_y.norm() " << snap_y.norm() << endl;

		/*	if (write_ROM_field || (qoi_dof >= 0))
			{
				locROM_qoi = recover_snapshot_data(reconstruct_solution, 0);
			}   */
	//		collected_qoi(iter_index) = locROM_qoi;
			if (use_fine_grid_VV_and_load_ref)
			{
				collected_relative_L2errors(iter_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				collected_relative_Linferrors(iter_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
		/*		if (use_non_unique_up_to_two)
				{
					collected_relative_L2errors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
					collected_relative_Linferrors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
				}    */
			}
/*			fine_grid_dir1_index++;
			if (fine_grid_dir1_index == fine_grid_dir1)
			{
				fine_grid_dir1_index = 0;
				fine_grid_dir0_index++;
			}			*/
		
			cout << "VV iteration " << iter_index << endl;
			cout << "relative L2 error " << collected_relative_L2errors(iter_index) << endl;
			cout << "relative Linf error " << collected_relative_Linferrors(iter_index) << endl;

		} // for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
	}	// if (use_fine_grid_VV)
	
	
    }
    
    
    
    void CoupledLinearNS_ROM::online_snapshot_check_with_smaller_basis(int reduction_int)
	{

	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	cout << "initiate online_snapshot_check_with_smaller_basis RB=" << RBsize - reduction_int <<  endl;
	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);  // is of size M_truth_size
	Eigen::VectorXd collected_relative_euclidean_errors = Eigen::VectorXd::Zero(Nmax);
	Eigen::VectorXd collected_relative_L2_errors = Eigen::VectorXd::Zero(Nmax);
	Eigen::VectorXd collected_relative_Linf_errors = Eigen::VectorXd::Zero(Nmax);
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double w;
		double current_geo;
		if (parameter_space_dimension == 1)
		{
			
			if (parameter_types[0] == 0)
			{
				current_nu = param_vector[current_index];
			}
			if (parameter_types[0] == 1)
			{
				current_geo = param_vector[current_index];
				current_nu = m_kinvis;
			}				
			
		}
		w = current_geo;
		if (debug_mode)
		{
			cout << " online phase current nu " << current_nu << endl;
			cout << " online phase current w " << w << endl;
		}
		Set_m_kinvis( current_nu );
		if (use_Newton)
		{
			DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		}

		Eigen::MatrixXd curr_xy_proj = project_onto_basis(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		Eigen::MatrixXd affine_mat_proj;
		Eigen::VectorXd affine_vec_proj;

		if (parameter_space_dimension == 1)
		{
			if (parameter_types[0] == 0)
			{
				affine_mat_proj = gen_affine_mat_proj(current_nu);
				affine_vec_proj = gen_affine_vec_proj(current_nu, w);
			}
			if (parameter_types[0] == 1)
			{
				affine_mat_proj = gen_affine_mat_proj_geo(current_nu, w);
				affine_vec_proj = gen_affine_vec_proj_geo(current_nu, w, current_index);
			}
		}
		else if (parameter_space_dimension == 2)
		{
			//cout << " VV online phase current nu " << current_nu << endl;
			//cout << " VV online phase current w " << w << endl;
//					affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
//					affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
		}

		Eigen::MatrixXd affine_mat_proj_cut = affine_mat_proj.block(0,0,RBsize-reduction_int,RBsize-reduction_int);
		Eigen::VectorXd affine_vec_proj_cut = affine_vec_proj.head(RBsize-reduction_int);
		Eigen::VectorXd solve_affine_cut = affine_mat_proj_cut.colPivHouseholderQr().solve(affine_vec_proj_cut);
		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
//		cout << "solve_affine " << solve_affine << endl;

//		Eigen::VectorXd repro_solve_affine = RB * solve_affine;
		Eigen::VectorXd repro_solve_affine = RB.leftCols(RBsize-reduction_int) * solve_affine_cut;

		Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
		mat_compare.col(0) = collect_f_all.col(current_index);
		if (debug_mode)
		{
			Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
			current_f_all = collect_f_all.col(current_index);
			Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
			Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
//				cout << "proj_current_f_all_wo_dbc " << proj_current_f_all_wo_dbc << endl;
			Eigen::VectorXd correctRHS = affine_mat_proj * proj_current_f_all_wo_dbc;
			Eigen::VectorXd correction_RHS = correctRHS - affine_vec_proj;
//				cout << "correctRHS " << correctRHS << endl;
//				cout << "correction_RHS " << correction_RHS << endl;
		}
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) - mat_compare.col(0);
//		cout << mat_compare << endl;
//		cout << "relative euclidean error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;

		// compare this to the projection error onto RB, cannot use collect_f_all, is numerically unstable
/*		cout << "mat_compare.cols() " << mat_compare.cols() << endl;
		cout << "mat_compare.rows() " << mat_compare.rows() << endl;
		cout << "collect_f_all.cols() " << collect_f_all.cols() << endl;
		cout << "collect_f_all.rows() " << collect_f_all.rows() << endl;
		cout << "RB.cols() " << RB.cols() << endl;
		cout << "RB.rows() " << RB.rows() << endl; */

		
//		cout << "curr_xy_proj.cols() " << curr_xy_proj.cols() << endl;
	//	cout << "curr_xy_proj.rows() " << curr_xy_proj.rows() << endl; 


//		Eigen::VectorXd proj_solution = collect_f_all.transpose() * mat_compare.col(0);
//		Eigen::VectorXd reproj_solution = collect_f_all * proj_solution;
		Eigen::VectorXd FOM_solution = mat_compare.col(0);
		Eigen::VectorXd FOM_solution_wo_dbc = remove_rows(FOM_solution, elem_loc_dbc);
		Eigen::VectorXd proj_FOM_solution_wo_dbc = RB.transpose() * FOM_solution_wo_dbc;
		Eigen::VectorXd reproj_FOM_solution_wo_dbc = RB * proj_FOM_solution_wo_dbc;
		Eigen::VectorXd reconstruct_FOM_solution = reconstruct_solution_w_dbc(reproj_FOM_solution_wo_dbc);

//		cout << "reconstruct_FOM_solution.norm(): " << reconstruct_FOM_solution.norm() << endl;
//		cout << "reconstruct_solution.norm(): " << reconstruct_solution.norm() << endl;
		Eigen::VectorXd diff_projection = reconstruct_FOM_solution - mat_compare.col(0);

//		cout << "relative euclidean projection error norm: " << diff_projection.norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;


		// now only in RB:
		Eigen::VectorXd diff_projection_RB = reproj_FOM_solution_wo_dbc - remove_rows(mat_compare.col(0), elem_loc_dbc);
		Eigen::VectorXd diff_RB = repro_solve_affine - remove_rows(mat_compare.col(0), elem_loc_dbc);
//		cout << "relative euclidean RB projection error norm: " << diff_projection_RB.norm() / remove_rows(mat_compare.col(0), elem_loc_dbc).norm() << " of snapshot number " << iter_index << endl;
//		cout << "relative euclidean RB error norm: " << diff_RB.norm() / remove_rows(mat_compare.col(0), elem_loc_dbc).norm() << " of snapshot number " << iter_index << endl;

		// have to use curr_xy_proj for better approximations


		if (debug_mode)
		{
			cout << "snapshot_x_collection.size() " << snapshot_x_collection.size() << " snapshot_x_collection[0].size() " << snapshot_x_collection[0].size() << endl;
		}

//		if (write_ROM_field || (qoi_dof >= 0))
//		{
//			recover_snapshot_data(reconstruct_solution, current_index); // this is setting the fields and fieldcoeffs
//		}

	  	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
		Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
		Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(GetNcoeffs());
		Array<OneD, double> field_1(GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
		int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
		int nz_loc = 1;
		int  nplanecoeffs = fields[0]->GetNcoeffs();
		int  nel  = m_fields[0]->GetNumElmts();
		for(int i = 0; i < nel; ++i) 
		{
		      int eid  = i;
		      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
		      fields[0]->GetExp(eid)->GetInteriorMap(imap);
		      int nbnd   = bmap.size();
		      int nint   = imap.size();
			      int offset = fields[0]->GetCoeff_Offset(eid);
		            
		      for(int j = 0; j < nvel; ++j)
		      {
		           for(int n = 0; n < nz_loc; ++n)
		           {
		                    for(int k = 0; k < nbnd; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
		                    }
		                    
		                    for(int k = 0; k < nint; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
		                    }
		                    cnt  += nbnd;
		                    cnt1 += nint;
		           }
		      }
		}
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);

        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.size()+1);
        int i;
        for(i = 0; i < m_fields.size(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
	    }
//		cout << "after recover_snapshot_data have fieldcoeffs[0].size() " << fieldcoeffs[0].size() << endl;

//		cout << "after recover_snapshot_data have curr_PhysBaseVec_x.size() " << curr_PhysBaseVec_x.size() << endl;
//		cout << "after recover_snapshot_data have snapshot_x_collection[current_index].size() " << snapshot_x_collection[current_index].size() << endl;

		// have the FOM_snapshot_solution projection available as curr_xy_proj
		// datastructure: 	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection;
		Eigen::VectorXd diff_x_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		Eigen::VectorXd snap_x = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		Eigen::VectorXd diff_y_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		Eigen::VectorXd snap_y = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.size());
		for (int index_recr = 0; index_recr < curr_PhysBaseVec_x.size(); ++index_recr)
		{
			snap_x(index_recr) = snapshot_x_collection[current_index][index_recr];
			snap_y(index_recr) = snapshot_y_collection[current_index][index_recr];
			diff_x_RB_solve(index_recr) = curr_PhysBaseVec_x[index_recr] - snapshot_x_collection[current_index][index_recr];
			diff_y_RB_solve(index_recr) = curr_PhysBaseVec_y[index_recr] - snapshot_y_collection[current_index][index_recr];
		}

		Eigen::MatrixXd curr_xy_reproj = reproject_from_basis(curr_xy_proj);

		double rel_x = diff_x_RB_solve.norm() / snap_x.norm();
		double rel_y = diff_y_RB_solve.norm() / snap_y.norm();

		cout << "relative euclidean error norm in x coords: " << rel_x << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean error norm in y coords: " << rel_y << " of snapshot number " << iter_index << endl;

		collected_relative_euclidean_errors(iter_index) = sqrt( rel_x * rel_x + rel_y * rel_y );

//		cout << "curr_xy_reproj.cols() " << curr_xy_reproj.cols() << endl;
//		cout << "curr_xy_reproj.rows() " << curr_xy_reproj.rows() << endl;
//		cout << "snap_x.rows() " << snap_x.rows() << endl;

		Eigen::VectorXd diff_x_proj = curr_xy_reproj.col(0) - snap_x;
		Eigen::VectorXd diff_y_proj = curr_xy_reproj.col(1) - snap_y;
		
		cout << "relative euclidean projection error norm in x coords: " << diff_x_proj.norm() / snap_x.norm() << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean projection error norm in y coords: " << diff_y_proj.norm() / snap_y.norm() << " of snapshot number " << iter_index << endl;

		double relative_L2_error = L2norm_abs_error_ITHACA(curr_PhysBaseVec_x, curr_PhysBaseVec_y, snapshot_x_collection[current_index], snapshot_y_collection[current_index]) / L2norm_ITHACA(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		double relative_Linf_error = Linfnorm_abs_error_ITHACA(curr_PhysBaseVec_x, curr_PhysBaseVec_y, snapshot_x_collection[current_index], snapshot_y_collection[current_index]) / Linfnorm_ITHACA(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);

		collected_relative_L2_errors(iter_index) = relative_L2_error;
		collected_relative_Linf_errors(iter_index) = relative_Linf_error;

	} //for (int iter_index = 0; iter_index < Nmax; ++iter_index)

	std::stringstream sstm;
	sstm << "ROM_cluster_reduc" << reduction_int << ".txt";
	std::string ROM_txt = sstm.str();
	const char* outname = ROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			myfile << std::setprecision(17) << collected_relative_euclidean_errors(iter_index) << "\n";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 


	std::stringstream sstmL2;
	sstmL2 << "ROM_cluster_L2reduc" << reduction_int << ".txt";
	std::string ROM_txtL2 = sstmL2.str();
	const char* outnameL2 = ROM_txtL2.c_str();
	ofstream myfileL2 (outnameL2);
	if (myfileL2.is_open())
	{
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			myfileL2 << std::setprecision(17) << collected_relative_L2_errors(iter_index) << "\n";
		}
		myfileL2.close();
	}
	else cout << "Unable to open file"; 



	std::stringstream sstmLinf;
	sstmLinf << "ROM_cluster_Linfreduc" << reduction_int << ".txt";
	std::string ROM_txtLinf = sstmLinf.str();
	const char* outnameLinf = ROM_txtLinf.c_str();
	ofstream myfileLinf (outnameLinf);
	if (myfileLinf.is_open())
	{
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			myfileLinf << std::setprecision(17) << collected_relative_Linf_errors(iter_index) << "\n";
		}
		myfileLinf.close();
	}
	else cout << "Unable to open file"; 




	}


    
    
    Eigen::VectorXd CoupledLinearNS_ROM::gen_affine_vec_proj(double current_nu, int current_index)
    {
	Eigen::VectorXd recovered_affine_adv_rhs_proj_xy = Eigen::VectorXd::Zero(RBsize); 
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_rhs_proj_xy -= adv_vec_proj_x[i] * curr_xy_projected(i,0) + adv_vec_proj_y[i] * curr_xy_projected(i,1);
	}	
	Eigen::VectorXd add_rhs_Newton = Eigen::VectorXd::Zero(RBsize); 
	Eigen::VectorXd recovered_affine_adv_rhs_proj_xy_newton = Eigen::VectorXd::Zero(RBsize);
	Eigen::MatrixXd recovered_affine_adv_rhs_proj_xy_newton_RB = Eigen::MatrixXd::Zero(RBsize,RBsize);  
	if (use_Newton)
	{
		// can I build the Newton-required term from recovered_affine_adv_mat_proj_xy and curr_xy_projected ?
		Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
		for (int i = 0; i < RBsize; ++i)
		{
			recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_xy_projected(i,0) + adv_mats_proj_y[i] * curr_xy_projected(i,1);
			recovered_affine_adv_rhs_proj_xy_newton_RB -= adv_vec_proj_x_newton_RB[i] * curr_xy_projected(i,0) + adv_vec_proj_y_newton_RB[i] * curr_xy_projected(i,1);
		}
		Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
		current_f_all = collect_f_all.col(current_index);
		Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
		Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
		proj_current_f_all_wo_dbc = PODmodes.transpose() * current_f_all;
		if (use_Newton)
		{
			add_rhs_Newton = recovered_affine_adv_rhs_proj_xy_newton_RB.transpose() * proj_current_f_all_wo_dbc;
			add_rhs_Newton = recovered_affine_adv_rhs_proj_xy_newton_RB * proj_current_f_all_wo_dbc;
		}
	}
	return -the_const_one_rhs_proj - current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy  -0.5*add_rhs_Newton ;  
    }

    Eigen::VectorXd CoupledLinearNS_ROM::reconstruct_solution_w_dbc(Eigen::VectorXd reprojected_solve)
    {
	Eigen::VectorXd reconstruct_solution = Eigen::VectorXd::Zero(f_bnd_dbc_full_size.rows());  // is of size M_truth_size
	int counter_wo_dbc = 0;
	for (int row_index=0; row_index < f_bnd_dbc_full_size.rows(); ++row_index)
	{
		if (!elem_loc_dbc.count(row_index))
		{
			reconstruct_solution(row_index) = reprojected_solve(counter_wo_dbc);
			counter_wo_dbc++;
		}
		else
		{
			reconstruct_solution(row_index) = f_bnd_dbc_full_size(row_index);
		}
	}
	return reconstruct_solution;
    }



    Eigen::MatrixXd CoupledLinearNS_ROM::gen_affine_mat_proj(double current_nu)
    {
	Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_xy_projected(i,0) + adv_mats_proj_y[i] * curr_xy_projected(i,1);
	}
	Eigen::MatrixXd affine_mat_proj = the_const_one_proj + current_nu * the_ABCD_one_proj + recovered_affine_adv_mat_proj_xy;
	return affine_mat_proj;
    }

    Eigen::MatrixXd CoupledLinearNS_ROM::project_onto_basis( Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y)
    {
	Eigen::VectorXd c_snapshot_x(GetNpoints());
	Eigen::VectorXd c_snapshot_y(GetNpoints());
	for (int i = 0; i < GetNpoints(); ++i)
	{
		c_snapshot_x(i) = snapshot_x[i];
		c_snapshot_y(i) = snapshot_y[i];
	}
	Eigen::VectorXd curr_x_proj = eigen_phys_basis_x.transpose() * c_snapshot_x;
	Eigen::VectorXd curr_y_proj = eigen_phys_basis_y.transpose() * c_snapshot_y;
	curr_xy_projected = Eigen::MatrixXd::Zero(curr_x_proj.rows(), 2);
	curr_xy_projected.col(0) = curr_x_proj;
	curr_xy_projected.col(1) = curr_y_proj;

	return curr_xy_projected;

    }

    Eigen::MatrixXd CoupledLinearNS_ROM::reproject_from_basis( Eigen::MatrixXd curr_xy_proj )
    {
		Eigen::VectorXd reproj_curr_x = eigen_phys_basis_x * curr_xy_proj.col(0);
		Eigen::VectorXd reproj_curr_y = eigen_phys_basis_y * curr_xy_proj.col(1);
		Eigen::MatrixXd curr_xy_reprojected = Eigen::MatrixXd::Zero(eigen_phys_basis_x.rows(), 2);
		curr_xy_reprojected.col(0) = reproj_curr_x;
		curr_xy_reprojected.col(1) = reproj_curr_y;

		return curr_xy_reprojected;

    }

    /** Virtual function to define if operator in DoSolve is negated
     * with regard to the strong form. This is currently only used in
     * Arnoldi solves. For Coupledd solver this is true since Stokes
     * operator is set up as a LHS rather than RHS operation
     */
    bool CoupledLinearNS_ROM::v_NegatedOp(void)
    {
        return true;
    }

    void CoupledLinearNS_ROM::Solve(void)
    {
        const unsigned int ncmpt = m_velocity.size();
        Array<OneD, Array<OneD, NekDouble> > forcing_phys(ncmpt);
        Array<OneD, Array<OneD, NekDouble> > forcing     (ncmpt);

        for(int i = 0; i < ncmpt; ++i)
        {
            forcing_phys[i] = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNpoints(), 0.0);
            forcing[i]      = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
        }

        for (auto &x : m_forcing)
        {
            const NekDouble time = 0;
            x->Apply(m_fields, forcing_phys, forcing_phys, time);
        }
        for (unsigned int i = 0; i < ncmpt; ++i)
        {
            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[i]->SetWaveSpace(true);
            m_fields[i]->IProductWRTBase(forcing_phys[i], forcing[i]);
            m_fields[i]->SetWaveSpace(waveSpace);
        }

        SolveLinearNS(forcing);
    }

    void CoupledLinearNS_ROM::DefineForcingTerm(void)
    {
        m_ForcingTerm = Array<OneD, Array<OneD, NekDouble> > (m_velocity.size());
        m_ForcingTerm_Coeffs = Array<OneD, Array<OneD, NekDouble> > (m_velocity.size());

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            m_ForcingTerm[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            m_ForcingTerm_Coeffs[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetNcoeffs(),0.0);
        }

        if(m_session->DefinesFunction("ForcingTerm"))
        {
            std::vector<std::string> fieldStr;
            for(int i = 0; i < m_velocity.size(); ++i)
            {
                fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
            }
            GetFunction( "ForcingTerm")->Evaluate(fieldStr,  m_ForcingTerm);
            for(int i = 0; i < m_velocity.size(); ++i)
            {
                m_fields[m_velocity[i]]->FwdTrans_IterPerExp(m_ForcingTerm[i], m_ForcingTerm_Coeffs[i]);
            }
        }
        else
        {
            cout << "'ForcingTerm' section has not been defined in the input file => forcing=0" << endl;
        }
    }

    void CoupledLinearNS_ROM::SolveSteadyNavierStokes(void)
    {
        LibUtilities::Timer Newtontimer;
        Newtontimer.Start();

        Array<OneD, Array<OneD, NekDouble> > RHS_Coeffs(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > RHS_Phys(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > delta_velocity_Phys(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> >Velocity_Phys(m_velocity.size());
        Array<OneD, NekDouble > L2_norm(m_velocity.size(), 1.0);
        Array<OneD, NekDouble > Inf_norm(m_velocity.size(), 1.0);

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            delta_velocity_Phys[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),1.0);
            Velocity_Phys[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
        }

        m_counter=1;

        L2Norm(delta_velocity_Phys, L2_norm);

        //while(max(Inf_norm[0], Inf_norm[1]) > m_tol)
        while(max(L2_norm[0], L2_norm[1]) > m_tol)
        {
            if(m_counter == 1)
                //At the first Newton step, we use the solution of the
                //Stokes problem (if Restart=0 in input file) Or the
                //solution of the .rst file (if Restart=1 in input
                //file)
            {
                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    RHS_Coeffs[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetNcoeffs(),0.0);
                    RHS_Phys[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }

                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    m_fields[m_velocity[i]]->BwdTrans_IterPerExp(m_fields[m_velocity[i]]->GetCoeffs(), Velocity_Phys[i]);
                }

                m_initialStep = true;
                EvaluateNewtonRHS(Velocity_Phys, RHS_Coeffs);
                SetUpCoupledMatrix(0.0, Velocity_Phys, true);
                SolveLinearNS(RHS_Coeffs);
                m_initialStep = false;
            }
            if(m_counter > 1)
            {
                EvaluateNewtonRHS(Velocity_Phys, RHS_Coeffs);

                if(m_counter%m_MatrixSetUpStep == 0) //Setting Up the matrix is expensive. We do it at each "m_MatrixSetUpStep" step.
                {
                    SetUpCoupledMatrix(0.0, Velocity_Phys, true);
                }
                SolveLinearNS(RHS_Coeffs);
            }

            for(int i = 0; i < m_velocity.size(); ++i)
            {
                m_fields[m_velocity[i]]->BwdTrans_IterPerExp(RHS_Coeffs[i], RHS_Phys[i]);
                m_fields[m_velocity[i]]->BwdTrans_IterPerExp(m_fields[m_velocity[i]]->GetCoeffs(), delta_velocity_Phys[i]);
            }

            for(int i = 0; i < m_velocity.size(); ++i)
            {
                Vmath::Vadd(Velocity_Phys[i].size(),Velocity_Phys[i], 1, delta_velocity_Phys[i], 1,
                            Velocity_Phys[i], 1);
            }

            //InfNorm(delta_velocity_Phys, Inf_norm);
            L2Norm(delta_velocity_Phys, L2_norm);

            if(max(Inf_norm[0], Inf_norm[1]) > 100)
            {
                cout<<"\nThe Newton method has failed at m_kinvis = "<<m_kinvis<<" (<=> Re = " << 1/m_kinvis << ")"<< endl;
                ASSERTL0(0, "The Newton method has failed... \n");
            }


            cout << "\n";
            m_counter++;
        }

        if (m_counter > 1) //We save u:=u+\delta u in u->Coeffs
        {
            for(int i = 0; i < m_velocity.size(); ++i)
            {
                m_fields[m_velocity[i]]->FwdTrans(Velocity_Phys[i], m_fields[m_velocity[i]]->UpdateCoeffs());
            }
        }

        Newtontimer.Stop();
        cout<<"We have done "<< m_counter-1 << " iteration(s) in " << Newtontimer.TimePerTest(1)/60 << " minute(s). \n\n";
    }


    void CoupledLinearNS_ROM::Continuation(void)
    {
        Array<OneD, Array<OneD, NekDouble> > u_N(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > tmp_RHS(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > RHS(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > u_star(m_velocity.size());

        cout << "We apply the continuation method: " <<endl;

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            u_N[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            m_fields[m_velocity[i]]->BwdTrans_IterPerExp(m_fields[m_velocity[i]]->GetCoeffs(), u_N[i]);

            RHS[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            tmp_RHS[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);

            m_fields[m_velocity[i]]->PhysDeriv(i, u_N[i], tmp_RHS[i]);
            Vmath::Smul(tmp_RHS[i].size(), m_kinvis, tmp_RHS[i], 1, tmp_RHS[i], 1);

            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[m_velocity[i]]->SetWaveSpace(true);
            m_fields[m_velocity[i]]->IProductWRTDerivBase(i, tmp_RHS[i], RHS[i]);
            m_fields[m_velocity[i]]->SetWaveSpace(waveSpace);
        }

        SetUpCoupledMatrix(0.0, u_N, true);
        SolveLinearNS(RHS);

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            u_star[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            m_fields[m_velocity[i]]->BwdTrans_IterPerExp(m_fields[m_velocity[i]]->GetCoeffs(), u_star[i]);

            //u_star(k+1) = u_N(k) + DeltaKinvis *  u_star(k)
            Vmath::Smul(u_star[i].size(), m_kinvis, u_star[i], 1, u_star[i], 1);
            Vmath::Vadd(u_star[i].size(), u_star[i], 1, u_N[i], 1, u_star[i], 1);

            m_fields[m_velocity[i]]->FwdTrans(u_star[i], m_fields[m_velocity[i]]->UpdateCoeffs());
        }

        m_kinvis -= m_kinvis*m_KinvisPercentage/100;
    }


    void  CoupledLinearNS_ROM::InfNorm(Array<OneD, Array<OneD, NekDouble> > &inarray,
                                   Array<OneD, NekDouble> &outarray)
    {
        for(int i = 0; i < m_velocity.size(); ++i)
        {
            outarray[i] = 0.0;
            for(int j = 0; j < inarray[i].size(); ++j)
            {
                if(inarray[i][j] > outarray[i])
                {
                    outarray[i] = inarray[i][j];
                }
            }
            cout << "InfNorm["<<i<<"] = "<< outarray[i] <<endl;
        }
    }

    void  CoupledLinearNS_ROM::L2Norm(Array<OneD, Array<OneD, NekDouble> > &inarray,
                                  Array<OneD, NekDouble> &outarray)
    {
        for(int i = 0; i < m_velocity.size(); ++i)
        {
            outarray[i] = 0.0;
            for(int j = 0; j < inarray[i].size(); ++j)
            {
                outarray[i] += inarray[i][j]*inarray[i][j];
            }
            outarray[i]=sqrt(outarray[i]);
            cout << "L2Norm["<<i<<"] = "<< outarray[i] <<endl;
        }
    }


    void CoupledLinearNS_ROM::EvaluateNewtonRHS(Array<OneD, Array<OneD, NekDouble> > &Velocity,
                                            Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
        Array<OneD, Array<OneD, NekDouble> > Eval_Adv(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > tmp_DerVel(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > AdvTerm(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > ViscTerm(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > Forc(m_velocity.size());

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            Eval_Adv[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            tmp_DerVel[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);

            AdvTerm[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            ViscTerm[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            Forc[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
            outarray[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);

            m_fields[m_velocity[i]]->PhysDeriv(i, Velocity[i], tmp_DerVel[i]);

            Vmath::Smul(tmp_DerVel[i].size(), m_kinvis, tmp_DerVel[i], 1, tmp_DerVel[i], 1);
        }

        EvaluateAdvectionTerms(Velocity, Eval_Adv, m_time);

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[m_velocity[i]]->SetWaveSpace(true);
            m_fields[m_velocity[i]]->IProductWRTBase(Eval_Adv[i], AdvTerm[i]); //(w, (u.grad)u)
            m_fields[m_velocity[i]]->IProductWRTDerivBase(i, tmp_DerVel[i], ViscTerm[i]); //(grad w, grad u)
            m_fields[m_velocity[i]]->IProductWRTBase(m_ForcingTerm[i], Forc[i]); //(w, f)
            m_fields[m_velocity[i]]->SetWaveSpace(waveSpace);

            Vmath::Vsub(outarray[i].size(), outarray[i], 1, AdvTerm[i], 1, outarray[i], 1);
            Vmath::Vsub(outarray[i].size(), outarray[i], 1, ViscTerm[i], 1, outarray[i], 1);

            Vmath::Vadd(outarray[i].size(), outarray[i], 1, Forc[i], 1, outarray[i], 1);
        }
    }
    
    const SpatialDomains::ExpansionInfoMap
    &CoupledLinearNS_ROM::GenPressureExp(const SpatialDomains::ExpansionInfoMap &VelExp)
    {
        int i;
        SpatialDomains::ExpansionInfoMapShPtr returnval;
        returnval = MemoryManager<SpatialDomains::ExpansionInfoMap>::
            AllocateSharedPtr();

        int nummodes;

        for (auto &expMapIter : VelExp)
        {
            LibUtilities::BasisKeyVector BasisVec;

            for(i = 0; i <  expMapIter.second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey B = expMapIter.second->m_basisKeyVector[i];
                nummodes = B.GetNumModes();
                ASSERTL0(nummodes > 3,"Velocity polynomial space not sufficiently high (>= 4)");
                // Should probably set to be an orthogonal basis.
                LibUtilities::BasisKey newB(B.GetBasisType(),nummodes-2,B.GetPointsKey());
                BasisVec.push_back(newB);
            }
            
            // Put new expansion into list. 
            SpatialDomains::ExpansionInfoShPtr expansionElementShPtr =
                MemoryManager<SpatialDomains::ExpansionInfo>::
                AllocateSharedPtr(expMapIter.second->m_geomShPtr, BasisVec);
            (*returnval)[expMapIter.first] = expansionElementShPtr;
        }
        
        // Save expansion into graph. 
        m_graph->SetExpansionInfo("p",returnval);
        
        return *returnval;
    }

    /**
     *  @param forcing A list of forcing functions for each velocity
     *  component
     *
     *  The routine involves two levels of static
     *  condensations. Initially we require a statically condensed
     *  forcing function which requires the following manipulation
     *
     * \f[ {F\_bnd} = {\bf f}_{bnd} -m\_B \,m\_Cinv\, {\bf f}_{int},
     * \hspace{1cm} F\_p = m\_D\_{int}\, m\_Cinv\, {\bf f}_{int} \f]
     *
     *  Where \f${\bf f}_{bnd}\f$ denote the forcing degrees of
     *  freedom of the elemental velocities on the boundary of the
     *  element, \f${\bf f}_{int}\f$ denote the forcing degrees of
     *  freedom of the elemental velocities on the interior of the
     *  element. (see detailed description for more details).
     *
     * This vector is further manipulated into
     *
     * \f[ Fh\_{bnd} = \left [ \begin{array}{c} f\_{bnd} -m\_B \,
     * m\_Cinv\, {\bf f}_{int}\\ \left [m\_D\_{int} \, m\_Cinv \,{\bf
     * f}_{int} \right]_0 \end{array}\right ] \hspace{1cm} [Fh\_p]_{i} =
     * \begin{array}{c} [m\_D\_{int} \, m\_Cinv \, {\bf
     * f}_{int}]_{i+1} \end{array} \f]
     *
     * where \f$-{[}m\_D\_{int}^T\, m\_Cinv \,{\bf f}_{int}]_0\f$
     * which is corresponds to the mean mode of the pressure degrees
     * of freedom and is now added to the boundary system and the
     * remainder of the block becomes the interior forcing for the
     * inner static condensation (see detailed description for more
     * details) which is set up in a GlobalLinSysDirectStaticCond
     * class.
     *
     * Finally we perform the final maniplation of the forcing to
     * using hte
     * \f[ Fh\_{bnd} = Fh\_{bnd} - m\_Bh \,m\_Chinv \, Fh\_p \f]
     *
     * We can now call the solver to the global coupled boundary
     * system (through the call to #m_CoupledBndSys->Solve) to obtain
     * the velocity boundary solution as the mean pressure solution,
     * i.e.
     *
     * \f[ {\cal A}^T(\hat{A} - \hat{C}^T \hat{D}^{-1} \hat{B} ){\cal
     * A} \, Bnd =  Fh\_{bnd} \f]
     *
     * Once we know the solution to the above the rest of the pressure
     * modes are recoverable thorugh
     *
     *  \f[ Ph = m\_Dhinv\, (Bnd  - m\_Ch^T \, Fh_{bnd}) \f]
     *
     * We can now unpack \f$ Fh\_{bnd} \f$ (last elemental mode) and
     * \f$ Ph \f$ into #m_pressure and \f$ F_p\f$ and \f$ Fh\_{bnd}\f$
     * into a closed pack list of boundary velocoity degrees of
     * freedom stored in \f$ F\_bnd\f$.
     *
     * Finally the interior velocity degrees of freedom are then
     * obtained through the relationship
     *
     *   \f[ F\_{int} = m\_Cinv\ ( F\_{int} + m\_D\_{int}^T\,
     *   F\_p - m\_Btilde^T\, Bnd) \f]
     *
     * We then unpack the solution back to the MultiRegion structures
     * #m_velocity and #m_pressure
     */
    void CoupledLinearNS_ROM::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing)
    {
        int i,n;
        Array<OneD,  MultiRegions::ExpListSharedPtr> vel_fields(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > force(m_velocity.size());

        // Impose Dirichlet conditions on velocity fields
        for(i = 0; i < m_velocity.size(); ++i)
        {
            Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(),1);
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
        }

        if(m_HomogeneousType == eHomogeneous1D)
        {
            int ncoeffsplane = m_fields[m_velocity[0]]->GetPlane(0)->GetNcoeffs();
            for(n = 0; n < m_npointsZ/2; ++n)
            {
                // Get the a Fourier mode of velocity and forcing components.
                for(i = 0; i < m_velocity.size(); ++i)
                {
                    vel_fields[i] = m_fields[m_velocity[i]]->GetPlane(2*n);
                    // Note this needs to correlate with how we pass forcing
                    force[i] = forcing[i] + 2*n*ncoeffsplane;
                }

                SolveLinearNS(force,vel_fields,m_pressure->GetPlane(2*n),n);
            }
            for(i = 0; i < m_velocity.size(); ++i)
            {
                m_fields[m_velocity[i]]->SetPhysState(false);
            }
            m_pressure->SetPhysState(false);
        }
        else
        {
            for(i = 0; i < m_velocity.size(); ++i)
            {
                vel_fields[i] = m_fields[m_velocity[i]];
                // Note this needs to correlate with how we pass forcing
                force[i] = forcing[i];
            }
            SolveLinearNS(force,vel_fields,m_pressure);
        }
    }

    void CoupledLinearNS_ROM::SolveLinearNS(Array<OneD, Array<OneD, NekDouble> > &forcing,
                                        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                        MultiRegions::ExpListSharedPtr &pressure,  const int mode)
    {
        int i,j,k,n,cnt,cnt1;
        int nbnd,nint,offset;
        int nvel = m_velocity.size();
        int nel  = fields[0]->GetNumElmts();
        Array<OneD, unsigned int> bmap, imap;

        Array<OneD, NekDouble > f_bnd(m_mat[mode].m_BCinv->GetRows());
        NekVector< NekDouble  > F_bnd(f_bnd.size(), f_bnd, eWrapper);
        Array<OneD, NekDouble > f_int(m_mat[mode].m_BCinv->GetColumns());
        NekVector< NekDouble  > F_int(f_int.size(),f_int, eWrapper);

        int nz_loc;
        int  nplanecoeffs = fields[m_velocity[0]]->GetNcoeffs();// this is fine since we pass the nplane coeff data.

        if(mode) // Homogeneous mode flag
        {
            nz_loc = 2;
        }
        else
        {
            if(m_singleMode)
            {
                nz_loc = 2;
            }
            else
            {
                nz_loc = 1;
                if(m_HomogeneousType == eHomogeneous1D)
                {
                    Array<OneD, NekDouble> tmp;
                    // Zero fields to set complex mode to zero;
                    for(i = 0; i < fields.size(); ++i)
                    {
                        Vmath::Zero(nplanecoeffs,tmp = fields[i]->UpdateCoeffs()+nplanecoeffs,1);
                    }
                    Vmath::Zero(2*pressure->GetNcoeffs(),pressure->UpdateCoeffs(),1);
                }
            }
        }
        
        for(k = 0; k < nvel; ++k)
        {
            MultiRegions::ContFieldSharedPtr cfield =
                std::dynamic_pointer_cast<MultiRegions::ContField>(fields[k]);

            Array<OneD, NekDouble> sign = cfield->GetLocalToGlobalMap()->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= cfield->GetLocalToGlobalMap()->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            // Add weak boundary conditions to forcing
            const Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                bndConds = fields[k]->GetBndConditions();
            Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp;

            if(m_HomogeneousType == eHomogeneous1D) 
            {
                bndCondExp = m_fields[k]->GetPlane(2*mode)->GetBndCondExpansions();
            }
            else
            {
                bndCondExp = m_fields[k]->GetBndCondExpansions();
            }
            
            for(n = 0; n < nz_loc; ++n)
            {
                int bndcnt = 0;
                for(i = 0; i < bndCondExp.size(); ++i)
                {
                    const Array<OneD, const NekDouble > bndcoeffs =
                        bndCondExp[i]->GetCoeffs();
                    
                    cnt = 0;
                    if(bndConds[i]->GetBoundaryConditionType() ==
                       SpatialDomains::eNeumann ||
                       bndConds[i]->GetBoundaryConditionType() ==
                       SpatialDomains::eRobin)
                    {
                        if(m_locToGloMap[mode]->GetSignChange())
                        {
                            for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                            {
                                forcing[k][n*nplanecoeffs + map[bndcnt+j]] += sign[bndcnt+j] *
                                    bndcoeffs[j]; 
                            }
                        }
                        else
                        {
                            for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                            {
                                forcing[k][n*nplanecoeffs + map[bndcnt+j]] += bndcoeffs[j]; 
                            }
                        }                    
                    }
                    
                    bndcnt += bndCondExp[i]->GetNcoeffs();
                }
            }
        }
        
        Array<OneD, NekDouble > bnd (m_locToGloMap[mode]->GetNumLocalCoeffs(),0.0);

        // Construct f_bnd and f_int and fill in bnd from inarray
        // (with Dirichlet BCs imposed)
        int bndoffset = 0;
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            fields[m_velocity[0]]->GetExp(i)->GetBoundaryMap(bmap);
            fields[m_velocity[0]]->GetExp(i)->GetInteriorMap(imap);
            nbnd   = bmap.size();
            nint   = imap.size();
            offset = fields[m_velocity[0]]->GetCoeff_Offset(i);

            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                Array<OneD, NekDouble> incoeffs = fields[j]->UpdateCoeffs();

                for(n = 0; n < nz_loc; ++n)
                {
                    for(k = 0; k < nbnd; ++k)
                    {
                        f_bnd[cnt+k] = forcing[j][n*nplanecoeffs + 
                                                  offset+bmap[k]];

                        bnd[bndoffset + (n + j*nz_loc)*nbnd + k] =
                            incoeffs[n*nplanecoeffs + offset + bmap[k]];
                    }
                    for(k = 0; k < nint; ++k)
                    {
                        f_int[cnt1+k] = forcing[j][n*nplanecoeffs +
                                                   offset+imap[k]];
                    }

                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
            bndoffset += nvel*nz_loc*nbnd + nz_loc*(pressure->GetExp(i)->GetNcoeffs()); 
        }

        Array<OneD, NekDouble > f_p(m_mat[mode].m_D_int->GetRows());
        NekVector<  NekDouble > F_p(f_p.size(),f_p,eWrapper);
        NekVector<  NekDouble > F_p_tmp(m_mat[mode].m_Cinv->GetRows());

        // fbnd does not currently hold the pressure mean
        F_bnd = F_bnd - (*m_mat[mode].m_BCinv)*F_int;
        F_p_tmp = (*m_mat[mode].m_Cinv)*F_int;
        F_p = (*m_mat[mode].m_D_int) * F_p_tmp;
        
        // construct inner forcing 
        Array<OneD, NekDouble > fh_bnd(m_locToGloMap[mode]->GetNumLocalCoeffs(),0.0);
        
        offset = cnt = 0; 
        for(i = 0; i < nel; ++i)
        {
            nbnd = nz_loc*fields[0]->GetExp(i)->NumBndryCoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    fh_bnd[offset + j*nbnd + k] = 
                        f_bnd[cnt+k];
                }
                cnt += nbnd;
            }
            
            nint    = pressure->GetExp(i)->GetNcoeffs();
            offset += nvel*nbnd + nint*nz_loc; 
        }

        offset = cnt1 = 0;
        for(i = 0; i <  nel; ++i)
        {
            nbnd = nz_loc*fields[0]->GetExp(i)->NumBndryCoeffs(); 
            nint = pressure->GetExp(i)->GetNcoeffs(); 
            
            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    fh_bnd[offset + nvel*nbnd + n*nint+j] = f_p[cnt1+j];
                }
                cnt1   += nint;
            }
            offset += nvel*nbnd + nz_loc*nint;
        }
        m_mat[mode].m_CoupledBndSys->Solve(fh_bnd,bnd,m_locToGloMap[mode]);

        // unpack pressure and velocity boundary systems.
        offset = cnt = 0;
        int totpcoeffs = pressure->GetNcoeffs();
        Array<OneD, NekDouble> p_coeffs = pressure->UpdateCoeffs();
        for(i = 0; i < nel; ++i)
        {
            nbnd = nz_loc*fields[0]->GetExp(i)->NumBndryCoeffs(); 
            nint = pressure->GetExp(i)->GetNcoeffs(); 
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    f_bnd[cnt+k] = bnd[offset + j*nbnd + k];
                }
                cnt += nbnd;
            }
            offset += nvel*nbnd + nint*nz_loc;
        }

        pressure->SetPhysState(false);

        offset = cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i)
        {
            nint = pressure->GetExp(i)->GetNcoeffs(); 
            nbnd = fields[0]->GetExp(i)->NumBndryCoeffs(); 
            cnt1 = pressure->GetCoeff_Offset(i);
            
            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    p_coeffs[n*totpcoeffs + cnt1+j] =
                    f_p[cnt+j] = bnd[offset +
                    nvel*nz_loc*nbnd +
                    n*nint + j];
                }
                cnt += nint;
            }
            offset += (nvel*nbnd + nint)*nz_loc;
        }

        // Back solve first level of static condensation for interior
        // velocity space and store in F_int
        F_int = F_int + Transpose(*m_mat[mode].m_D_int)*F_p
        - Transpose(*m_mat[mode].m_Btilde)*F_bnd;
        F_int = (*m_mat[mode].m_Cinv)*F_int;
        
	curr_f_bnd = Eigen::VectorXd::Zero(f_bnd.size());
	for (int i_phys_dof = 0; i_phys_dof < f_bnd.size(); i_phys_dof++)
	{
		curr_f_bnd(i_phys_dof) = f_bnd[i_phys_dof]; 
	}
	curr_f_p = Eigen::VectorXd::Zero(f_p.size()); 
	for (int i_phys_dof = 0; i_phys_dof < f_p.size(); i_phys_dof++)
	{
		curr_f_p(i_phys_dof) = f_p[i_phys_dof]; 
	}
	curr_f_int = Eigen::VectorXd::Zero(f_int.size()); 
	for (int i_phys_dof = 0; i_phys_dof < f_int.size(); i_phys_dof++)
	{
		curr_f_int(i_phys_dof) = f_int[i_phys_dof]; 
	}        

        // Unpack solution from Bnd and F_int to v_coeffs
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            fields[0]->GetExp(i)->GetBoundaryMap(bmap);
            fields[0]->GetExp(i)->GetInteriorMap(imap);
            nbnd   = bmap.size();
            nint   = imap.size();
            offset = fields[0]->GetCoeff_Offset(i);
            
            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(n = 0; n < nz_loc; ++n)
                {
                    for(k = 0; k < nbnd; ++k)
                    {
                        fields[j]->SetCoeff(n*nplanecoeffs +
                        offset+bmap[k],
                        f_bnd[cnt+k]);
                    }

                    for(k = 0; k < nint; ++k)
                    {
                        fields[j]->SetCoeff(n*nplanecoeffs +
                        offset+imap[k],
                        f_int[cnt1+k]);                       
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }

        for(j = 0; j < nvel; ++j)
        {
            fields[j]->SetPhysState(false);
        }
    }

    void CoupledLinearNS_ROM::v_Output(void)
    {
        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.size()+1);
        std::vector<std::string> variables(m_fields.size()+1);
        int i;

        for(i = 0; i < m_fields.size(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
        }

        fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());
        // project pressure field to velocity space
        if(m_singleMode==true)
        {
            Array<OneD, NekDouble > tmpfieldcoeffs (m_fields[0]->GetNcoeffs()/2);
            m_pressure->GetPlane(0)->BwdTrans_IterPerExp(m_pressure->GetPlane(0)->GetCoeffs(), m_pressure->GetPlane(0)->UpdatePhys());
            m_pressure->GetPlane(1)->BwdTrans_IterPerExp(m_pressure->GetPlane(1)->GetCoeffs(), m_pressure->GetPlane(1)->UpdatePhys());
            m_fields[0]->GetPlane(0)->FwdTrans_IterPerExp(m_pressure->GetPlane(0)->GetPhys(),fieldcoeffs[i]);
            m_fields[0]->GetPlane(1)->FwdTrans_IterPerExp(m_pressure->GetPlane(1)->GetPhys(),tmpfieldcoeffs);
            for(int e=0; e<m_fields[0]->GetNcoeffs()/2; e++)
            {
                fieldcoeffs[i][e+m_fields[0]->GetNcoeffs()/2] = tmpfieldcoeffs[e];
            }
        }
        else
        {
            m_pressure->BwdTrans_IterPerExp(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
            m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);
        }
        variables[i] = "p";

        std::string outname = m_sessionName + ".fld";

        WriteFld(outname,m_fields[0],fieldcoeffs,variables);
    }

    int CoupledLinearNS_ROM::v_GetForceDimension()
    {
        return m_session->GetVariables().size();
    }
    

    Array<OneD, Array<OneD, NekDouble> > CoupledLinearNS_ROM::trafo_current_para(Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y, Array<OneD, NekDouble> parameter_of_interest, Eigen::VectorXd & ref_f_bnd, Eigen::VectorXd & ref_f_p, Eigen::VectorXd & ref_f_int)
    {
//	cout << "starting trafo_current_para " << endl;

	double w = parameter_of_interest[0];	
	double mKinvis = parameter_of_interest[1];
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 

	// verify and transform to bnd / p / int the snapshot data

        int nz_loc;
        nz_loc = 1;

        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nbndry = nsize_bndry;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
	int nint = nsize_int;
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
//	cout << "nsize_int " << nsize_int << endl;
//	cout << "nsize_p " << nsize_p << endl;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;

	Array<OneD, Eigen::MatrixXd > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1
	
	// set the newton forcing in terms of the current iterate
	myAdvField_Newton = Array<OneD, Array<OneD, NekDouble> > (2);
	myAdvField_Newton[0] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	myAdvField_Newton[1] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	for (int i = 0; i<m_fields[0]->GetTotPoints(); ++i)
	{
		myAdvField_Newton[0][i] = snapshot_x[i];
		myAdvField_Newton[1][i] = snapshot_y[i];
	}


	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
		double detT = Geo_T(w, curr_elem_pos, 0);
		double Ta = Geo_T(w, curr_elem_pos, 1);
		double Tb = Geo_T(w, curr_elem_pos, 2);
		double Tc = Geo_T(w, curr_elem_pos, 3);
		double Td = Geo_T(w, curr_elem_pos, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td;
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int npoints = locExp->GetTotPoints();
                int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(curr_elem);
                Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
//                Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
                Array<OneD, NekDouble> tmpphys = Array<OneD, NekDouble>(npoints,0.0);
//		cout << "ncoeffs " << ncoeffs << endl;
//		cout << "nphys " << nphys << endl;
//		cout << "pqsize " << pqsize << endl;   // pqsize == nphys and ncoeffs == nphys / 2 when?

                int nbmap = bmap.size();
                int nimap = imap.size();
		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = snapshot_x[curr_elem*nphys + i];
			curr_snap_y_part[i] = snapshot_y[curr_elem*nphys + i];
		}

                // Calculate derivative of base flow 
                if(use_Newton)
                {
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(int nv = 0; nv < nvel; ++nv)
                    {
                        for(int nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],myAdvField_Newton[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                        }
                    }
                }
                

		Array<OneD, double> Ah_ele_vec(Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec(nsize_int*nsize_int, 0.0);
		Array<OneD, double> Dbnd_ele_vec(nsize_p*nsize_bndry, 0.0);
		Array<OneD, double> Dint_ele_vec(nsize_p*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] );
				}
				for (int j = 0; j < nimap; ++j)
				{
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] );
				}
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
/*					for (int il = 0; il < nphys; ++il)
					{
						cout << "curr_snap_x_part[il] " << curr_snap_x_part[il] << endl;
						cout << "deriv_0[il] " << deriv_0[il] << endl;
						cout << "tmpphys_x[il] " << tmpphys_x[il] << endl;
					} */
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Ta * adv_x_coeffs[int(bmap[j])] + Tc * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Ta * adv_x_coeffs[int(imap[j])] + Tc * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = detT * (Ta * pcoeffs_x[il] + Tc * pcoeffs_y[il]);
					}
				}
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
/*					for (int il = 0; il < nphys; ++il)
					{
						cout << "curr_snap_x_part[il] " << curr_snap_x_part[il] << endl;
						cout << "deriv_0[il] " << deriv_0[il] << endl;
						cout << "tmpphys_x[il] " << tmpphys_x[il] << endl;
					} */
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Tb * adv_x_coeffs[int(bmap[j])] + Td * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Tb * adv_x_coeffs[int(imap[j])] + Td * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = detT * (Tb * pcoeffs_x[il] + Td * pcoeffs_y[il]);
					}


				}


                        if(use_Newton)
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1, AdvDeriv[k*nvel+nv], 1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1) // n1 is only zero
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        Ah_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nbmap)*Ahrows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        C_ele_vec[i+(nv*nz_loc+n1)*nbmap + (j+(k*nz_loc+n1)*nimap)*nbndry] += coeffs[imap[j]];
                                    }
                                }
                            } 
			}

			} //for (int k = 0; k < 2; ++k)


		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					C_ele_vec[ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] );
				}
				for (int j = 0; j < nimap; ++j)
				{
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] );
				}
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Ta * adv_x_coeffs[int(bmap[j])] + Tc * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += detT * (Ta * adv_x_coeffs[int(imap[j])] + Tc * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = detT * (Ta * pcoeffs_x[il] + Tc * pcoeffs_y[il]);
					}
				}
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Tb * adv_x_coeffs[int(bmap[j])] + Td * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += detT * (Tb * adv_x_coeffs[int(imap[j])] + Td * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = detT * (Tb * pcoeffs_x[il] + Td * pcoeffs_y[il]);
					}
				}

                        if(use_Newton)
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1, AdvDeriv[k*nvel+nv], 1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1) // n1 is only zero
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        B_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nimap)*nbndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        D_ele_vec[j+(k*nz_loc+n1)*nimap + (i+(nv*nz_loc+n1)*nimap)*nint] += coeffs[imap[j]];
                                    }
                                }
                            } 
			}



			} //for (int k = 0; k < 2; ++k)


		}

 		// chosen choice: redo the nektar++ approach
		// possible alternatives:
		// build instead the sing_* matrices directly
		// or copy for test purposes from setupcoupledmats??

		Ah_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 
		Dbnd_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_bndry);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				Dbnd_elem[curr_elem](i,j) = Dbnd_ele_vec[ i + j*nsize_p];
			}
		} 
		Dint_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_int);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				Dint_elem[curr_elem](i,j) = Dint_ele_vec[ i + j*nsize_p];
			}
		} 

		
		D_elem[curr_elem] = D_elem[curr_elem].inverse();
		B_elem[curr_elem] = B_elem[curr_elem] * D_elem[curr_elem];
		Ah_elem[curr_elem].block(0,0,nsize_bndry,nsize_bndry) = Ah_elem[curr_elem].block(0,0,nsize_bndry,nsize_bndry) - B_elem[curr_elem] * C_elem[curr_elem].transpose();
		Eigen::MatrixXd Elem_Cinv = D_elem[curr_elem];
		Eigen::MatrixXd Elem_BCinv = B_elem[curr_elem];
		Eigen::MatrixXd Elem_Btilde = C_elem[curr_elem];
		Eigen::MatrixXd Elem_DintCinvDTint = Dint_elem[curr_elem] * Elem_Cinv * Dint_elem[curr_elem].transpose();
		Eigen::MatrixXd Elem_BCinvDTint_m_DTbnd = Elem_BCinv * Dint_elem[curr_elem].transpose() - Dbnd_elem[curr_elem].transpose();
		Eigen::MatrixXd Elem_DintCinvBTtilde_m_Dbnd = Dint_elem[curr_elem] * Elem_Cinv * Elem_Btilde.transpose() - Dbnd_elem[curr_elem];

		Eigen::MatrixXd Bh_curr_ele = Eigen::MatrixXd::Zero( nsize_bndry_p1, nsize_p_m1 );
		Eigen::MatrixXd Ch_curr_ele = Eigen::MatrixXd::Zero( nsize_p_m1, nsize_bndry_p1 );
		Eigen::MatrixXd Dh_curr_ele = Eigen::MatrixXd::Zero( nsize_p_m1, nsize_p_m1 );
		
		for (int i = 0; i < nsize_p_m1; ++i)
		{
			for (int j = 0; j < nsize_p_m1; ++j)
			{
				Dh_curr_ele(i,j) = -Elem_DintCinvDTint(i+1,j+1);
			}
		}

		for (int i = 0; i < nsize_bndry; ++i)
		{
			Ah_elem[curr_elem](i,nsize_bndry) = Elem_BCinvDTint_m_DTbnd(i,0);
			Ah_elem[curr_elem](nsize_bndry,i) = Elem_DintCinvBTtilde_m_Dbnd(0,i);
		}		
		Ah_elem[curr_elem](nsize_bndry,nsize_bndry) = -Elem_DintCinvDTint(0,0);

		for (int j = 0; j < nsize_p_m1; ++j)
		{
			for (int i = 0; i < nsize_bndry; ++i)
			{
				Bh_curr_ele(i,j) = Elem_BCinvDTint_m_DTbnd(i,j+1);
				Ch_curr_ele(j,i) = Elem_DintCinvBTtilde_m_Dbnd(j+1,i);
			}
		}

		for (int j = 0; j < nsize_p_m1; ++j)
		{
			Bh_curr_ele(nsize_bndry, j) = - Elem_DintCinvDTint(0, j+1);
			Ch_curr_ele(j, nsize_bndry) = - Elem_DintCinvDTint(j+1, 0);
		}

		Dh_curr_ele = Dh_curr_ele.inverse();
		Bh_curr_ele = Bh_curr_ele * Dh_curr_ele;
		Ah_elem[curr_elem] = Ah_elem[curr_elem] - Bh_curr_ele * Ch_curr_ele;
		Bh_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_p_m1  );
		Bh_elem[curr_elem] = Bh_curr_ele;
		Ch_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p_m1, nsize_bndry_p1 );
		Ch_elem[curr_elem] = Ch_curr_ele;
		Dh_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p_m1, nsize_p_m1 );
		Dh_elem[curr_elem] = Dh_curr_ele;

		///////////////////////////
		// temporary debugging
//		cout << "Dh_elem[curr_elem](0,0) " << Dh_elem[curr_elem](0,0) << endl;
//		cout << "Dh_elem[curr_elem](1,0) " << Dh_elem[curr_elem](1,0) << endl;
//		cout << "Dh_elem[curr_elem](0,1) " << Dh_elem[curr_elem](0,1) << endl;
//		cout << "Dh_elem[curr_elem](1,1) " << Dh_elem[curr_elem](1,1) << endl;
		/////////////////////////
		


//		sing_A[]



//		cout << "finished curr_elem " << curr_elem << endl;


	} // loop over curr_elem

        
       	int num_elem = m_fields[0]->GetNumElmts();

//        Array<OneD, NekDouble > f_p(num_elem*nsize_p);
//        NekVector<  NekDouble > F_p(f_p.size(),f_p,eWrapper);
//        NekVector<  NekDouble > F_p_tmp(num_elem*nsize_int);
        
         
         

	// nBndDofs already defined
	// NumDirBCs from GetNumGlobalDirBndCoeffs ??

	int NumDirBCs = m_locToGloMap[0]->GetNumGlobalDirBndCoeffs();
	nBndDofs = m_locToGloMap[0]->GetNumGlobalBndCoeffs();  // number of global bnd dofs, also in the .h
	int nGlobHomBndDofs = nBndDofs - NumDirBCs;
	int rows = nGlobHomBndDofs;
	if (debug_mode)
	{
	//	cout << "rows " << rows << endl;
	//	cout << "nBndDofs " << nBndDofs << endl;
	//	cout << "NumDirBCs " << NumDirBCs << endl;
	}
	Eigen::MatrixXd my_Gmat = Eigen::MatrixXd::Zero(rows, rows);

	Eigen::MatrixXd M_trafo = Eigen::MatrixXd::Zero(num_elem*nsize_bndry + num_elem, nGlobHomBndDofs);
	Eigen::MatrixXd M_trafo_no_pp = Eigen::MatrixXd::Zero(num_elem*nsize_bndry, nGlobHomBndDofs);

        const Array<OneD,const int>& loctoglobndmap = m_locToGloMap[0]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign = m_locToGloMap[0]->GetLocalToGlobalBndSign();

	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		int cnt_no_pp = curr_elem*nsize_bndry;
		Eigen::MatrixXd loc_Ah = Ah_elem[curr_elem];
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			int gid1 = loctoglobndmap[cnt + i] - NumDirBCs;
			int sign1 = loctoglobndsign[cnt + i];
			if (gid1 >= 0)
			{
				M_trafo(cnt+i, gid1) = sign1;
				if (i < nsize_bndry)
				{
					M_trafo_no_pp(cnt_no_pp + i, gid1) = sign1;
				}
				for (int j = 0; j < nsize_bndry_p1; ++j)
				{
					int gid2 = loctoglobndmap[cnt + j] - NumDirBCs;
					int sign2 = loctoglobndsign[cnt + j];
					if (gid2 >= 0)
					{
						my_Gmat(gid1,gid2) += sign1*sign2*loc_Ah(i,j);
					}
				}
			}
		}
	}

	int nGlobBndDofs = nBndDofs;

        Array<OneD, NekDouble > f_bnd(num_elem*nsize_bndry);
        NekVector< NekDouble  > F_bnd(f_bnd.size(), f_bnd, eWrapper);
        Array<OneD, NekDouble > f_int(num_elem*nsize_int);
        NekVector< NekDouble  > F_int(f_int.size(),f_int, eWrapper);

	Array<OneD, Array<OneD, NekDouble> > forcing(2); // local dofs
	forcing[0] = Array<OneD, NekDouble>(m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
	forcing[1] = Array<OneD, NekDouble>(m_fields[m_velocity[0]]->GetNcoeffs(),0.0);

	// set the newton forcing in terms of the current iterate
//	myAdvField_Newton = Array<OneD, Array<OneD, NekDouble> > (2);
//	myAdvField_Newton[0] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
//	myAdvField_Newton[1] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
        Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.size());
	Array<OneD, Array<OneD, NekDouble> > Eval_Adv(m_velocity.size());
	Array<OneD, Array<OneD, NekDouble> > AdvTerm(m_velocity.size());
        for(int il = 0; il < m_velocity.size(); ++il)
        {
                AdvField[il] = Array<OneD, NekDouble> (m_fields[m_velocity[il]]->GetTotPoints(),0.0);
		Eval_Adv[il] = Array<OneD, NekDouble> (m_fields[m_velocity[il]]->GetTotPoints(),0.0);
		AdvTerm[il] = Array<OneD, NekDouble> (m_fields[m_velocity[il]]->GetNcoeffs(),0.0);
        }
/*	for (int i = 0; i<m_fields[0]->GetTotPoints(); ++i)
	{
		myAdvField_Newton[0][i] = snapshot_x[i];
		myAdvField_Newton[1][i] = snapshot_y[i];
	}
*/
	EvaluateAdvectionTerms(myAdvField_Newton, Eval_Adv, m_time);
        for (unsigned int i = 0; i < m_velocity.size(); ++i)
        {
            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[i]->SetWaveSpace(true);
//            m_fields[i]->IProductWRTBase(forcing_phys[i], forcing[i]);
            m_fields[m_velocity[i]]->IProductWRTBase(Eval_Adv[i], AdvTerm[i]); //(w, (u.grad)u)
	    if (use_Newton)
	    {
		for (unsigned int il = 0; il < forcing[i].size(); ++il)
	        {
			forcing[i][il] = forcing[i][il] - AdvTerm[i][il];
//			cout << Eval_Adv[i][il] << endl;
		}
	    }
            m_fields[i]->SetWaveSpace(waveSpace);
        }
        
        Array<OneD,  MultiRegions::ExpListSharedPtr> vel_fields(m_velocity.size());
        Array<OneD, Array<OneD, NekDouble> > force(m_velocity.size());

        // Impose Dirichlet conditions on velocity fields
        for(int i = 0; i < m_velocity.size(); ++i)
        {
            Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(),1);
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
        }
        
        // start of adjusted SolveLinearNS function

        int i,j,k,n,cnt,cnt1;
        int nbnd;
        int mode = 0;
      	int nplanecoeffs = m_locToGloMap[0]->GetNumLocalCoeffs();
        Array<OneD, NekDouble > bnd   (m_locToGloMap[0]->GetNumLocalCoeffs(),0.0);
        int offset;
                
        		
	for(k = 0; k < nvel; ++k)
        {
            MultiRegions::ContFieldSharedPtr cfield =
                std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[k]);

            Array<OneD, NekDouble> sign = cfield->GetLocalToGlobalMap()->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= cfield->GetLocalToGlobalMap()->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            // Add weak boundary conditions to forcing
            const Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                bndConds = m_fields[k]->GetBndConditions();
            Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp;

            if(m_HomogeneousType == eHomogeneous1D) 
            {
                bndCondExp = m_fields[k]->GetPlane(2*mode)->GetBndCondExpansions();
            }
            else
            {
                bndCondExp = m_fields[k]->GetBndCondExpansions();
            }
            
            for(n = 0; n < nz_loc; ++n)
            {
                int bndcnt = 0;
                for(i = 0; i < bndCondExp.size(); ++i)
                {
                    const Array<OneD, const NekDouble > bndcoeffs =
                        bndCondExp[i]->GetCoeffs();
                    
                    cnt = 0;
                    if(bndConds[i]->GetBoundaryConditionType() ==
                       SpatialDomains::eNeumann ||
                       bndConds[i]->GetBoundaryConditionType() ==
                       SpatialDomains::eRobin)
                    {
                        if(m_locToGloMap[mode]->GetSignChange())
                        {
                            for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                            {
                                forcing[k][n*nplanecoeffs + map[bndcnt+j]] += sign[bndcnt+j] *
                                    bndcoeffs[j]; 
                            }
                        }
                        else
                        {
                            for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                            {
                                forcing[k][n*nplanecoeffs + map[bndcnt+j]] += bndcoeffs[j]; 
                            }
                        }                    
                    }
                    
                    bndcnt += bndCondExp[i]->GetNcoeffs();
                }
            }
        }
        
        
        // Construct f_bnd and f_int and fill in bnd from inarray
        // (with Dirichlet BCs imposed)
        int bndoffset = 0;
        cnt = cnt1 = 0;

        for(i = 0; i < nel; ++i) // loop over elements
        {
            m_fields[m_velocity[0]]->GetExp(i)->GetBoundaryMap(bmap);
            m_fields[m_velocity[0]]->GetExp(i)->GetInteriorMap(imap);
            nbnd   = bmap.size();
            nint   = imap.size();
            offset = m_fields[m_velocity[0]]->GetCoeff_Offset(i);

	    
            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                Array<OneD, NekDouble> incoeffs = m_fields[j]->UpdateCoeffs();


                for(n = 0; n < nz_loc; ++n)
                {
                    
                    for(k = 0; k < nbnd; ++k)
                    {
                        f_bnd[cnt+k] = forcing[j][n*nplanecoeffs + offset+bmap[k]];
//                        cout << "for incoeffs " << n*nplanecoeffs + offset + bmap[k] << endl;
  //                      cout << "m_fields[j]->GetNcoeffs() " << m_fields[j]->GetNcoeffs() << endl;
  //                      cout << "for bnd " << bndoffset + (n + j*nz_loc)*nbnd + k << endl;
  //                      cout << "m_locToGloMap[0]->GetNumLocalCoeffs() " << m_locToGloMap[0]->GetNumLocalCoeffs() << endl;
                        bnd[bndoffset + (n + j*nz_loc)*nbnd + k] = incoeffs[n*nplanecoeffs + offset + bmap[k]];
                    }
                    
                    for(k = 0; k < nint; ++k)
                    {
                        f_int[cnt1+k] = forcing[j][n*nplanecoeffs + offset+imap[k]];
                    }
                    

                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
            
            bndoffset += nvel*nz_loc*nbnd + nz_loc*(m_pressure->GetExp(i)->GetNcoeffs()); 
        }
        

//        Array<OneD, NekDouble > f_p(m_mat[mode].m_D_int->GetRows());
	int f_p_size = num_elem*nsize_p;
	//cout << "f_p_size " << f_p_size << endl;
        Array<OneD, NekDouble > f_p(f_p_size);
        // cout << "find malloc 2" << endl;
        NekVector<  NekDouble > F_p(f_p.size(),f_p,eWrapper);
        // cout << "find malloc 3" << endl;
        NekVector<  NekDouble > F_p_tmp(m_mat[0].m_Cinv->GetRows());
        // cout << "find malloc 4" << endl;


	Eigen::VectorXd f_bnd_rhs_eigen = Eigen::VectorXd::Zero(f_bnd.size());
	for (int i = 0; i < f_bnd.size(); ++i)
	{
		f_bnd_rhs_eigen(i) = f_bnd[i];
	}
	Eigen::VectorXd f_p_rhs_eigen = Eigen::VectorXd::Zero(f_p.size());
	for (int i = 0; i < f_p.size(); ++i)
	{
	//	f_p_rhs_eigen(i) = f_p[i]; // should be undefined at this stage
	}
	Eigen::VectorXd f_int_rhs_eigen = Eigen::VectorXd::Zero(f_int.size()); // num_elem*nsize_int
	Eigen::VectorXd f_p_tmp_eigen = Eigen::VectorXd::Zero(f_int.size()); // num_elem*nsize_int
	for (int i = 0; i < f_int.size(); ++i)
	{
		f_int_rhs_eigen(i) = f_int[i];
	}



        // fbnd does not currently hold the pressure mean
//        F_bnd = F_bnd - (*m_mat[mode].m_BCinv)*F_int;
//        F_p_tmp = (*m_mat[mode].m_Cinv)*F_int;
//        F_p = (*m_mat[mode].m_D_int) * F_p_tmp;

//	 should do these operations in Eigen
	for (int i = 0; i < num_elem; ++i)
	{
		f_bnd_rhs_eigen.segment(i*nsize_bndry, nsize_bndry) = f_bnd_rhs_eigen.segment(i*nsize_bndry, nsize_bndry) - B_elem[i] * f_int_rhs_eigen.segment(i*nsize_int, nsize_int);
		f_p_tmp_eigen.segment(i*nsize_int, nsize_int) = D_elem[i] * f_int_rhs_eigen.segment(i*nsize_int, nsize_int);
		f_p_rhs_eigen.segment(i*nsize_p, nsize_p) = Dint_elem[i] * f_p_tmp_eigen.segment(i*nsize_int, nsize_int);

	}


	////////////////////
	// temporary debugging
/*	cout << "f_bnd_rhs_eigen.size() " << f_bnd_rhs_eigen.size() << endl;
	cout << "f_bnd_rhs_eigen.norm() " << f_bnd_rhs_eigen.norm() << endl;
	cout << "f_p_rhs_eigen.size() " << f_p_rhs_eigen.size() << endl;
	cout << "f_p_rhs_eigen.norm() " << f_p_rhs_eigen.norm() << endl;

*/
	///////////////////


	// construct inner forcing
        Array<OneD, NekDouble > fh_bnd(m_locToGloMap[0]->GetNumLocalCoeffs(),0.0);

        const Array<OneD,const int>& loctoglomap = m_locToGloMap[0]->GetLocalToGlobalMap();
        const Array<OneD,const NekDouble>& loctoglosign = m_locToGloMap[0]->GetLocalToGlobalSign();
        
	offset = 0;
	cnt = 0; 
        for(int i = 0; i < num_elem; ++i)
        {
            int eid  = i;
            int nbnd = nz_loc*m_fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            
            for(int j = 0; j < nvel; ++j)
            {
                for(int k = 0; k < nbnd; ++k)
                {
// old                    fh_bnd[loctoglomap[offset+j*nbnd+k]] += loctoglosign[offset+j*nbnd+k]*f_bnd_rhs_eigen(cnt+k);
                    fh_bnd[offset+j*nbnd+k] = f_bnd_rhs_eigen(cnt+k);
                }
                cnt += nbnd;
            }
            
            int nint    = m_pressure->GetExp(eid)->GetNcoeffs();
            offset += nvel*nbnd + nint*nz_loc; 
        }
        
        offset = cnt1 = 0; 
        for(int i = 0; i <  num_elem; ++i)
        {
            int eid  = i;
            int nbnd = nz_loc*m_fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            int nint = m_pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(int n = 0; n < nz_loc; ++n)
            {
                for(int j = 0; j < nint; ++j)
                {
                    fh_bnd[offset + nvel*nbnd + n*nint+j] = f_p_rhs_eigen(cnt1+j);
                }
                cnt1   += nint;
            }
            offset += nvel*nbnd + nz_loc*nint; 
        }





            bool atLastLevel       = m_locToGloMap[mode]->AtLastLevel();
            int  scLevel           = m_locToGloMap[mode]->GetStaticCondLevel();
            
            int nGlobDofs    = m_locToGloMap[mode]->GetNumGlobalCoeffs();
            int nLocBndDofs  = m_locToGloMap[mode]->GetNumLocalBndCoeffs();
        //    int nGlobBndDofs = m_locToGloMap[mode]->GetNumGlobalBndCoeffs();
            int nDirBndDofs  = m_locToGloMap[mode]->GetNumGlobalDirBndCoeffs();
            int nIntDofs     = nGlobDofs - nGlobBndDofs;
            int nLocalBnd    = m_locToGloMap[mode]->GetNumLocalBndCoeffs();

	
	Array<OneD, NekDouble> m_wsp = Array<OneD, NekDouble> (3*nLocalBnd + nIntDofs, 0.0);
	
	
            Array<OneD, NekDouble> F_bnd_glssc, F_bnd1_glssc, F_int_glssc, V_bnd_glssc; 
            Array<OneD, NekDouble> tmp_glssc;

            F_bnd_glssc  = m_wsp;
            F_bnd1_glssc = m_wsp + nLocBndDofs;
            V_bnd_glssc  = m_wsp + 2*nLocBndDofs;
            F_int_glssc  = m_wsp + 3*nLocBndDofs; 
            
            m_locToGloMap[mode]->LocalToLocalBnd(bnd,V_bnd_glssc);

            NekVector<NekDouble> F_Int_glssc(nIntDofs, F_int_glssc,eWrapper);
            NekVector<NekDouble> V_Bnd_glssc(nLocBndDofs,V_bnd_glssc,eWrapper);            
            
            if(nIntDofs)
            {
                m_locToGloMap[mode]->LocalToLocalInt(fh_bnd, F_int_glssc);
            }

	Eigen::VectorXd F_Bnd_glssc_eigen = Eigen::VectorXd::Zero(F_bnd_glssc.size());

            // Boundary system solution
            if(nGlobBndDofs-nDirBndDofs)
            {
                m_locToGloMap[mode]->LocalToLocalBnd(fh_bnd,F_bnd_glssc);                
//                v_PreSolve(scLevel, F_bnd_glssc);  // set up normalisation factor for right hand side on first SC level                
                NekVector<NekDouble> F_Bnd_glssc(nLocBndDofs,F_bnd1_glssc,eWrapper); // Gather boundary expansison into locbnd                 





                if(nIntDofs)   // construct boundary forcing
                {
           //         DNekScalBlkMat &BinvD      = *m_BinvD;
           //         F_Bnd_glssc = BinvD*F_Int_glssc; 
                    

	Eigen::VectorXd F_Int_glssc_eigen = Eigen::VectorXd::Zero(F_int_glssc.size());
	for (int F_Int_glssc_eigen_index = 0; F_Int_glssc_eigen_index < F_int_glssc.size(); ++F_Int_glssc_eigen_index)
	{
		F_Int_glssc_eigen(F_Int_glssc_eigen_index) = F_Int_glssc[F_Int_glssc_eigen_index];
	}                    
	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		Eigen::MatrixXd loc_Ah = Ah_elem[curr_elem];
		Eigen::MatrixXd loc_Bh = Bh_elem[curr_elem]; // of size (nsize_bndry_p1, nsize_p_m1) 
		Eigen::VectorXd loc_F_Int_glssc_eigen = F_Int_glssc_eigen.segment(curr_elem*nsize_p_m1,nsize_p_m1);
		F_Bnd_glssc_eigen.segment(curr_elem*nsize_bndry_p1,nsize_bndry_p1) = loc_Bh * loc_F_Int_glssc_eigen;
	//	cout << "curr_elem*nsize_p_m1 " << curr_elem*nsize_p_m1 << endl;
	//	cout << "curr_elem*nsize_bndry_p1 " << curr_elem*nsize_bndry_p1 << endl;
	}
	for (int F_Bnd_glssc_eigen_index = 0; F_Bnd_glssc_eigen_index < nLocBndDofs; ++F_Bnd_glssc_eigen_index)
	{
		F_Bnd_glssc[F_Bnd_glssc_eigen_index] = F_Bnd_glssc_eigen(F_Bnd_glssc_eigen_index);
	}                    
	
	                    
                    Vmath::Vsub(nLocBndDofs, F_bnd_glssc,1, F_bnd1_glssc,1, F_bnd_glssc,1);
               }
               
//		cout << "value of atLastLevel " << atLastLevel << endl;               	
               	
//                if(atLastLevel)
                if(1)
                {                    
			//                    v_BasisFwdTransform(F_bnd_glssc);  // Transform to new basis if required 
			//                    DNekScalBlkMat &SchurCompl = *m_schurCompl;
			//                    v_CoeffsFwdTransform(V_bnd_glssc,V_bnd_glssc);                    
		  	// need to do this in local coordinates,so use the Ah                 F_Bnd_glssc = my_Gmat*V_Bnd_glssc;   // subtract dirichlet boundary forcing
		  	Eigen::VectorXd V_Bnd_glssc_eigen = Eigen::VectorXd::Zero(nLocBndDofs);
			for (int V_Bnd_glssc_eigen_index = 0; V_Bnd_glssc_eigen_index < nLocBndDofs; ++V_Bnd_glssc_eigen_index)
			{
				V_Bnd_glssc_eigen(V_Bnd_glssc_eigen_index) = V_Bnd_glssc[V_Bnd_glssc_eigen_index];
			}  
			for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
			{
				int cnt = curr_elem*nsize_bndry_p1;
				Eigen::MatrixXd loc_Ah = Ah_elem[curr_elem]; // of size (nsize_bndry_p1, nsize_bndry_p1)
				Eigen::MatrixXd loc_Bh = Bh_elem[curr_elem]; // of size (nsize_bndry_p1, nsize_p_m1) 
				Eigen::VectorXd loc_V_Bnd_glssc_eigen = V_Bnd_glssc_eigen.segment(curr_elem*nsize_bndry_p1,nsize_bndry_p1);
				F_Bnd_glssc_eigen.segment(curr_elem*nsize_bndry_p1,nsize_bndry_p1) = loc_Ah * loc_V_Bnd_glssc_eigen;
			//	cout << "curr_elem*nsize_p_m1 " << curr_elem*nsize_p_m1 << endl;
			//	cout << "curr_elem*nsize_bndry_p1 " << curr_elem*nsize_bndry_p1 << endl;
			}
			for (int F_Bnd_glssc_eigen_index = 0; F_Bnd_glssc_eigen_index < nLocBndDofs; ++F_Bnd_glssc_eigen_index)
			{
				F_Bnd_glssc[F_Bnd_glssc_eigen_index] = F_Bnd_glssc_eigen(F_Bnd_glssc_eigen_index);
			}	
	                Vmath::Vsub(nLocBndDofs, F_bnd_glssc,1, F_bnd1_glssc, 1, F_bnd_glssc,1);
        	        Array<OneD, NekDouble> F_hom, pert(nGlobBndDofs,0.0);
                	m_locToGloMap[mode]->AssembleBnd(F_bnd_glssc, F_bnd1_glssc);                    
              		//      SolveLinearSystem(nGlobBndDofs, F_bnd1_glssc, pert, pLocToGloMap, nDirBndDofs);   // Solve for difference from initial solution given inout;
			int nHomDofs = nGlobBndDofs - nDirBndDofs ;
			Eigen::VectorXd my_invec = Eigen::VectorXd::Zero(nHomDofs);
			for (int my_invec_index = 0; my_invec_index < nHomDofs; ++my_invec_index)
			{
				my_invec(my_invec_index) = F_bnd1_glssc[nDirBndDofs + my_invec_index]; 
			}
			/////////////////////// actual solve here ////////////////////////////////
			Eigen::VectorXd my_Asolution = my_Gmat.colPivHouseholderQr().solve(my_invec);
			//////////////////////////////////////////////////////////////////////////
			// write my_Asolution back into pert
			for (int my_Asolution_index = 0; my_Asolution_index < nHomDofs; ++my_Asolution_index)
			{
				pert[nDirBndDofs + my_Asolution_index] = my_Asolution(my_Asolution_index); 
			}
	                Array<OneD, NekDouble> outloc = F_bnd_glssc; 
        	        m_locToGloMap[mode]->GlobalToLocalBnd(pert,outloc);                    
                	Vmath::Vadd(nLocBndDofs, V_bnd_glssc, 1, outloc, 1, V_bnd_glssc,1);   // Add back initial conditions onto difference                    
			//    v_CoeffsBwdTransform(V_bnd_glssc);    // Transform back to original basis                    
                	m_locToGloMap[mode]->LocalBndToLocal(V_bnd_glssc,bnd);    // put final bnd solution back in output array
                }
                else // Process next level of recursion for multi level SC
                {
                	cout << "Process next level of recursion for multi level SC undefined" << endl;
                	Nektar::MultiRegions::AssemblyMapSharedPtr nextLevLocToGloMap = m_locToGloMap[mode]->GetNextLevelLocalToGlobalMap();
                    // partially assemble F for next level and
                    // reshuffle V_bnd
             	        nextLevLocToGloMap->PatchAssemble     (F_bnd_glssc,F_bnd_glssc);
                	nextLevLocToGloMap->PatchLocalToGlobal(V_bnd_glssc,V_bnd_glssc);
               //         m_recursiveSchurCompl->Solve(F_bnd_glssc,V_bnd_glssc, nextLevLocToGloMap);
                    // unpack V_bnd
                //    nextLevLocToGloMap->PatchGlobalToLocal(V_bnd_glssc,V_bnd_glssc);
                    // place V_bnd in output array
                //    m_locToGloMap[mode]->LocalBndToLocal(V_bnd_glssc, bnd);
                }
            }



            // solve interior system
            if(nIntDofs)
            {
                Array<OneD, NekDouble> V_int(nIntDofs);
                NekVector<NekDouble>   V_Int(nIntDofs, V_int ,eWrapper);

                // get array of local solutions
//                DNekScalBlkMat &invD  = *m_invD;
//                DNekScalBlkMat &C     = *m_C;
//                F_Int = F_Int - C*V_Bnd;
//                Multiply(V_Int, invD, F_Int); // is doing V_Int = invD * F_int 

//	change F_int_glssc and V_bnd_glssc into Eigen datastructures
	Eigen::VectorXd F_Int_glssc_eigen = Eigen::VectorXd::Zero(F_int_glssc.size());
	for (int F_Int_glssc_eigen_index = 0; F_Int_glssc_eigen_index < F_int_glssc.size(); ++F_Int_glssc_eigen_index)
	{
		F_Int_glssc_eigen(F_Int_glssc_eigen_index) = F_Int_glssc[F_Int_glssc_eigen_index];
	}
	Eigen::VectorXd V_Bnd_glssc_eigen = Eigen::VectorXd::Zero(nLocBndDofs);
	for (int V_Bnd_glssc_eigen_index = 0; V_Bnd_glssc_eigen_index < nLocBndDofs; ++V_Bnd_glssc_eigen_index)
	{
		V_Bnd_glssc_eigen(V_Bnd_glssc_eigen_index) = V_Bnd_glssc[V_Bnd_glssc_eigen_index];
	}  
	Eigen::VectorXd V_Int_glssc_eigen = Eigen::VectorXd::Zero(nIntDofs);
	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		Eigen::MatrixXd loc_Ch = Ch_elem[curr_elem];
//		F_int[curr_elem*nsize_p_m1:curr_elem*nsize_p_m1+nsize_p_m1] = fh_bnd[nGlobHomBndDofs + NumDirBCs + curr_elem*nsize_p_m1: nGlobHomBndDofs + NumDirBCs + curr_elem*nsize_p_m1+nsize_p_m1] - np.dot(loc_Ch,  loc_dbc[cnt:cnt+nsize_bndry_p1])
		F_Int_glssc_eigen.segment(curr_elem*nsize_p_m1,nsize_p_m1) = F_Int_glssc_eigen.segment(curr_elem*nsize_p_m1,nsize_p_m1).eval() - loc_Ch * V_Bnd_glssc_eigen.segment(cnt, nsize_bndry_p1);  // actually an fh_bnd in case of body forcing as well
		Eigen::MatrixXd loc_Dh = Dh_elem[curr_elem];
		V_Int_glssc_eigen.segment(curr_elem*nsize_p_m1,nsize_p_m1) = loc_Dh * F_Int_glssc_eigen.segment(curr_elem*nsize_p_m1,nsize_p_m1).eval();
	}
	for (int V_Int_glssc_eigen_index = 0; V_Int_glssc_eigen_index < nIntDofs; ++V_Int_glssc_eigen_index)
	{
		V_int[V_Int_glssc_eigen_index] = V_Int_glssc_eigen(V_Int_glssc_eigen_index);
	}  


	
                m_locToGloMap[mode]->LocalIntToLocal(V_int, bnd);
            }


	
	
	

	
//	int nLocBndDofs = num_elem * nsize_bndry_p1;
	// ??	nGlobDofs = ??
	// nGlobBndDofs defined??, si, si, si come nBndDofs
//	int nDirBndDofs = NumDirBCs;
/*	Eigen::VectorXd loc_dbc = Eigen::VectorXd::Zero(nLocBndDofs);
	Eigen::VectorXd loc_dbc_pt1 = Eigen::VectorXd::Zero(nLocBndDofs);
	Eigen::VectorXd loc_dbc_pt2 = Eigen::VectorXd::Zero(nLocBndDofs);
	Eigen::VectorXd V_GlobBnd = Eigen::VectorXd::Zero(nGlobBndDofs);
*/
/*	cout << "bnd.size() is the same as m_locToGloMap[0]->GetNumGlobalCoeffs() " << bnd.size() << endl;
	cout << "nGlobBndDofs " << nGlobBndDofs << endl;
	cout << "m_locToGloMap[0]->AtLastLevel() " << m_locToGloMap[0]->AtLastLevel() << endl;	
	cout << "m_locToGloMap[0]->GetStaticCondLevel() " << m_locToGloMap[0]->GetStaticCondLevel() << endl;
	cout << "m_locToGloMap[0]->GetLowestStaticCondLevel() " << m_locToGloMap[0]->GetLowestStaticCondLevel() << endl;
*/


/*	for (int i = 0; i < nGlobBndDofs; ++i)
	{
		V_GlobBnd(i) = bnd[i];
	}
	
	for (int i = 0; i < nLocBndDofs; ++i)
	{
		loc_dbc(i) = loctoglobndsign[i] * V_GlobBnd(loctoglobndmap[i]);
	}
*/

	// can I now copy-paste from Nektar ?? -- try to do so:

//        const Array<OneD,const int>& loctoglomap = m_locToGloMap[0]->GetLocalToGlobalMap();
//        const Array<OneD,const NekDouble>& loctoglosign = m_locToGloMap[0]->GetLocalToGlobalSign();


        // unpack pressure and velocity boundary systems. 
        offset = cnt = 0;
        int totpcoeffs = m_pressure->GetNcoeffs();
        Array<OneD, NekDouble> p_coeffs = m_pressure->UpdateCoeffs();
        for(i = 0; i < nel; ++i)
        {
            nbnd = nz_loc*m_fields[0]->GetExp(i)->NumBndryCoeffs(); 
            nint = m_pressure->GetExp(i)->GetNcoeffs(); 
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    f_bnd[cnt+k] = bnd[offset + j*nbnd + k];
                }
                cnt += nbnd;
            }
            offset += nvel*nbnd + nint*nz_loc;
        }

        m_pressure->SetPhysState(false);
        
  //      Array<OneD, NekDouble > f_p(num_elem*nsize_p);
//        NekVector<  NekDouble > F_p(f_p.size(),f_p,eWrapper);

        offset = cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i)
        {
            nint = m_pressure->GetExp(i)->GetNcoeffs(); 
            nbnd = m_fields[0]->GetExp(i)->NumBndryCoeffs(); 
            cnt1 = m_pressure->GetCoeff_Offset(i);
            
            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    p_coeffs[n*totpcoeffs + cnt1+j] =
                    f_p[cnt+j] = bnd[offset +
                    nvel*nz_loc*nbnd +
                    n*nint + j];
                }
                cnt += nint;
            }
            offset += (nvel*nbnd + nint)*nz_loc;
        }

        // Back solve first level of static condensation for interior
        // velocity space and store in F_int
	// F_int should initially be empty without body forcing

//        F_int = F_int + Transpose(*m_mat[mode].m_D_int)*F_p - Transpose(*m_mat[mode].m_Btilde)*F_bnd;
//        F_int = (*m_mat[mode].m_Cinv)*F_int;
		// like to do the MatVec in Eigen and the transform into Nektar?? for f_p and bnd

	Eigen::VectorXd f_bnd_eigen = Eigen::VectorXd::Zero(f_bnd.size());
	for (int i = 0; i < f_bnd.size(); ++i)
	{
		f_bnd_eigen(i) = f_bnd[i];
	}
	Eigen::VectorXd f_p_eigen = Eigen::VectorXd::Zero(f_p.size());
	for (int i = 0; i < f_p.size(); ++i)
	{
		f_p_eigen(i) = f_p[i];
	}
	Eigen::VectorXd f_int_eigen = Eigen::VectorXd::Zero(f_int.size()); // num_elem*nsize_int

	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt_Dint = curr_elem*nsize_p;
		int cnt_C = curr_elem*nsize_bndry;
		int cnt_D = curr_elem*nsize_int;
		Eigen::MatrixXd loc_Dint = Dint_elem[curr_elem];
		Eigen::MatrixXd loc_C = C_elem[curr_elem];
		Eigen::MatrixXd loc_D = D_elem[curr_elem];
		// f_int_rhs should come in as well...
//		cout << "f_p_eigen.rows, f_p_eigen.cols " << f_p_eigen.rows() << " " << f_p_eigen.cols() << endl;
//		cout << "f_bnd_eigen.rows, f_bnd_eigen.cols " << f_bnd_eigen.rows() << " " << f_bnd_eigen.cols() << endl;
//		Eigen::VectorXd nn =  f_p_eigen.segment(cnt_Dint, nsize_p);
//		cout << "nn.rows, nn.cols " << nn.rows() << " " << nn.cols() << endl;
		f_int_eigen.segment(curr_elem*nsize_int, nsize_int) = f_int_rhs_eigen.segment(curr_elem*nsize_int, nsize_int) + loc_Dint.transpose() * f_p_eigen.segment(cnt_Dint, nsize_p) - loc_C.transpose() * f_bnd_eigen.segment(cnt_C, nsize_bndry);  // also f_int_rhs here
		f_int_eigen.segment(curr_elem*nsize_int, nsize_int) = loc_D * f_int_eigen.segment(curr_elem*nsize_int, nsize_int);
	}

	for (int i = 0; i < f_int.size(); ++i)
	{
		f_int[i] = f_int_eigen(i);
	}

	int nlc = m_locToGloMap[0]->GetNumLocalCoeffs();
	nplanecoeffs = nlc;

        // Unpack solution from Bnd and F_int to v_coeffs 
        cnt = 0; 
	cnt1 = 0;
        for(int i = 0; i < nel; ++i) // loop over elements
        {
            int eid  = i;
            m_fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
            m_fields[0]->GetExp(eid)->GetInteriorMap(imap);
            int nbnd   = bmap.size();
            int nint   = imap.size();
            int offset = m_fields[0]->GetCoeff_Offset(eid);
            
            for(int j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(int n = 0; n < nz_loc; ++n)
                {
                    for(int k = 0; k < nbnd; ++k)
                    {
                        m_fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd[cnt+k]);
                    }
                    
                    for(int k = 0; k < nint; ++k)
                    {
                        m_fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int[cnt1+k]);
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }
        
        for(int j = 0; j < nvel; ++j) 
        {
            m_fields[j]->SetPhysState(false);
        }

	// the fields should be the same as the input:  Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = UpdateFields();
	m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys());
	m_fields_t[1]->BwdTrans(m_fields_t[1]->GetCoeffs(), m_fields_t[1]->UpdatePhys());
	Array<OneD, NekDouble> out_field_trafo_x(GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_trafo_y(GetNpoints(), 0.0);

	Eigen::VectorXd csx0_trafo(GetNpoints());
	Eigen::VectorXd csy0_trafo(GetNpoints());
	Eigen::VectorXd csx0(GetNpoints());
	Eigen::VectorXd csy0(GetNpoints());

	CopyFromPhysField(0, out_field_trafo_x); 
	CopyFromPhysField(1, out_field_trafo_y);
	for( int index_conv = 0; index_conv < GetNpoints(); ++index_conv)
	{
		csx0_trafo(index_conv) = out_field_trafo_x[index_conv];
		csy0_trafo(index_conv) = out_field_trafo_y[index_conv];
		csx0(index_conv) = snapshot_x[index_conv];
		csy0(index_conv) = snapshot_y[index_conv];
	}

	if (debug_mode)
	{
//		cout << "csx0.norm() " << csx0.norm() << endl;
//		cout << "csx0_trafo.norm() " << csx0_trafo.norm() << endl;
//		cout << "csy0.norm() " << csy0.norm() << endl;
//		cout << "csy0_trafo.norm() " << csy0_trafo.norm() << endl;
	}

	Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y(2);
	snapshot_result_phys_velocity_x_y[0] = Array<OneD, NekDouble>(GetNpoints());
	snapshot_result_phys_velocity_x_y[1] = Array<OneD, NekDouble>(GetNpoints());
	snapshot_result_phys_velocity_x_y[0] = out_field_trafo_x;
	snapshot_result_phys_velocity_x_y[1] = out_field_trafo_y;

	curr_f_bnd = Eigen::VectorXd::Zero(f_bnd.size());
	for (int i_phys_dof = 0; i_phys_dof < f_bnd.size(); i_phys_dof++)
	{
		curr_f_bnd(i_phys_dof) = f_bnd[i_phys_dof]; 
	}
	curr_f_p = Eigen::VectorXd::Zero(f_p.size()); 
	for (int i_phys_dof = 0; i_phys_dof < f_p.size(); i_phys_dof++)
	{
		curr_f_p(i_phys_dof) = f_p[i_phys_dof]; 
	}
	curr_f_int = Eigen::VectorXd::Zero(f_int.size()); 
	for (int i_phys_dof = 0; i_phys_dof < f_int.size(); i_phys_dof++)
	{
		curr_f_int(i_phys_dof) = f_int[i_phys_dof]; 
	}
	
	ref_f_bnd = curr_f_bnd;
	ref_f_p = curr_f_p;
	ref_f_int = curr_f_int;

	return snapshot_result_phys_velocity_x_y;



/*

		for i in range(0, nimap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(imap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			deriv_0 = np.dot(phys, (loc_cm0))
			deriv_1 = np.dot(phys, (loc_cm1))
			coeffs_0_0 = np.dot(deriv_0, loc_IP_d0)
			coeffs_0_1 = np.dot(deriv_0, loc_IP_d1)
			coeffs_1_0 = np.dot(deriv_1, loc_IP_d0)
			coeffs_1_1 = np.dot(deriv_1, loc_IP_d1)
			for k in range(0, 2):
				for j in range(0, nbmap):
					C_ele_vec[j+k*nbmap + (i+k*nimap)*nsize_bndry] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] )
				for j in range(0, nimap):
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] )
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					deriv_y = np.dot(phys, (loc_cm1))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					tmpphys_y = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv_y
					coeffs = np.dot(tmpphys, loc_IP)
					coeffs_y = np.dot(tmpphys_y, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Ta * coeffs[int(bmap[j])] + Tc * coeffs_y[int(bmap[j])])
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += detT * (Ta * coeffs[int(imap[j])] + Tc * coeffs_y[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp)
					pcoeff_y = np.dot(deriv_y, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = detT * (Ta * pcoeff + Tc * pcoeff_y)
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Tb * coeffs[int(bmap[j])] + Td * coeffs[int(bmap[j])])
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += detT * (Tb * coeffs[int(imap[j])] + Td * coeffs[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = detT * (Tb * pcoeff + Td * pcoeff)

		for i in range(0, nbmap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(bmap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			deriv_0 = np.dot(phys, (loc_cm0))
			deriv_1 = np.dot(phys, (loc_cm1))
			coeffs_0_0 = np.dot(deriv_0, loc_IP_d0)
			coeffs_0_1 = np.dot(deriv_0, loc_IP_d1)
			coeffs_1_0 = np.dot(deriv_1, loc_IP_d0)
			coeffs_1_1 = np.dot(deriv_1, loc_IP_d1)
			for k in range(0, 2):
				for j in range(0, nbmap):
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] )
				for j in range(0, nimap):
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] )
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					deriv_y = np.dot(phys, (loc_cm1))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					tmpphys_y = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv_y
					coeffs = np.dot(tmpphys, loc_IP)
					coeffs_y = np.dot(tmpphys_y, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Ta * coeffs[int(bmap[j])] + Tc * coeffs_y[int(bmap[j])])
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Ta * coeffs[int(imap[j])] + Tc * coeffs_y[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp)
					pcoeff_y = np.dot(deriv_y, loc_IPp)
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = detT * (Ta * pcoeff + Tc * pcoeff_y)
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Tb * coeffs[int(bmap[j])] + Td * coeffs[int(bmap[j])])
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Tb * coeffs[int(imap[j])] + Td * coeffs[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp) 
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = detT * (Tb * pcoeff + Td * pcoeff)											


*/

/*	for curr_elem in range(0, num_elem):
		Ah_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		H1_bnd_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		B_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		C_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		D_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		DnoK_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		Dadv_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		H1_int_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		D1_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		Dbnd_ele_vec = 0*np.arange(nsize_bndry*nsize_p*1.0)
		Dint_ele_vec = 0*np.arange(nsize_int*nsize_p*1.0)
		loc_bwd_mat = bwdtrans[(curr_elem*ncoeffs) : (curr_elem*ncoeffs + ncoeffs) ,:]
		loc_cm0 = cartmap0[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_cm1 = cartmap1[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP = IP[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP_d0 = IP_d0[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP_d1 = IP_d1[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IPp = IPp[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		HelmMat = HelmMats[curr_elem,:]
		L00 = Lapl00[curr_elem,:]
		L11 = Lapl11[curr_elem,:]
		mimic_L00 = 0*Lapl00[curr_elem,:]
		curr_elem_pos = get_curr_elem_pos(curr_elem, geo_trafo_list)
		detT = Geo_T(w, curr_elem_pos, 0)
		Ta = Geo_T(w, curr_elem_pos, 1)
		Tb = Geo_T(w, curr_elem_pos, 2)
		Tc = Geo_T(w, curr_elem_pos, 3)
		Td = Geo_T(w, curr_elem_pos, 4)
		c00 = Ta*Ta + Tb*Tb
		c01 = Ta*Tc + Tb*Td
		c11 = Tc*Tc + Td*Td			
*/
    }    // closing bracket for trafo_current_para
 
    int CoupledLinearNS_ROM::get_curr_elem_pos(int curr_elem)
    {
	// find within Array<OneD, std::set<int> > elements_trafo; the current entry
	for(int i = 0; i < elements_trafo.size(); ++i)
	{
		if (elements_trafo[i].count(curr_elem))
		{
			return i;
		}
	}
	return -1;
    }     
    
  /*  double CoupledLinearNS_ROM::Geo_T(double w, int elemT, int index)
    {
	Eigen::Matrix2d T;
	if (elemT == 0) 
	{
		T << 1, 0, 0, 1/w;
	}
	else if (elemT == 1) 
	{
		T << 1, 0, 0, 2/(3-w);
	}
	else if (elemT == 2) 
	{
		T << 1, 0, -(1-w)/2, 1;
	}
	else if (elemT == 3) 
	{
		T << 1, 0, -(w-1)/2, 1;
	}
	else if (elemT == 4) 
	{
		T << 1, 0, 0, 1;
	}
	if (index == 0) 
	{
		return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
	}
	else if (index == 1) 
	{
		return T(0,0);
	}
	else if (index == 2) 
	{
		return T(0,1);
	}
	else if (index == 3) 
	{
		return T(1,0);
	}
	else if (index == 4) 
	{
		return T(1,1);
	}
	return 0;
    }
*/

    void CoupledLinearNS_ROM::write_curr_field(std::string filename)
    {

        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.size()+1);
        std::vector<std::string> variables(m_fields.size()+1);
        int i;
        
        for(i = 0; i < m_fields.size(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
	    if (debug_mode)
	    {
//		    cout << "variables[i] " << variables[i] << endl;
	    }
        }

	if (debug_mode)
	{
//		cout << "m_singleMode " << m_singleMode << endl;	
	}

	fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs(), 0.0);  
        variables[i] = "p"; 

	WriteFld(filename,m_fields[0],fieldcoeffs,variables);

    }

    Eigen::MatrixXd CoupledLinearNS_ROM::DoTrafo_1p_geo()
    {
 
 
     	cout << " transform snapshots to different format for single geometric parameter ... " << endl;
    	cout << " need to set general_param_vector " << endl;
     	
     	// general_param_vector expects first the geo parameter and then the kinvis
     	// here the kinvis is constant and given by m_kinvis
     	
     	Array<OneD, NekDouble> current_parameter_vector = Array<OneD, NekDouble>(2);
     	current_parameter_vector[0] = 0;
     	current_parameter_vector[1] = m_kinvis;
     	
        for(int i = 0; i < Nmax; ++i)
	{
		//cout << "param_vector[i] " << param_vector[i] << endl;
	}     	
     	
 	Eigen::MatrixXd transformed_snapshots;
	// setting transformed_snapshots, making use of snapshot_x_collection, snapshot_y_collection

//	cout << "checking if all quantities set: " << endl;
//	cout << "Nmax " << Nmax << endl;
//	cout << "curr_f_bnd.size() " << curr_f_bnd.size() << endl;
//	cout << "curr_f_p.size() " << curr_f_p.size() << endl;
//	cout << "curr_f_int.size() " << curr_f_int.size() << endl;

	Eigen::MatrixXd collect_f_bnd( curr_f_bnd.size() , Nmax );
	Eigen::MatrixXd collect_f_p(   curr_f_p.size()   , Nmax );
	Eigen::MatrixXd collect_f_int( curr_f_int.size() , Nmax );

	Array<OneD, NekDouble> collected_qoi = Array<OneD, NekDouble> (Nmax);

        for(int i = 0; i < Nmax; ++i)
	{
		current_parameter_vector[0] = param_vector[i];
		
//		cout << "general_param_vector[i][0] " << general_param_vector[i][0] << endl;
//		cout << "general_param_vector[i][1] " << general_param_vector[i][1] << endl;
//		cout << "curr_f_bnd.size() " << curr_f_bnd.size() << endl;
//		cout << "curr_f_p.size() " << curr_f_p.size() << endl;
//		cout << "curr_f_int.size() " << curr_f_int.size() << endl;
		// identify the geometry one - introduce in the .h some vars keeping the para_type
		// here now [0] is geometry 'w' and [1] is k_invis

		//for (int repeat_i = 0; repeat_i < 10; ++repeat_i)

		// take a timing for trafo_current_para

		Eigen::VectorXd ref_f_bnd;
		Eigen::VectorXd ref_f_p;
		Eigen::VectorXd ref_f_int;

		// if using Newton should set myAdvField

		Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_x_collection[i], snapshot_y_collection[i], current_parameter_vector, ref_f_bnd, ref_f_p, ref_f_int); 

//		cout << "curr_f_bnd.norm() " << curr_f_bnd.norm() << endl;
//		cout << "ref_f_bnd.norm() " << ref_f_bnd.norm() << endl;

//		Eigen::VectorXd trafo_f_bnd = curr_f_bnd;
//		Eigen::VectorXd trafo_f_p = curr_f_p;
//		Eigen::VectorXd trafo_f_int = curr_f_int;

		if (do_trafo_check)
		{
			double L2error = 1;
			do
			{
				Array<OneD, Array<OneD, NekDouble> > prev_snapshot_result_phys_velocity_x_y = snapshot_result_phys_velocity_x_y;
				snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_result_phys_velocity_x_y[0], snapshot_result_phys_velocity_x_y[1], current_parameter_vector, ref_f_bnd, ref_f_p, ref_f_int);
//				cout << "curr_f_bnd.norm() " << curr_f_bnd.norm() << endl;
//				cout << "ref_f_bnd.norm() " << ref_f_bnd.norm() << endl;

//				trafo_f_bnd = curr_f_bnd;
//				trafo_f_p = curr_f_p;
//				trafo_f_int = curr_f_int;

				double L2error_x = L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]);
				double L2error_y = L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]);
				double L2error_x_ref = L2Error(0);
				double L2error_y_ref = L2Error(1);
				L2error = sqrt(L2error_x*L2error_x + L2error_y*L2error_y) / sqrt(L2error_x_ref*L2error_x_ref + L2error_y_ref*L2error_y_ref);
				cout << "relative L2error w.r.t. current iterate " << L2error << endl;
//				cout << " L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]) " << L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]) << endl;
//				cout << " L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]) " << L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]) << endl;
			}
			while ((L2error > do_trafo_check_relative_error) && (!load_cO_snapshot_data_from_files));
		}

//		snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_result_phys_velocity_x_y[0], snapshot_result_phys_velocity_x_y[1], general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int);

//		cout << "curr_f_bnd.norm() " << curr_f_bnd.norm() << endl;
//		cout << "ref_f_bnd.norm() " << ref_f_bnd.norm() << endl;

		if (qoi_dof >= 0)
		{
//			cout << "converged qoi dof " << snapshot_result_phys_velocity_x_y[1][qoi_dof] << endl;
			collected_qoi[i] = snapshot_result_phys_velocity_x_y[1][qoi_dof];
		}

		Eigen::VectorXd trafo_f_bnd = curr_f_bnd;
		Eigen::VectorXd trafo_f_p = curr_f_p;
		Eigen::VectorXd trafo_f_int = curr_f_int;

		collect_f_bnd.col(i) = trafo_f_bnd;
		collect_f_p.col(i) = trafo_f_p;
		collect_f_int.col(i) = trafo_f_int;

		// need to replace the snapshot data with the converged one for error computations
		// that means replace data in the snapshot_x_collection and snapshot_y_collection
		if (replace_snapshot_with_transformed)
		{
			snapshot_x_collection[i] = snapshot_result_phys_velocity_x_y[0];
			snapshot_y_collection[i] = snapshot_result_phys_velocity_x_y[1];
		}
		
//		cout << "collect_f_bnd.col(i).norm() " << collect_f_bnd.col(i).norm() << endl;

		// generate the correct string
		std::stringstream sstm;
		sstm << "Conv_Oseen_param" << i << ".fld";
		std::string filename = sstm.str();
		if (!load_cO_snapshot_data_from_files)
		{
			write_curr_field(filename);
		}

	}

	std::stringstream sstm;
	sstm << "FOM_qoi.txt";
	std::string LocROM_txt = sstm.str();
	const char* outname = LocROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < Nmax; i0++)
		{
			myfile << collected_qoi[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 	


	transformed_snapshots = Eigen::MatrixXd::Zero( curr_f_bnd.size()+curr_f_p.size()+curr_f_int.size() , Nmax );
	transformed_snapshots.block(0,0,collect_f_bnd.rows(),collect_f_bnd.cols()) = collect_f_bnd;
	transformed_snapshots.block(collect_f_bnd.rows(),0,collect_f_p.rows(),collect_f_p.cols()) = collect_f_p;
	transformed_snapshots.block(collect_f_bnd.rows()+collect_f_p.rows(),0,collect_f_int.rows(),collect_f_int.cols()) = collect_f_int;

	// do the same for VV reference solutions, if that is being used
	if (use_fine_grid_VV)
	{
    		for(int i = 0; i < fine_grid_dir0; ++i)
		{
			current_parameter_vector[0] = fine_general_param_vector[i][0];
			cout << "\n attempting trafo for VV reference solutions \n \n";
			double err_threshold = 1e-9;
			int num_iter = 0;
			Eigen::VectorXd ref_f_bnd;
			Eigen::VectorXd ref_f_p;
			Eigen::VectorXd ref_f_int;
			Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_x_collection_VV[i], snapshot_y_collection_VV[i], current_parameter_vector, ref_f_bnd, ref_f_p, ref_f_int);  
			if (1)
			{
				double L2error = 1;
				do
				{
					Array<OneD, Array<OneD, NekDouble> > prev_snapshot_result_phys_velocity_x_y = snapshot_result_phys_velocity_x_y;
					snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_result_phys_velocity_x_y[0], snapshot_result_phys_velocity_x_y[1], current_parameter_vector, ref_f_bnd, ref_f_p, ref_f_int); 
					double L2error_x = L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]);
					double L2error_y = L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]);
					double L2error_x_ref = L2Error(0);
					double L2error_y_ref = L2Error(1);
					L2error = sqrt(L2error_x*L2error_x + L2error_y*L2error_y) / sqrt(L2error_x_ref*L2error_x_ref + L2error_y_ref*L2error_y_ref);
					cout << "relative L2error w.r.t. current iterate " << L2error << endl;
					num_iter++;
					if (num_iter == 100)
					{
						num_iter = 0;
						err_threshold = err_threshold*10;
					}	
				}
				while ((L2error > err_threshold));
			}

		std::stringstream sstm;
		sstm << "Conv_ref_Oseen_param" << i << ".fld";
		std::string filename = sstm.str();
		if (1)
		{
			write_curr_field(filename);
		}

	}
	}
	
     	cout << " ... finished transform snapshots " << endl;

	return transformed_snapshots;
 

    }


    Array<OneD, Array<OneD, Eigen::MatrixXd > > CoupledLinearNS_ROM::gen_adv_mats_proj_x_geo(Array<OneD, double> curr_PhysBaseVec_x, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_x_2d)
    {
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  curr_adv_mats_proj_x_2d(number_elem_trafo);
	adv_vec_proj_x_2d = Array<OneD, Array<OneD, Eigen::VectorXd > >(number_elem_trafo);
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		curr_adv_mats_proj_x_2d[i] = Array<OneD, Eigen::MatrixXd > (2);
		curr_adv_mats_proj_x_2d[i][0] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		curr_adv_mats_proj_x_2d[i][1] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		adv_vec_proj_x_2d[i] = Array<OneD, Eigen::VectorXd > (2);
		adv_vec_proj_x_2d[i][0] = Eigen::VectorXd::Zero(RBsize);
		adv_vec_proj_x_2d[i][1] = Eigen::VectorXd::Zero(RBsize);
	}
	// not what I wanted:   cout << "sizeof curr_adv_mats_proj_x_2d " << sizeof(curr_adv_mats_proj_x_2d) << endl;
	
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc;
        nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 

	Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
//	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
//	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
//	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.size();
                int nimap = imap.size();
		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = curr_PhysBaseVec_x[curr_elem*nphys + i];
		}
		Array<OneD, Array<OneD, double> > Ah_ele_vec(2);
		Ah_ele_vec[0] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[1] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, Array<OneD, double> > B_ele_vec(2);
		B_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > C_ele_vec(2); 
		C_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > D_ele_vec(2);
		D_ele_vec[0] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		D_ele_vec[1] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[0][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_x_coeffs[int(bmap[j])];
							Ah_ele_vec[1][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[0][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(imap[j])];
							C_ele_vec[1][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(imap[j])];
						}
					}

				}
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[0][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(bmap[j])];
							B_ele_vec[1][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[0][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_x_coeffs[int(imap[j])];
							D_ele_vec[1][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_y_coeffs[int(imap[j])];
						}
					}
				}
			}

		} // for (int i = 0; i < nimap; ++i)

		Ah_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		Ah_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		Ah_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem][0](i,j) = Ah_ele_vec[0][ i + j*Ahrows ];
				Ah_elem[curr_elem][1](i,j) = Ah_ele_vec[1][ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		B_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		B_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem][0](i,j) = B_ele_vec[0][ i + j*nsize_bndry ];
				B_elem[curr_elem][1](i,j) = B_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		C_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		C_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem][0](i,j) = C_ele_vec[0][ i + j*nsize_bndry ];
				C_elem[curr_elem][1](i,j) = C_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		D_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		D_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem][0](i,j) = D_ele_vec[0][ i + j*nsize_int];
				D_elem[curr_elem][1](i,j) = D_ele_vec[1][ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	// would need a function which (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects
	for (int i = 0; i < number_elem_trafo; ++i)
	{

		time_t timer_1;
		time_t timer_2;
		time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

		curr_adv_mats_proj_x_2d[i][0] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 0, adv_vec_proj_x_2d[i][0]);

		time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
		double seconds = difftime(timer_1, timer_2);
//		cout << "time for a single adv_geo_mat_projector in seconds " << seconds << endl;


		curr_adv_mats_proj_x_2d[i][1] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 1, adv_vec_proj_x_2d[i][1]);
	}


	return curr_adv_mats_proj_x_2d;

    }
    
    Array<OneD, Array<OneD, Eigen::MatrixXd > > CoupledLinearNS_ROM::gen_adv_mats_proj_y_geo(Array<OneD, double> curr_PhysBaseVec_y, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_y_2d)
    {
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  curr_adv_mats_proj_y_2d(number_elem_trafo);
	adv_vec_proj_y_2d = Array<OneD, Array<OneD, Eigen::VectorXd > >(number_elem_trafo);
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		curr_adv_mats_proj_y_2d[i] = Array<OneD, Eigen::MatrixXd > (2);
		curr_adv_mats_proj_y_2d[i][0] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		curr_adv_mats_proj_y_2d[i][1] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		adv_vec_proj_y_2d[i] = Array<OneD, Eigen::VectorXd > (2);
		adv_vec_proj_y_2d[i][0] = Eigen::VectorXd::Zero(RBsize);
		adv_vec_proj_y_2d[i][1] = Eigen::VectorXd::Zero(RBsize);
	}
	// not what I wanted:   cout << "sizeof curr_adv_mats_proj_x_2d " << sizeof(curr_adv_mats_proj_x_2d) << endl;
	
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc;
        nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 

	Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
//	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
//	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
//	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.size();
                int nimap = imap.size();
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_y_part[i] = curr_PhysBaseVec_y[curr_elem*nphys + i];
		}
		Array<OneD, Array<OneD, double> > Ah_ele_vec(2);
		Ah_ele_vec[0] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[1] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, Array<OneD, double> > B_ele_vec(2);
		B_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > C_ele_vec(2); 
		C_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > D_ele_vec(2);
		D_ele_vec[0] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		D_ele_vec[1] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[0][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_x_coeffs[int(bmap[j])];
							Ah_ele_vec[1][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[0][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(imap[j])];
							C_ele_vec[1][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(imap[j])];
						}
					}

				}
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[0][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(bmap[j])];
							B_ele_vec[1][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[0][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_x_coeffs[int(imap[j])];
							D_ele_vec[1][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_y_coeffs[int(imap[j])];
						}
					}
				}
			}

		} // for (int i = 0; i < nimap; ++i)

		Ah_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		Ah_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		Ah_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem][0](i,j) = Ah_ele_vec[0][ i + j*Ahrows ];
				Ah_elem[curr_elem][1](i,j) = Ah_ele_vec[1][ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		B_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		B_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem][0](i,j) = B_ele_vec[0][ i + j*nsize_bndry ];
				B_elem[curr_elem][1](i,j) = B_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		C_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		C_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem][0](i,j) = C_ele_vec[0][ i + j*nsize_bndry ];
				C_elem[curr_elem][1](i,j) = C_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		D_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		D_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem][0](i,j) = D_ele_vec[0][ i + j*nsize_int];
				D_elem[curr_elem][1](i,j) = D_ele_vec[1][ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	// would need a function which (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		curr_adv_mats_proj_y_2d[i][0] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 0, adv_vec_proj_y_2d[i][0]);
		curr_adv_mats_proj_y_2d[i][1] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 1, adv_vec_proj_y_2d[i][1]);
	}


	return curr_adv_mats_proj_y_2d;

    }    
 
    Eigen::MatrixXd CoupledLinearNS_ROM::adv_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &adv_vec_proj)
    {

	// eats up all the compute time...

	// this function (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
        int nz_loc;
        nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int*nel, nsize_int*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel, nsize_int*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel, nsize_int*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel, nsize_bndry*nel );

	time_t timer_1;
	time_t timer_2;
	time_t timer_3;
	time_t timer_4;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
//		cout << " i " << i << endl;
		int curr_elem_pos = get_curr_elem_pos(i);
		if (curr_elem_pos == curr_elem_trafo)
		{
			Eigen::MatrixXd curr_A_elem = Ah_elem[i][deriv_index].block(0, 0, Ah_elem[i][deriv_index].rows() - 1, Ah_elem[i][deriv_index].cols() - 1);
			Eigen::MatrixXd curr_B_elem = B_elem[i][deriv_index];
			Eigen::MatrixXd curr_C_elem = C_elem[i][deriv_index];
			Eigen::MatrixXd curr_D_elem = D_elem[i][deriv_index];

			A_adv_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = curr_A_elem;
			B_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_B_elem;
			C_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_C_elem;
			D_adv_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = curr_D_elem;
	
		}
	}

	time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
	double seconds = difftime(timer_2, timer_1);
//	cout << "time within adv_geo_mat_projector for block write in seconds " << seconds << endl;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

/*	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_adv_all * Mtrafo;
			adv_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_adv_all;
			adv_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_adv_all.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 2:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
	}
*/
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;


	time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_2, timer_1);
//	cout << "time within adv_geo_mat_projector for adv_matrix write in seconds " << seconds << endl;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

	// also need to create appropriate vector right-hand-side contribution
	Eigen::VectorXd add_to_rhs_adv(M_truth_size); // probably need this for adv and non-adv
	add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;

	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
	Eigen::VectorXd adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc); // or is this eating up the time?
	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for remove_rows(add_to_rhs_adv, elem_loc_dbc) in seconds " << seconds << endl;

	adv_vec_proj = RB.transpose() * adv_rhs_add;

	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
	Eigen::MatrixXd adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);
	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for remove_cols_and_rows in seconds " << seconds << endl;

	// it is limiting the Eigen MatVec operation -- could compare to the Nektar++ MatVec operation -- I rather think one has to employ the sparsity structure of A

	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
	Eigen::MatrixXd adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB; // is this eating up the time?
	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for RB.transpose() * adv_matrix_simplified * RB in seconds " << seconds << endl;

//	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
//	Eigen::MatrixXd nn = (adv_matrix_simplified * RB.col(2)); // is this eating up the time?
//	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
//	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for (adv_matrix_simplified * RB.col(2)) " << seconds << endl;

//	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
//	Eigen::MatrixXd nn2 = RB.transpose(); // is this eating up the time?
//	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
//	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for RB.transpose() in seconds " << seconds << endl;

//	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
//	Eigen::MatrixXd nn3 = nn2 * nn; // is this eating up the time?
//	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
//	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for nn2 * nn in seconds " << seconds << endl; 

	time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_2, timer_1);
//	cout << "time within adv_geo_mat_projector for some MatVecs in seconds " << seconds << endl;

//	cout << " finished " << endl;

	return adv_mat_proj;

    }
 
    void CoupledLinearNS_ROM::gen_proj_adv_terms_geo()
    {
	RBsize = RB.cols();

	adv_mats_proj_x = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_mats_proj_y = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_x = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_y = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_x_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_y_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);

	// to be superseded by:
	adv_mats_proj_x_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > (RBsize); // should be RBsize x number_elem_trafo x 2 x RBsize x RBsize
	adv_mats_proj_y_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > (RBsize);
	adv_vec_proj_x_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > (RBsize);
	adv_vec_proj_y_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > (RBsize);

	Array<OneD, double> PhysBase_zero(GetNpoints(), 0.0);
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		time_t timer_1;
		time_t timer_2;
		time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];

		adv_mats_proj_x_2d[trafo_iter] = gen_adv_mats_proj_x_geo(curr_PhysBaseVec_x, adv_vec_proj_x_2d[trafo_iter]);

		time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
		double seconds = difftime(timer_2, timer_1);
		if (debug_mode)
		{
			cout << "time for a single gen_adv_mats_proj_x_2d in seconds " << seconds << endl;
			cout << "have 2 times RBsize of that " << endl;
		}

		adv_mats_proj_y_2d[trafo_iter] = gen_adv_mats_proj_y_geo(curr_PhysBaseVec_y, adv_vec_proj_y_2d[trafo_iter]);
	
	}
	

//	Array<OneD, double> PhysBase_zero(GetNpoints(), 0.0);
/*
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];

		InitObject();

		DoInitialiseAdv(curr_PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
		// needs to be replaced with a more gen. term		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(RB_A.rows() + RB_Dbnd.rows() + RB_C.cols(), RB_A.cols() + RB_Dbnd.rows() + RB_B.cols() );
		Eigen::MatrixXd adv_matrix;
		adv_matrix = Get_advection_matrix();
		Eigen::VectorXd add_to_rhs_adv(M_truth_size); // probably need this for adv and non-adv
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;
		Eigen::MatrixXd adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);

		adv_vec_proj_x_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);

		if (use_Newton)
		{
			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
				adv_vec_proj_x_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}


		Eigen::VectorXd adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		Eigen::MatrixXd adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		Eigen::VectorXd adv_rhs_proj = RB.transpose() * adv_rhs_add;

		adv_mats_proj_x[trafo_iter] = adv_mat_proj;
		adv_vec_proj_x[trafo_iter] = adv_rhs_proj;

		adv_vec_proj_y_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);
		DoInitialiseAdv(PhysBase_zero , curr_PhysBaseVec_y ); // call with parameter in phys state
		adv_matrix = Get_advection_matrix();
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   
		adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);
		adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		adv_rhs_proj = RB.transpose() * adv_rhs_add;
		adv_mats_proj_y[trafo_iter] = adv_mat_proj;
		adv_vec_proj_y[trafo_iter] = adv_rhs_proj;

		if (use_Newton)
		{
			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
				adv_vec_proj_y_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}
	} */
    }

    void CoupledLinearNS_ROM::gen_reference_matrices_geo()
    {
	// should also loop through the structures, doing an elementwise assembly
	// in principle similar to the advection business		adv_mats_proj_x_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > (RBsize); // should be RBsize x number_elem_trafo x 2 x RBsize x RBsize
	the_const_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo); // should be number_elem_trafo x 4 x RBsize x RBsize
	the_ABCD_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo);
	the_ABCD_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);	
	the_const_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		the_const_one_proj_2d[i] = Array<OneD, Eigen::MatrixXd > (4);
		the_const_one_rhs_proj_2d[i] = Array<OneD, Eigen::VectorXd > (4);
		the_ABCD_one_proj_2d[i] = Array<OneD, Eigen::MatrixXd > (4);
		the_ABCD_one_rhs_proj_2d[i] = Array<OneD, Eigen::VectorXd > (4);
		for (int j = 0; j < 4; ++j)
		{
			the_const_one_proj_2d[i][j] = Eigen::MatrixXd::Zero(RBsize, RBsize);
			the_const_one_rhs_proj_2d[i][j] = Eigen::VectorXd::Zero(RBsize);
			the_ABCD_one_proj_2d[i][j] = Eigen::MatrixXd::Zero(RBsize, RBsize);
			the_ABCD_one_rhs_proj_2d[i][j] = Eigen::VectorXd::Zero(RBsize);
		}
	}	

	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;
	
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  A_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_bndry
//	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
//	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{

//		cout << "curr_elem d " << curr_elem << endl;

		A_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		B_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		C_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		D_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		Dbnd_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		Dint_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		for (int i = 0; i < 4; i++)
		{
			A_elem[curr_elem][i] = Eigen::MatrixXd::Zero(Ahrows-1, Ahrows-1);
			B_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
			C_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
			D_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
			Dbnd_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_p, nsize_bndry);
			Dint_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_p, nsize_int);
		}
	}

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{

//		cout << "curr_elem " << curr_elem << endl;

/*		int curr_elem_pos = get_curr_elem_pos(curr_elem);
		double detT = Geo_T(w, curr_elem_pos, 0);
		double Ta = Geo_T(w, curr_elem_pos, 1);
		double Tb = Geo_T(w, curr_elem_pos, 2);
		double Tc = Geo_T(w, curr_elem_pos, 3);
		double Td = Geo_T(w, curr_elem_pos, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td; */
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.size();
                int nimap = imap.size();
/*		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = snapshot_x[curr_elem*nphys + i];
			curr_snap_y_part[i] = snapshot_y[curr_elem*nphys + i];
		} */

		Array<OneD, Array<OneD, double> > Ah_ele_vec(4);
		Ah_ele_vec[0] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[1] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[2] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[3] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, Array<OneD, double> > B_ele_vec(4);
		B_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[2] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[3] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > C_ele_vec(4);
		C_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[2] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[3] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > D_ele_vec(4);
		D_ele_vec[0] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		D_ele_vec[1] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		D_ele_vec[2] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		D_ele_vec[3] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > Dbnd_ele_vec(4);
		Dbnd_ele_vec[0] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Dbnd_ele_vec[1] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Dbnd_ele_vec[2] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Dbnd_ele_vec[3] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Array<OneD, Array<OneD, double> > Dint_ele_vec(4);
		Dint_ele_vec[0] = Array<OneD, double> (nsize_p*nsize_int, 0.0);
		Dint_ele_vec[1] = Array<OneD, double> (nsize_p*nsize_int, 0.0);
		Dint_ele_vec[2] = Array<OneD, double> (nsize_p*nsize_int, 0.0);
		Dint_ele_vec[3] = Array<OneD, double> (nsize_p*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					Ah_ele_vec[0][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_0_0[int(bmap[j])];
					Ah_ele_vec[1][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_0_1[int(bmap[j])];
					Ah_ele_vec[2][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_1_0[int(bmap[j])];
					Ah_ele_vec[3][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					B_ele_vec[0][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_0_0[int(imap[j])];
					B_ele_vec[1][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_0_1[int(imap[j])];
					B_ele_vec[2][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_1_0[int(imap[j])];
					B_ele_vec[3][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_1_1[int(imap[j])];
				}
				if (k == 0)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[0][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dbnd_ele_vec[1][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[2][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dbnd_ele_vec[3][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);
			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					C_ele_vec[0][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_0_0[int(bmap[j])];
					C_ele_vec[1][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_0_1[int(bmap[j])];
					C_ele_vec[2][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_1_0[int(bmap[j])];
					C_ele_vec[3][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					D_ele_vec[0][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_0_0[int(imap[j])];
					D_ele_vec[1][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_0_1[int(imap[j])];
					D_ele_vec[2][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_1_0[int(imap[j])];
					D_ele_vec[3][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_1_1[int(imap[j])];
				}
				if (k == 0)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[0][ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dint_ele_vec[1][ (k*nimap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[2][ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dint_ele_vec[3][ (k*nimap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		}

//		cout << "starting writing vec -> mat" << endl;

		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				A_elem[curr_elem][0](i,j) = Ah_ele_vec[0][ i + j*Ahrows ];
				A_elem[curr_elem][1](i,j) = Ah_ele_vec[1][ i + j*Ahrows ];
				A_elem[curr_elem][2](i,j) = Ah_ele_vec[2][ i + j*Ahrows ];
				A_elem[curr_elem][3](i,j) = Ah_ele_vec[3][ i + j*Ahrows ];
			}
		} 
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem][0](i,j) = B_ele_vec[0][ i + j*nsize_bndry ];
				B_elem[curr_elem][1](i,j) = B_ele_vec[1][ i + j*nsize_bndry ];
				B_elem[curr_elem][2](i,j) = B_ele_vec[2][ i + j*nsize_bndry ];
				B_elem[curr_elem][3](i,j) = B_ele_vec[3][ i + j*nsize_bndry ];
			}
		} 
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem][0](i,j) = C_ele_vec[0][ i + j*nsize_bndry ];
				C_elem[curr_elem][1](i,j) = C_ele_vec[1][ i + j*nsize_bndry ];
				C_elem[curr_elem][2](i,j) = C_ele_vec[2][ i + j*nsize_bndry ];
				C_elem[curr_elem][3](i,j) = C_ele_vec[3][ i + j*nsize_bndry ];
			}
		} 
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem][0](i,j) = D_ele_vec[0][ i + j*nsize_int];
				D_elem[curr_elem][1](i,j) = D_ele_vec[1][ i + j*nsize_int];
				D_elem[curr_elem][2](i,j) = D_ele_vec[2][ i + j*nsize_int];
				D_elem[curr_elem][3](i,j) = D_ele_vec[3][ i + j*nsize_int];
			}
		} 
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				Dbnd_elem[curr_elem][0](i,j) = Dbnd_ele_vec[0][ i + j*nsize_p];
				Dbnd_elem[curr_elem][1](i,j) = Dbnd_ele_vec[1][ i + j*nsize_p];
				Dbnd_elem[curr_elem][2](i,j) = Dbnd_ele_vec[2][ i + j*nsize_p];
				Dbnd_elem[curr_elem][3](i,j) = Dbnd_ele_vec[3][ i + j*nsize_p];
			}
		} 
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				Dint_elem[curr_elem][0](i,j) = Dint_ele_vec[0][ i + j*nsize_p];
				Dint_elem[curr_elem][1](i,j) = Dint_ele_vec[1][ i + j*nsize_p];
				Dint_elem[curr_elem][2](i,j) = Dint_ele_vec[2][ i + j*nsize_p];
				Dint_elem[curr_elem][3](i,j) = Dint_ele_vec[3][ i + j*nsize_p];
			}
		}

	} // loop over curr_elem

//	to be build in a seperate projector function:
//	the_const_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo); // should be number_elem_trafo x 4 x RBsize x RBsize
//	the_ABCD_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo);
//	the_ABCD_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);	
//	the_const_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);
	
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			the_const_one_proj_2d[i][j] = press_geo_mat_projector(Dbnd_elem, Dint_elem, i, j, the_const_one_rhs_proj_2d[i][j]);
//			cout << "finished the_const_one_proj_2d " << endl;
			the_ABCD_one_proj_2d[i][j] = ABCD_geo_mat_projector(A_elem, B_elem, C_elem, D_elem, i, j, the_ABCD_one_rhs_proj_2d[i][j]);
		}
	}
	


    }

     Eigen::MatrixXd CoupledLinearNS_ROM::press_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Dbnd_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > Dint_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &press_vec_proj)
    {

	// this function (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects

//	cout << "entering press_geo_mat_projector" << endl;

	Eigen::MatrixXd press_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
        int nz_loc = 1;
        int nel = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
	Eigen::MatrixXd Dbnd_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_bndry*nel );
	Eigen::MatrixXd Dint_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_int*nel );

	time_t timer_1;
	time_t timer_2;
	time_t timer_3;
	time_t timer_4;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		int curr_elem_pos = get_curr_elem_pos(i);
		if (curr_elem_pos == curr_elem_trafo)
		{
			Eigen::MatrixXd curr_Dbnd_elem = Dbnd_elem[i][deriv_index];
			Eigen::MatrixXd curr_Dint_elem = Dint_elem[i][deriv_index];

			Dbnd_all.block(i*nsize_p, i*nsize_bndry, nsize_p, nsize_bndry) = curr_Dbnd_elem;
			Dint_all.block(i*nsize_p, i*nsize_int, nsize_p, nsize_int) = curr_Dint_elem;
	
		}
	}

	// better to use:
//	int f_bnd_size;
//	int f_p_size;
//	int f_int_size;

/*	switch(globally_connected) {
		case 0:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -MtM * Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 1:
			press_matrix.block(0, nBndDofs, nBndDofs, f_p_size) = -Mtrafo.transpose() * Dbnd_all.transpose();
			press_matrix.block(nBndDofs, 0, f_p_size, nBndDofs) = -Dbnd_all * Mtrafo;
			press_matrix.block(nBndDofs, nBndDofs + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(nBndDofs + f_p_size, nBndDofs, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 2:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
	}	
*/
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();

	// also need to create appropriate vector right-hand-side contribution

	Eigen::VectorXd add_to_rhs_press(M_truth_size); // probably need this for adv and non-adv
	add_to_rhs_press = press_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd press_rhs_add = remove_rows(add_to_rhs_press, elem_loc_dbc);
	press_vec_proj = RB.transpose() * press_rhs_add;
	Eigen::MatrixXd press_matrix_simplified = remove_cols_and_rows(press_matrix, elem_loc_dbc);
	Eigen::MatrixXd press_mat_proj = RB.transpose() * press_matrix_simplified * RB;

	return press_mat_proj;
    }

    Eigen::MatrixXd CoupledLinearNS_ROM::ABCD_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > A_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &ABCD_vec_proj)
    {
	Eigen::MatrixXd ABCD_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
        int nz_loc = 1;
        int nel = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel = m_velocity.size();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
	Eigen::MatrixXd A_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	Eigen::MatrixXd B_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd D_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		int curr_elem_pos = get_curr_elem_pos(i);
		if (curr_elem_pos == curr_elem_trafo)
		{
			Eigen::MatrixXd curr_A_elem = A_elem[i][deriv_index];
			Eigen::MatrixXd curr_B_elem = B_elem[i][deriv_index];
			Eigen::MatrixXd curr_C_elem = C_elem[i][deriv_index];
			Eigen::MatrixXd curr_D_elem = D_elem[i][deriv_index];
			A_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = curr_A_elem;
			B_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_B_elem;
			C_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_C_elem;
			D_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = curr_D_elem;
		}
	}
/*	switch(globally_connected) {
		case 0:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -MtM * Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 1:
			press_matrix.block(0, nBndDofs, nBndDofs, f_p_size) = -Mtrafo.transpose() * Dbnd_all.transpose();
			press_matrix.block(nBndDofs, 0, f_p_size, nBndDofs) = -Dbnd_all * Mtrafo;
			press_matrix.block(nBndDofs, nBndDofs + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(nBndDofs + f_p_size, nBndDofs, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 2:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
	}	*/
	/* switch(globally_connected) {
		case 0:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 1:
			ABCD_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_all * Mtrafo;
			ABCD_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_all;
			ABCD_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_all.transpose() * Mtrafo;
			ABCD_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 2:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
	}
	*/
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;

	
	Eigen::VectorXd add_to_rhs_ABCD(M_truth_size);
	add_to_rhs_ABCD = ABCD_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd ABCD_rhs_add = remove_rows(add_to_rhs_ABCD, elem_loc_dbc);
	ABCD_vec_proj = RB.transpose() * ABCD_rhs_add;
	Eigen::MatrixXd ABCD_matrix_simplified = remove_cols_and_rows(ABCD_matrix, elem_loc_dbc);
	Eigen::MatrixXd ABCD_mat_proj = RB.transpose() * ABCD_matrix_simplified * RB;
	return ABCD_mat_proj;
    }

    Eigen::MatrixXd CoupledLinearNS_ROM::gen_affine_mat_proj_geo(double current_nu, double w)
    {
	// avail datastructure	Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > adv_mats_proj_x_2d;
//	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_const_one_proj_2d;
//	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_ABCD_one_proj_2d;
	Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
	Eigen::MatrixXd recovered_press_proj = Eigen::MatrixXd::Zero(RBsize, RBsize);
	Eigen::MatrixXd recovered_ABCD_proj = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int index_elem = 0; index_elem < number_elem_trafo; ++index_elem)
	{
		double detT = Geo_T(w, index_elem, 0);
		double Ta = Geo_T(w, index_elem, 1);
		double Tb = Geo_T(w, index_elem, 2);
		double Tc = Geo_T(w, index_elem, 3);
		double Td = Geo_T(w, index_elem, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td;
		for (int i = 0; i < RBsize; ++i)
		{
//			recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_xy_projected(i,0) + adv_mats_proj_y[i] * curr_xy_projected(i,1);
			recovered_affine_adv_mat_proj_xy += detT * curr_xy_projected(i,0) * (Ta * adv_mats_proj_x_2d[i][index_elem][0] + Tc * adv_mats_proj_x_2d[i][index_elem][1]) + detT * curr_xy_projected(i,1) * (Tb * adv_mats_proj_y_2d[i][index_elem][0] + Td * adv_mats_proj_y_2d[i][index_elem][1]);

		}
		recovered_ABCD_proj += detT * (c00 * the_ABCD_one_proj_2d[index_elem][0] + c01*(the_ABCD_one_proj_2d[index_elem][1] + the_ABCD_one_proj_2d[index_elem][2]) + c11*the_ABCD_one_proj_2d[index_elem][3]);
		recovered_press_proj += detT * (Ta * the_const_one_proj_2d[index_elem][0] + Tc * the_const_one_proj_2d[index_elem][1] + Tb * the_const_one_proj_2d[index_elem][2] + Td * the_const_one_proj_2d[index_elem][3]);

	}

	Eigen::MatrixXd affine_mat_proj = recovered_press_proj + current_nu * recovered_ABCD_proj + recovered_affine_adv_mat_proj_xy;
//	Eigen::MatrixXd affine_mat_proj = the_const_one_proj + current_nu * the_ABCD_one_proj + recovered_affine_adv_mat_proj_xy;

/*	if (debug_mode)
	{
		Eigen::MatrixXd affine_mat = Get_complete_matrix();
		cout << "affine_mat.rows() " << affine_mat.rows() << " affine_mat.cols() " << affine_mat.cols() << endl;
		cout << "affine_mat_proj.rows() " << affine_mat_proj.rows() << " affine_mat_proj.cols() " << affine_mat_proj.cols() << endl;
		cout << "RB.rows() " << RB.rows() << " RB.cols() " << RB.cols() << endl;
		cout << "PODmodes.rows() " << PODmodes.rows() << " PODmodes.cols() " << PODmodes.cols() << endl;
		cout << "affine_mat.norm() "  << affine_mat.norm() << endl;
		cout << "affine_mat_proj.norm() " << affine_mat_proj.norm() << endl;
		Eigen::MatrixXd reproj_affine_mat = Eigen::MatrixXd::Zero(RB.rows(), RB.rows());
		reproj_affine_mat = RB * affine_mat_proj * RB.transpose();
		cout << "reproj_affine_mat.norm() " << reproj_affine_mat.norm() << endl;
		Eigen::MatrixXd affine_matrix_simplified = remove_cols_and_rows(affine_mat, elem_loc_dbc);
		cout << "affine_matrix_simplified.norm() "  << affine_matrix_simplified.norm() << endl;
		Eigen::MatrixXd reduced_affine_mat = Eigen::MatrixXd::Zero(RB.cols(), RB.cols());
		reduced_affine_mat = RB.transpose() * affine_matrix_simplified * RB;
		cout << "reduced_affine_mat.norm() "  << reduced_affine_mat.norm() << endl;

	} */

	return affine_mat_proj;
    }


    Eigen::VectorXd CoupledLinearNS_ROM::gen_affine_vec_proj_geo(double current_nu, double w, int current_index)
    {
	Eigen::VectorXd recovered_affine_adv_rhs_proj_xy = Eigen::VectorXd::Zero(RBsize); 
	Eigen::VectorXd recovered_press_proj = Eigen::VectorXd::Zero(RBsize);
	Eigen::VectorXd recovered_ABCD_proj = Eigen::VectorXd::Zero(RBsize);
	for (int index_elem = 0; index_elem < number_elem_trafo; ++index_elem)
	{
		double detT = Geo_T(w, index_elem, 0);
		double Ta = Geo_T(w, index_elem, 1);
		double Tb = Geo_T(w, index_elem, 2);
		double Tc = Geo_T(w, index_elem, 3);
		double Td = Geo_T(w, index_elem, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td;
		for (int i = 0; i < RBsize; ++i)
		{
	//		recovered_affine_adv_rhs_proj_xy -= adv_vec_proj_x[i] * curr_xy_projected(i,0) + adv_vec_proj_y[i] * curr_xy_projected(i,1);
			recovered_affine_adv_rhs_proj_xy -= detT * curr_xy_projected(i,0) * (Ta * adv_vec_proj_x_2d[i][index_elem][0] + Tc * adv_vec_proj_x_2d[i][index_elem][1]) + detT * curr_xy_projected(i,1) * (Tb * adv_vec_proj_y_2d[i][index_elem][0] + Td * adv_vec_proj_y_2d[i][index_elem][1]);
		}	
		recovered_ABCD_proj += detT * (c00 * the_ABCD_one_rhs_proj_2d[index_elem][0] + c01*(the_ABCD_one_rhs_proj_2d[index_elem][1] + the_ABCD_one_rhs_proj_2d[index_elem][2]) + c11*the_ABCD_one_rhs_proj_2d[index_elem][3]);
		recovered_press_proj += detT * (Ta * the_const_one_rhs_proj_2d[index_elem][0] + Tc * the_const_one_rhs_proj_2d[index_elem][1] + Tb * the_const_one_rhs_proj_2d[index_elem][2] + Td * the_const_one_rhs_proj_2d[index_elem][3]);

	}


	Eigen::VectorXd add_rhs_Newton = Eigen::VectorXd::Zero(RBsize); 

//	cout << " recovered_press_proj " <<  recovered_press_proj  << endl;
//	cout << " recovered_ABCD_proj " <<  recovered_ABCD_proj  << endl;
//	cout << " recovered_affine_adv_rhs_proj_xy " <<  recovered_affine_adv_rhs_proj_xy  << endl;


	return -recovered_press_proj - current_nu * recovered_ABCD_proj + recovered_affine_adv_rhs_proj_xy  -0.5*add_rhs_Newton ;  
//	return -the_const_one_rhs_proj - current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy  -0.5*add_rhs_Newton ;  
    }

    void CoupledLinearNS_ROM::recover_snapshot_loop(Eigen::VectorXd reconstruct_solution, Array<OneD, double> & field_x, Array<OneD, double> & field_y)
    {
	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
	Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
	Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
	Array<OneD, unsigned int> bmap, imap; 
	Array<OneD, double> field_0(GetNcoeffs());
	Array<OneD, double> field_1(GetNcoeffs());
	Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
	Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
	int cnt = 0;
	int cnt1 = 0;
	int nvel = 2;
	int nz_loc = 1;
	int  nplanecoeffs = fields[0]->GetNcoeffs();
	int  nel  = m_fields[0]->GetNumElmts();
	for(int i = 0; i < nel; ++i) 
	{
	      int eid  = i;
	      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
	      fields[0]->GetExp(eid)->GetInteriorMap(imap);
	      int nbnd   = bmap.size();
	      int nint   = imap.size();
	      int offset = fields[0]->GetCoeff_Offset(eid);
	            
	      for(int j = 0; j < nvel; ++j)
	      {
	           for(int n = 0; n < nz_loc; ++n)
	           {
	                    for(int k = 0; k < nbnd; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
	                    }
	                    
	                    for(int k = 0; k < nint; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
	                    }
	                    cnt  += nbnd;
	                    cnt1 += nint;
	           }
	      }
	}
	Array<OneD, double> test_nn = fields[0]->GetCoeffs();
	fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
	fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
	field_x = curr_PhysBaseVec_x;
	field_y = curr_PhysBaseVec_y;
    }
    
    double CoupledLinearNS_ROM::L2norm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y )
    {
	Array< OneD, NekDouble > x_difference(component1_x.size());
	Array< OneD, NekDouble > y_difference(component1_y.size());
	for (int i = 0; i < component1_y.size(); ++i)
	{
		x_difference[i] = component1_x[i] - component2_x[i];
		y_difference[i] = component1_y[i] - component2_y[i];
	}
	double result = L2norm_ITHACA(x_difference, y_difference);
	return result;
    }

    double CoupledLinearNS_ROM::Linfnorm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y )
    {
	Array< OneD, NekDouble > x_difference(component1_x.size());
	Array< OneD, NekDouble > y_difference(component1_y.size());
	for (int i = 0; i < component1_y.size(); ++i)
	{
		x_difference[i] = component1_x[i] - component2_x[i];
		y_difference[i] = component1_y[i] - component2_y[i];
	}
	double result = Linfnorm_ITHACA(x_difference, y_difference);
	return result;
    }    

    double CoupledLinearNS_ROM::L2norm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y )
    {
	// the input comes in phys coords
        NekDouble L2norm = -1.0;
//	cout << "m_NumQuadPointsError " << m_NumQuadPointsError << endl; // should be 0
//	cout << "component_x.size() " << component_x.size() << endl; // should be nphys
        if (m_NumQuadPointsError == 0)
        {
		double L2norm_x = m_fields[0]->L2(component_x);
		double L2norm_y = m_fields[1]->L2(component_y);
		L2norm = sqrt( L2norm_x*L2norm_x + L2norm_y*L2norm_y );
        }
	return L2norm;
    }

    double CoupledLinearNS_ROM::Linfnorm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y )
    {
	// the input comes in phys coords
        NekDouble Linfnorm = -1.0;
//	cout << "m_NumQuadPointsError " << m_NumQuadPointsError << endl; // should be 0
//	cout << "component_x.size() " << component_x.size() << endl; // should be nphys
        if (m_NumQuadPointsError == 0)
        {
		double Linfnorm_x = m_fields[0]->L2(component_x);
		double Linfnorm_y = m_fields[1]->L2(component_y);
		Linfnorm = max(Linfnorm_x, Linfnorm_y);
        }
	return Linfnorm;
    }




   
}
