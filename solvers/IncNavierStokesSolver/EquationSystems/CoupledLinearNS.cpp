///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS.cpp
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
// Description: Coupled  Solver for the Linearised Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS.h>

namespace Nektar
{

    string CoupledLinearNS::className = GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearisedNS", CoupledLinearNS::create);


    /**
     *  @class CoupledLinearNS 
     *
     * Set up expansion field for velocity and pressure, the local to
     * global mapping arrays and the basic memory definitions for
     * coupled matrix system
     */ 
    CoupledLinearNS::CoupledLinearNS(LibUtilities::CommSharedPtr& pComm,
                            LibUtilities::SessionReaderSharedPtr &pSession):
        IncNavierStokes(pComm, pSession)
    {
        int  n,i,j,k,eid;
        int  expdim = m_graph->GetMeshDimension();
        int  n_exp  = m_fields[m_velocity[0]]->GetExpSize();
        int  nvel   = m_velocity.num_elements();

        // Get Expansion list for orthogonal expansion at p-2
        const SpatialDomains::ExpansionMap &pressure_exp = GenPressureExp(m_graph->GetExpansions("u"));

        m_nConvectiveFields = m_fields.num_elements();
        if(NoCaseStringCompare(m_boundaryConditions->GetVariable(m_nConvectiveFields-1),"p") == 0)
        {
            ASSERTL0(false,"Last field is defined as pressure but this is not suitable for this solver, please remove this field as it is implicitly defined");
        }
        // Decide how to declare explist for pressure. 
        if(expdim == 2)
        {
            m_pressure = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(pressure_exp);
            SetUp2DExpansionC0ContMap();
            //SetUp2DExpansionC0ContMap(MultiRegions::eDirectStaticCond);
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

    void CoupledLinearNS::SetUpCoupledMatrix(const NekDouble lambda, const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation)
    {
        int  n,i,j,k,eid;
        int  expdim = m_graph->GetMeshDimension();
        int  n_exp  = m_fields[m_velocity[0]]->GetExpSize();
        int  nvel   = m_velocity.num_elements();

        // if Advfield is defined can assume it is an Oseen or LinearNS equation
        bool AddAdvectionTerms = (Advfield ==  NullNekDoubleArrayofArray)? false: true; 
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

        Array<OneD,unsigned int> nsize_bndry   (n_exp);
        Array<OneD,unsigned int> nsize_bndry_p1(n_exp);
        Array<OneD,unsigned int> nsize_int     (n_exp);
        Array<OneD,unsigned int> nsize_p       (n_exp);
        Array<OneD,unsigned int> nsize_p_m1    (n_exp);
        
        // Set up block matrix sizes - 
        for(n = 0; n < n_exp; ++n)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(n);
            nsize_bndry[n] = nvel*m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs();
            nsize_bndry_p1[n] = nsize_bndry[n]+1;
            nsize_int  [n] = nvel*m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs() - nsize_bndry[n];
            nsize_p[n] = m_pressure->GetExp(eid)->GetNcoeffs();
            nsize_p_m1[n] = nsize_p[n]-1;
        }
            
        MatrixStorage blkmatStorage = eDIAGONAL;
        DNekScalBlkMatSharedPtr pAh = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_bndry_p1,nsize_bndry_p1,blkmatStorage);
        m_BCinv = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
        m_Btilde = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
        m_Cinv = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_int,nsize_int,blkmatStorage);
        
        m_D_bnd = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_p,nsize_bndry,blkmatStorage);

        m_D_int = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_p,nsize_int,blkmatStorage);

        // Final level static condensation matrices. 
        DNekScalBlkMatSharedPtr pBh = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_bndry_p1,nsize_p_m1,blkmatStorage);
        DNekScalBlkMatSharedPtr pCh = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_p_m1,nsize_bndry_p1,blkmatStorage);
        DNekScalBlkMatSharedPtr pDh = MemoryManager<DNekScalBlkMat>
            ::AllocateSharedPtr(nsize_p_m1,nsize_p_m1,blkmatStorage);


        for(n = 0; n < n_exp; ++n)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(n);
            nbndry = nsize_bndry[n];
            nint   = nsize_int[n];
            k = nbndry+1;
            DNekMatSharedPtr Ah = MemoryManager<DNekMat>::AllocateSharedPtr(k,k,zero);
            DNekMatSharedPtr B  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            DNekMatSharedPtr C  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            DNekMatSharedPtr D  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);
            
            
            locExp = m_fields[m_velocity[0]]->GetExp(eid);
            locExp->GetBoundaryMap(bmap);
            locExp->GetInteriorMap(imap);
            LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
                                            locExp->DetExpansionType(),
                                            *locExp,
                                            lambda/m_kinvis);


            int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
            int nbmap = bmap.num_elements();
            int nimap = imap.num_elements(); 

            Array<OneD, NekDouble> coeffs  = m_fields[m_velocity[0]]->GetExp(eid)->UpdateCoeffs();
            Array<OneD, NekDouble> phys    = m_fields[m_velocity[0]]->GetExp(eid)->UpdatePhys();
            int psize   = m_pressure->GetExp(eid)->GetNcoeffs();

            Array<OneD, NekDouble> deriv   = m_pressure->GetExp(eid)->UpdatePhys();
            Array<OneD, NekDouble> pcoeffs = m_pressure->GetExp(eid)->UpdateCoeffs();
            
            
            DNekMatSharedPtr Dbnd = MemoryManager<DNekMat>::AllocateSharedPtr(psize,(k= nvel*bmap.num_elements()));

            k=nvel*imap.num_elements();
            DNekMatSharedPtr Dint = MemoryManager<DNekMat>::AllocateSharedPtr(psize,k);
           

            if(AddAdvectionTerms == false) // use static condensed managed matrices
            {
                // construct velocity matrices using statically
                // condensed elemental matrices and then construct
                // pressure matrix systems
                DNekScalBlkMatSharedPtr CondMat; 
                CondMat = locExp->GetLocStaticCondMatrix(helmkey);
                
                for(k = 0; k < nvel; ++k)
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
                
                for(k = 0; k < nvel; ++k)
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
                
                for(k = 0; k < nvel; ++k)
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
                for(i = 0; i < bmap.num_elements(); ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[bmap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);
                
                    // Differentiation & Inner product wrt base. 
                    for(j = 0; j < nvel; ++j)
                    {
                        m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(j,phys,deriv);
                        StdRegions::StdExpansionSharedPtr locExp = m_pressure->GetExp(eid);
                        locExp->IProductWRTBase(deriv,pcoeffs);
                        // copy into column major storage. 
                        Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                    Dbnd->GetRawPtr() + (j*bmap.num_elements() + i)*psize,1);
                    }
                }

                for(i = 0; i < imap.num_elements(); ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[imap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);
                    
                    // Differentiation & Inner product wrt base. 
                    for(j = 0; j < nvel; ++j)
                    {
                        m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(j,phys,
                                                                        deriv);
                        m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                        // copy into column major storage. 
                        Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                    Dint->GetRawPtr() + (j*imap.num_elements() + i)*psize,1);
                    }
                }
            }
            else
            {
                // construct velocity matrices and pressure systems at
                // the same time resusing differential of velocity
                // space

                DNekScalMat &HelmMat = *locExp->GetLocMatrix(helmkey);
                
                Array<OneD, NekDouble> Advtmp;
                Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel);
                // Use ExpList phys array for temporaary storage
                Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
                int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(eid);
                int nv;
                int npoints = locExp->GetTotPoints();

                // Calculate derivative of base flow 
                if(IsLinearNSEquation)
                {
                    for(nv = 0; nv < nvel; ++nv)
                    {
                        if(nv == 0)
                        {
                            AdvDeriv[nv] = Array<OneD, NekDouble>(npoints);
                        }
                        else
                        {
                            AdvDeriv[nv] = m_fields[nv]->UpdatePhys(); 
                        }
                        m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(nv,Advfield[nv] + phys_offset, AdvDeriv[nv]);
                    }
                }

                for(i = 0; i < nbmap; ++i)
                {
                    
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[bmap[i]] = 1.0;
                    locExp->BwdTrans(coeffs,phys);
                    
                    for(k = 0; k < nvel; ++k)
                    {
                        
                        // Differentiation & Inner product wrt base. 
                        m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(k,phys,
                                                                        deriv);
                     
                        for(j = 0; j < nbmap; ++j)
                        {
                            (*Ah)(i+k*nbmap,j+k*nbmap) += m_kinvis*HelmMat(bmap[i],bmap[j]);
                        }
                        
                        for(j = 0; j < nimap; ++j)
                        {
                            (*B)(i+k*nbmap,j+k*nimap) += m_kinvis*HelmMat(bmap[i],imap[j]);
                        }
                        
                        // Advfield[k] *d/dx_k to all velocity
                        // components on diagonal
                        Vmath::Vmul(npoints, Advtmp = Advfield[k] + phys_offset,
                                    1,deriv,1,tmpphys,1);
                        locExp->IProductWRTBase(tmpphys,coeffs);


                        for(nv = 0; nv < nvel; ++nv)
                        {
                            for(j = 0; j < nbmap; ++j)
                            {
                                (*Ah)(j+nv*nbmap,i+nv*nbmap) +=
                                    coeffs[bmap[j]];
                            }
                            
                            for(j = 0; j < nimap; ++j)
                            {
                                (*C)(i+nv*nbmap,j+nv*nimap) += 
                                    coeffs[imap[j]];
                            }
                        }
  
                        if(IsLinearNSEquation)
                        {
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1, AdvDeriv[nv],1,
                                            tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(j = 0; j < nbmap; ++j)
                                {
                                    (*Ah)(j+k*nbmap,i+nv*nbmap) +=
                                        coeffs[bmap[j]];
                                }
                                
                                for(j = 0; j < nimap; ++j)
                                {
                                    (*C)(i+nv*nbmap,j+k*nimap) += 
                                        coeffs[imap[j]];
                                }
                            }                            
                        }
  
                        // copy into column major storage. 
                        m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                        Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                    Dbnd->GetRawPtr() + 
                                    (k*bmap.num_elements() + i)*psize,1);
                    }
                }
                    

                for(i = 0; i < nimap; ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[imap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);
                    
                    for(k = 0; k < nvel; ++k)
                    {
                        // Differentiation & Inner product wrt base. 
                        m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(k,phys,deriv);
                        
                        for(j = 0; j < nbmap; ++j) // C set up as transpose
                        {
                            (*C)(j+k*nbmap,i+k*nimap) += m_kinvis*HelmMat(imap[i],bmap[j]);
                        }
                        
                        for(j = 0; j < nimap; ++j)
                        {
                            (*D)(i+k*nimap,j+k*nimap) += m_kinvis*HelmMat(imap[i],imap[j]);
                        }

                        // Advfield[k] *d/dx_k to all velocity
                        // components on diagonal
                        Vmath::Vmul(npoints, Advtmp = Advfield[k] + phys_offset,
                                    1,deriv,1,tmpphys,1);
                        locExp->IProductWRTBase(tmpphys,coeffs);

                        for(nv = 0; nv < nvel; ++nv)
                        {
                            for(j = 0; j < nbmap; ++j)
                            {
                                (*B)(j+nv*nbmap,i+nv*nimap) += 
                                    coeffs[bmap[j]];
                            }

                            for(j = 0; j < nimap; ++j)
                            {
                                (*D)(j+nv*nimap,i+nv*nimap) += 
                                    coeffs[imap[j]];
                            }
                        }

                        if(IsLinearNSEquation)
                        {
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1, AdvDeriv[nv],1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);

                                for(j = 0; j < nbmap; ++j)
                                {
                                    (*B)(j+k*nbmap,i+nv*nimap) += 
                                        coeffs[bmap[j]];
                                }

                                for(j = 0; j < nimap; ++j)
                                {
                                    (*D)(j+k*nimap,i+nv*nimap) += 
                                        coeffs[imap[j]];
                                }
                            }
                        }

                        // copy into column major storage. 
                        m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                        Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                    Dint->GetRawPtr() +
                                    (k*imap.num_elements() + i)*psize,1);
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
            
            m_BCinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
            m_Btilde->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,C));
            m_Cinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,D));
            m_D_bnd->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dbnd));
            m_D_int->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dint));
            
            // Do matrix manipulations and get final set of block matries    
            // reset boundary to put mean mode into boundary system. 
            
            DNekMatSharedPtr Cinv,BCinv,Btilde; 
            DNekMat  DintCinvDTint, BCinvDTint_m_DTbnd, DintCinvBTtilde_m_Dbnd;

            Cinv   = D;
            BCinv  = B;  
            Btilde = C; 

            DintCinvDTint = (*Dint)*(*Cinv)*Transpose(*Dint);
            
            BCinvDTint_m_DTbnd =  (*BCinv)*Transpose(*Dint) - Transpose(*Dbnd);
            
            // This could be transpose of BCinvDint in some cases
            DintCinvBTtilde_m_Dbnd =(*Dint)*(*Cinv)*Transpose(*Btilde) -  (*Dbnd); 
            
            // Set up final set of matrices. 
            DNekMatSharedPtr Bh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_bndry_p1[n],nsize_p_m1[n]);
            DNekMatSharedPtr Ch = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n],nsize_bndry_p1[n]);
            DNekMatSharedPtr Dh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n], nsize_p_m1[n]);
            

            // Copy matrices into final structures. 
            for(i = 0; i < nsize_p_m1[n]; ++i)
            {
                for(j = 0; j < nsize_p_m1[n]; ++j)
                {
                    (*Dh)(i,j) = -DintCinvDTint(i+1,j+1);
                }
            }
            
            for(i = 0; i < nsize_bndry_p1[n]-1; ++i)
            {
                (*Ah)(i,nsize_bndry_p1[n]-1) = BCinvDTint_m_DTbnd(i,0);
                (*Ah)(nsize_bndry_p1[n]-1,i) = DintCinvBTtilde_m_Dbnd(0,i);
            }
            
            (*Ah)(nsize_bndry_p1[n]-1,nsize_bndry_p1[n]-1) = -DintCinvDTint(0,0);

            
            for(j = 0; j < nsize_p_m1[n]; ++j)
            {
                for(i = 0; i < nsize_bndry_p1[n]-1; ++i)
                {
                    (*Bh)(i,j) = BCinvDTint_m_DTbnd(i,j+1);
                    (*Ch)(j,i) = DintCinvBTtilde_m_Dbnd(j+1,i);
                }
                (*Bh)(nsize_bndry_p1[n]-1,j) = -DintCinvDTint(0,j+1);
                (*Ch)(j,nsize_bndry_p1[n]-1) = -DintCinvDTint(j+1,0);
            }

            // Do static condensation
            Dh->Invert();
            (*Bh) = (*Bh)*(*Dh);
            (*Ah) = (*Ah) - (*Bh)*(*Ch);

            // Set matrices for later inversion. Probably do not need to be 
            // attached to class
            pAh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ah));
            pBh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Bh));
            pCh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ch));
            pDh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Dh));    
        
        }

        // Set up global coupled boundary sovler. 
        // This is a key to define the solution matrix type
        // currently we are giving it a argument of eLInearAdvectionReaction 
        // since this then makes the matrix storage of type eFull
        MultiRegions::GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,m_locToGloMap);
        m_CoupledBndSys = MemoryManager<MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(key,m_fields[0],pAh,pBh,pCh,pDh,m_locToGloMap);
    }

    void CoupledLinearNS::v_PrintSummary(std::ostream &out)
    {
        cout <<  "\tSolver Type     : Coupled Linearised NS" <<endl;
    }

    void CoupledLinearNS::v_DoInitialise(void)
    {
        switch(m_equationType)
        {
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            {

                LibUtilities::TimeIntegrationMethod intMethod;
                std::string TimeIntStr = m_boundaryConditions->GetSolverInfo("TIMEINTEGRATIONMETHOD");
                int i;
                for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
                {
                    if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
                    {
                        intMethod = (LibUtilities::TimeIntegrationMethod)i; 
                            break;
                    }
                }
                
                ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
                
                switch(intMethod)
                {
                case LibUtilities::eIMEXOrder1: 
                    {
                        m_intSteps = 1;
                        m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                        LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
                        m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    }
                    break;
                case LibUtilities::eIMEXOrder2: 
                    {
                        m_intSteps = 2;
                        m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                        LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXOrder1);
                        m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                        LibUtilities::TimeIntegrationSchemeKey       IntKey1(intMethod);
                        m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    }
                    break;
                default:
                    ASSERTL0(0,"Integration method not setup: Options include ImexOrder1, ImexOrder2");
                    break;
                }
                
                // Could defind this from IncNavierStokes class? 
                m_integrationOps.DefineOdeRhs(&CoupledLinearNS::EvaluateAdvection, this);

                m_integrationOps.DefineImplicitSolve(&CoupledLinearNS::SolveUnsteadyStokesSystem,this);

            }
        case eSteadyStokes:
            SetUpCoupledMatrix();
            break;
        case eSteadyOseen:
            {
                
                Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.num_elements());
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }

                EvaluateUserDefinedEqn(AdvField); 
                
                SetUpCoupledMatrix(0.0,AdvField,false);
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS");
        }
    }

    void CoupledLinearNS::EvaluateAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                            Array<OneD, Array<OneD, NekDouble> > &outarray,
                                            const NekDouble time)
    {
        // evaluate convection terms
        EvaluateAdvectionTerms(inarray,outarray);
        int nqtot  =m_fields[0]->GetTotPoints(); 
        //add the force
	if(m_bforce)
	{	
	  for(int i = 0; i < m_nConvectiveFields; ++i)
          {
			 Vmath::Vadd(nqtot,outarray[i],1,(m_forces[i]->GetPhys()),1,outarray[i],1);
          }        
        }
    }

    void CoupledLinearNS::SolveUnsteadyStokesSystem(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                    const NekDouble time, 
                                                    const NekDouble aii_Dt)
    {

        int i,n;
        Array<OneD, Array< OneD, NekDouble> > F(m_nConvectiveFields);
        NekDouble  lambda = 1.0/aii_Dt; 
        static NekDouble lambda_store; 
        Array <OneD, Array<OneD, NekDouble> > forcing(m_velocity.num_elements());
        // Matrix solution
        if(fabs(lambda_store - lambda) > 1e-10)
        {
            cout << "Setting up Stokes matrix problem [.";
            fflush(stdout);
            SetUpCoupledMatrix(lambda);
            cout << "]" << endl;
            lambda_store = lambda;
        }
        
        SetBoundaryConditions(time);		

        // Forcing for advection solve 
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
            m_fields[m_velocity[i]]->IProductWRTBase(inarray[i],m_fields[m_velocity[i]]->UpdateCoeffs());
            Vmath::Smul(m_fields[m_velocity[i]]->GetNcoeffs(),lambda,m_fields[m_velocity[i]]->GetCoeffs(), 1,m_fields[m_velocity[i]]->UpdateCoeffs(),1);
            forcing[i] = m_fields[m_velocity[i]]->GetCoeffs();
        }

        SolveLinearNS(forcing);

        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),outarray[i]);
        }
    }
    
    void CoupledLinearNS::v_DoSolve(void)
    {
        switch(m_equationType)
        {
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            // Set initial condition using time t=0
            SetInitialConditions(0.0);

            AdvanceInTime(m_steps);
            break;
        case eSteadyStokes:
        case eSteadyOseen:
            Solve();
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS");
        }
    }

    void CoupledLinearNS::Solve(void)
    {
        Array <OneD, Array<OneD, NekDouble> > forcing(m_velocity.num_elements());

        // Should put in read forcing in here. 
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
            forcing[i] = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
        }
        
        SolveLinearNS(forcing);
    }

    const SpatialDomains::ExpansionMap &CoupledLinearNS::GenPressureExp(const SpatialDomains::ExpansionMap &VelExp)
    {
        int i;
        SpatialDomains::ExpansionMapShPtr returnval;

        returnval = MemoryManager<SpatialDomains::ExpansionMap>::AllocateSharedPtr();
        
        SpatialDomains::ExpansionMap::const_iterator  expMapIter;
        int nummodes;

        for (expMapIter = VelExp.begin(); expMapIter != VelExp.end(); ++expMapIter)
        {
            LibUtilities::BasisKeyVector BasisVec;

            for(i = 0; i <  expMapIter->second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey B = expMapIter->second->m_basisKeyVector[i];
                nummodes = B.GetNumModes();
                ASSERTL0(nummodes > 3,"Velocity polynomial space not sufficiently high (>= 4)");
                // Should probably set to be an orthogonal basis. 
                LibUtilities::BasisKey newB(B.GetBasisType(),nummodes-2,B.GetPointsKey());
                BasisVec.push_back(newB);
            }

            // Put new expansion into list. 
            SpatialDomains::ExpansionShPtr expansionElementShPtr =
                MemoryManager<SpatialDomains::Expansion>::AllocateSharedPtr(expMapIter->second->m_geomShPtr, BasisVec);
            (*returnval)[expMapIter->first] = expansionElementShPtr;
        }
        
        // Save expansion into graph. 
        m_graph->SetExpansions("p",returnval);
            
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
    void CoupledLinearNS::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing)
    {
        int i,j,k,eid,cnt,cnt1;
        int nbnd,nint,offset;
        int nvel = m_velocity.num_elements();
        int nel = m_fields[m_velocity[0]]->GetExpSize();
        Array<OneD, unsigned int> bmap, imap; 

        Array<OneD, NekDouble> f_bnd(m_BCinv->GetRows());
        NekVector< NekDouble > F_bnd(f_bnd.num_elements(), f_bnd, eWrapper);

        Array<OneD, NekDouble> f_int(m_BCinv->GetColumns());
        NekVector< NekDouble > F_int(f_int.num_elements(),f_int, eWrapper);

        // Assemble f_bnd and f_int
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            m_fields[m_velocity[0]]->GetExp(eid)->GetBoundaryMap(bmap);
            m_fields[m_velocity[0]]->GetExp(eid)->GetInteriorMap(imap);
            nbnd   = bmap.num_elements();
            nint   = imap.num_elements();
            offset = m_fields[m_velocity[0]]->GetCoeff_Offset(eid);

            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(k = 0; k < nbnd; ++k)
                {
                    f_bnd[cnt+k] = forcing[m_velocity[j]][offset+bmap[k]];
                }
                for(k = 0; k < nint; ++k)
                {
                    f_int[cnt1+k] = forcing[m_velocity[j]][offset+imap[k]];
                }
                cnt  += nbnd;
                cnt1 += nint;
            }
        }

        Array<OneD, NekDouble > f_p(m_D_int->GetRows());
        NekVector<  NekDouble > F_p(f_p.num_elements(),f_p,eWrapper);

        // fbnd does not currently hold the pressure mean
        F_bnd = F_bnd - (*m_BCinv)*F_int;
        F_p   = (*m_D_int)*((*m_Cinv)*F_int);
        
        // construct inner forcing 
        Array<OneD, NekDouble > bnd   (m_locToGloMap->GetNumGlobalCoeffs(),0.0);
        Array<OneD, NekDouble > fh_bnd(m_locToGloMap->GetNumGlobalCoeffs(),0.0);
        
        const Array<OneD,const int>& loctoglomap
            = m_locToGloMap->GetLocalToGlobalMap();
        const Array<OneD,const NekDouble>& loctoglosign
            = m_locToGloMap->GetLocalToGlobalSign();

        // no sign change set up for this type of expansion so define
        // default as a vector of 1.0
        if(loctoglosign == NullNekDouble1DArray)
        {
            m_locToGloMap->SetLocalToGlobalSign(Array<OneD, const NekDouble>(loctoglomap.num_elements(),1.0));
        }
        
        offset = cnt1 = cnt = 0; 
        for(i = 0; i < nel; ++i)
        {
            eid  = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            nbnd = m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    fh_bnd[loctoglomap[offset+j*nbnd+k]] += 
                        loctoglosign[offset+j*nbnd+k]*f_bnd[cnt+k];
                }
                cnt += nbnd;
            }
            fh_bnd[loctoglomap[offset + nvel*nbnd]] = f_p[cnt1];
            nint    = m_pressure->GetExp(eid)->GetNcoeffs();
            cnt1   += nint;
            offset += nvel*nbnd + nint; 
        }

        offset = cnt1 = 0; 
        for(i = 0; i <  nel; ++i)
        {
            eid  = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            nbnd = m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs(); 
            nint = m_pressure->GetExp(eid)->GetNcoeffs(); 

            for(j = 0; j < nint-1; ++j)
            {
                fh_bnd[loctoglomap[offset + nvel*nbnd+1 + j]] += f_p[cnt1+1+j];
            }
            cnt1 += nint;
            offset += nvel*nbnd + nint; 
        }

        //  Set Weak BC into f_bnd and Dirichlet Dofs in bnd
        const Array<OneD,const int>& bndmap
            = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
        
        // Forcing function with weak boundary conditions
        int bndcnt=0;
        
        for(k = 0; k < nvel; ++k)
        {
            const Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp = m_fields[m_velocity[k]]->GetBndCondExpansions();
            const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConds = m_fields[m_velocity[k]]->GetBndConditions();
            
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                if(bndConds[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                    {
                        bnd[bndmap[bndcnt++]]
                            = (bndCondExp[i]->GetCoeffs())[j];
                    }
                }
                else
                {                    
                    for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                    {
                        fh_bnd[bndmap[bndcnt++]]
                            += (bndCondExp[i]->GetCoeffs())[j];
                    }
                }
            }
        }
        
        m_CoupledBndSys->Solve(fh_bnd,bnd,m_locToGloMap);

        // unpack pressure and velocity boundary systems. 
        Array<OneD, NekDouble> p_coeffs = m_pressure->UpdateCoeffs();
        offset = cnt1 = cnt = 0; 
        for(i = 0; i <  nel; ++i)
        {
            eid  = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            nbnd = m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs(); 
            nint = m_pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                                   	
                    f_bnd[cnt+k] = loctoglosign[offset+j*nbnd+k]*bnd[loctoglomap[offset + j*nbnd + k]];
                }
                cnt += nbnd;
            }
            
            p_coeffs[m_pressure->GetCoeff_Offset(eid)] = f_p[cnt1] = bnd[loctoglomap[offset + nvel*nbnd]];
            cnt1 += nint;
            offset += nvel*nbnd + nint;
        }

        m_pressure->SetPhysState(false);

        offset = cnt = cnt1 = 0;
        for(i = 0; i <  nel; ++i)
        {
            eid  = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            nint = m_pressure->GetExp(eid)->GetNcoeffs(); 
            nbnd = m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs(); 
            cnt1 = m_pressure->GetCoeff_Offset(eid);
            for(j = 0; j < nint-1; ++j)
            {
            	
            	    
                p_coeffs[cnt1+1+j] = 
                    f_p[cnt+1+j] = bnd[loctoglomap[offset + nvel*nbnd + 1 +j]];

            }
            cnt    += nint;
            offset += nvel*nbnd + nint; 
        }

        // Back solve first level of static condensation for interior
        // velocity space and store in F_int
        F_int = (*m_Cinv)*(F_int + Transpose(*m_D_int)*F_p - Transpose(*m_Btilde)*F_bnd);
                   
              
        // Unpack solution from Bnd and  F_int to m_fields 
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            eid  = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            m_fields[m_velocity[0]]->GetExp(eid)->GetBoundaryMap(bmap);
            m_fields[m_velocity[0]]->GetExp(eid)->GetInteriorMap(imap);
            nbnd   = bmap.num_elements();
            nint   = imap.num_elements();
            offset = m_fields[m_velocity[0]]->GetCoeff_Offset(eid);

            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(k = 0; k < nbnd; ++k)
                {
                    m_fields[m_velocity[j]]->SetCoeff(offset+bmap[k],f_bnd[cnt+k]);
                }
                for(k = 0; k < nint; ++k)
                {
                    m_fields[m_velocity[j]]->SetCoeff(offset+imap[k],f_int[cnt1+k]);
                }
                cnt  += nbnd;
                cnt1 += nint;
            }
        }
        
        for(j = 0; j < nvel; ++j) 
        {
            m_fields[m_velocity[j]]->SetPhysState(false);
        }
    }
    
    /** 
     * This is an vector extension of
     * MultiRegion::LocalToGlobalC0BaseMap::SetUp2DExpansionC0ContMap
     * related to the Linearised Navier Stokes problem
     */
    void CoupledLinearNS::SetUp2DExpansionC0ContMap(const MultiRegions::GlobalSysSolnType solnType)
    {

        int i,j,k;
        int cnt = 0,offset=0;
        int meshVertId, meshVertId2;
        int meshEdgeId, meshEdgeId2;
        int bndEdgeCnt;
        int globalId;
        int nEdgeCoeffs;
        int nEdgeInteriorCoeffs;
        int firstNonDirGraphVertId;
        int nLocBndCondDofs = 0;
        int nLocDirBndCondDofs = 0;
        StdRegions::StdExpansion2DSharedPtr locExpansion;
        LocalRegions::SegExpSharedPtr       bndSegExp;
        LibUtilities::BasisType             bType;
        StdRegions::EdgeOrientation         edgeOrient;
        Array<OneD, unsigned int>           edgeInteriorMap;
        Array<OneD, int>                    edgeInteriorSign;
        bool signChange = false;
        int nvel = m_velocity.num_elements();
        
        const StdRegions::StdExpansionVector &locExpVector = *(m_fields[m_velocity[0]]->GetExp());
        int eid, id, diff;
        int nel = locExpVector.size();
        
        map<int,int> periodicEdges;
        vector<map<int,int> > periodicVertices;
        Array<OneD, map<int,int> > ReorderedGraphVertId(2);
        MultiRegions::BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
        int staticCondLevel = 0;
        bool IsSingular = true; 

        /**
         * STEP 1: Wrap boundary conditions vector in an array
         * (since routine is set up for multiple fields) and call
         * the graph re-odering subroutine to obtain the reordered
         * values
         */

        // Obtain any periodic information and allocate default mapping array
        SpatialDomains::MeshGraph2DSharedPtr graph2D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D> (m_graph);
        m_fields[m_velocity[0]]->GetPeriodicEdges(*graph2D,*m_boundaryConditions,
                                                  m_boundaryConditions->GetVariable(0),
                                                  periodicVertices,periodicEdges);
        
        MultiRegions::LocalToGlobalC0ContMapSharedPtr locToGloMap;
        m_locToGloMap = locToGloMap = MemoryManager<MultiRegions::LocalToGlobalC0ContMap>::AllocateSharedPtr(m_comm);

        const Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp = m_fields[m_velocity[0]]->GetBndCondExpansions();
        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> > bndConditionsVec(nvel);

        map<int,int> IsDirVertDof;
        map<int,int> IsDirEdgeDof;
        map<int,int>::iterator mapIt;

        for(i = 0; i < nvel; ++i)
        {
            bndConditionsVec[i] = m_fields[m_velocity[i]]->GetBndConditions(); 
        }
        
        for(j = 0; j < bndCondExp.num_elements(); ++j)
        {
            map<int,int> BndExpVids;
            // collect unique list of vertex ids for this expansion
            for(k = 0; k < bndCondExp[j]->GetExpSize(); ++k)
            {
                BndExpVids[bndCondExp[j]->GetExp(k)->GetGeom1D()->GetVid(0)] = bndCondExp[j]->GetExp(k)->GetGeom1D()->GetVid(0);
                BndExpVids[bndCondExp[j]->GetExp(k)->GetGeom1D()->GetVid(1)] = bndCondExp[j]->GetExp(k)->GetGeom1D()->GetVid(1);
            }
            
            for(i = 0; i < nvel; ++i)
            {
                if(bndConditionsVec[i][j]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                {
                    // set number of Dirichlet conditions along edge
                    for(k = 0; k < bndCondExp[j]->GetExpSize(); ++k)
                    {
                        IsDirEdgeDof[bndCondExp[j]->GetExp(k)->GetGeom1D()->GetEid()] += 1;
                    }

                            
                    // Set number of Dirichlet conditions at vertices
                    // with a clamp on its maximum value being nvel to
                    // handle corners between expansions
                    for(mapIt = BndExpVids.begin(); mapIt !=  BndExpVids.end(); mapIt++)
                    {
                        id = IsDirVertDof[mapIt->second]+1; 
                        IsDirVertDof[mapIt->second] = (id > nvel)?nvel:id;
                    }
                }
                else
                { // if not all dirichlet then not singular
                    IsSingular = false;
                }
            }
        }

        Array<OneD, map<int,int> >Dofs(2);

        Array<OneD, int> AddMeanPressureToEdgeId(locExpVector.size(),-1);
        int edgeId,vertId;


        // special case of singular problem - need to fix one
        // pressure dof to and dirichlet edge 
        if(IsSingular)
        {
            id = bndCondExp[0]->GetExp(0)->GetGeom1D()->GetEid(); 

            // determine element with this edge id. There may be a
            // more direct way of getting element from spatialDomains
            for(i = 0; i < locExpVector.size(); ++i)
            {
                for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
                {
                    edgeId = (locExpVector[i]->GetGeom2D())->GetEid(j);
                    
                    if(edgeId == id)
                    {
                        AddMeanPressureToEdgeId[i] = id;
                        break;
                    }
                }
                
                if(AddMeanPressureToEdgeId[i] != -1)
                {
                    break;
                }
            }
        }


        for(i = 0; i < locExpVector.size(); ++i)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            for(j = 0; j < locExpVector[eid]->GetNverts(); ++j)
            {
                vertId = (locExpVector[eid]->GetGeom2D())->GetVid(j);
                if(Dofs[0].count(vertId) == 0)
                {
                    Dofs[0][vertId] = nvel;
                    
                    // Adjust for a mixed boundary condition
                    if(IsDirVertDof.count(vertId) != 0)
                    {
                        int diff = nvel-IsDirVertDof[vertId];
                        Dofs[0][vertId] -= diff; 
                    }
                }

                edgeId = (locExpVector[eid]->GetGeom2D())->GetEid(j);
                if(Dofs[1].count(edgeId) == 0)
                {
                    Dofs[1][edgeId] = nvel*(locExpVector[eid]->GetEdgeNcoeffs(j)-2);
                }
                
                // count how many vertices 
                if(IsDirEdgeDof.count(edgeId) != 0)
                {
                    int diff = nvel-IsDirEdgeDof[edgeId];
                    Dofs[1][edgeId] -= diff*(locExpVector[i]->GetEdgeNcoeffs(j)-2);
                }
            }
        }
        
        locToGloMap->SetUp2DGraphC0ContMap(*m_fields[m_velocity[0]],
                                           solnType,
                                           bndCondExp,
                                           bndConditionsVec,
                                           periodicVertices,
                                           periodicEdges,
                                           Dofs,
                                           ReorderedGraphVertId,
                                           firstNonDirGraphVertId,
                                           bottomUpGraph,
                                           false,
                                           4);


        /**
         * STEP 2a: Set the mean pressure modes to edges depending on
         * type of direct solver technique;
         */

        // determine which edge to add mean pressure dof based on
        // ensuring that at least one pressure dof from an internal
        // patch is associated with its boundary system
        if((solnType == MultiRegions::eDirectMultiLevelStaticCond))
        {
            // Should put this in a separate function!!

            // Make list of homogeneous graph edges to elmt mappings
            Array<TwoD, int> EdgeIdToElmts(ReorderedGraphVertId[1].size(),2,-1);
            map<int,int> HomGraphEdgeIdToEdgeId;
            
            for(i = 0; i < locExpVector.size(); ++i)
            {
                for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
                {
                    edgeId = (locExpVector[i]->GetGeom2D())->GetEid(j);
                  
                    if(ReorderedGraphVertId[1][edgeId] >= firstNonDirGraphVertId)
                    {
                        HomGraphEdgeIdToEdgeId[ReorderedGraphVertId[1][edgeId]-firstNonDirGraphVertId] = edgeId;
                        if(EdgeIdToElmts[edgeId][0] == -1)
                        {
                            EdgeIdToElmts[edgeId][0] = i;
                        }
                        else
                        {
                            EdgeIdToElmts[edgeId][1] = i;
                        }
                    }
                }
            }
        
            map<int,int>::iterator mapIt;
            
            // Start at second to last level and find edge on boundary
            // to attach element
            int nlevels = bottomUpGraph->GetNlevels();

            // determine a default edge to attach pressure modes to
            // which is part of the inner solve;
            int defedge = -1;
            vector<MultiRegions::SubGraphSharedPtr> bndgraphs = bottomUpGraph->GetInteriorBlocks(nlevels);
            for(i = 0; i < bndgraphs.size(); ++i)
            {
                int GlobIdOffset = bndgraphs[i]->GetIdOffset();
                
                for(j = 0; j < bndgraphs[i]->GetNverts(); ++j)
                {
                    // find edge in graph vert list
                    if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
                    {
                        edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];
                        if(defedge == -1)
                        {
                            defedge = edgeId;
                            break;
                        }
                    }
                }
                if(defedge != -1)
                {
                    break;
                }
            }

            for(int n = 1; n < nlevels; ++n)
            {
                // produce a map with a key that is the element id
                // that contains which next level patch it belongs to
                vector<MultiRegions::SubGraphSharedPtr> bndgraphs = bottomUpGraph->GetInteriorBlocks(n+1);

                // Fill next level graph  of adjacent elements and their level
                map<int,int> ElmtInBndry;

                for(i = 0; i < bndgraphs.size(); ++i)
                {
                    int GlobIdOffset = bndgraphs[i]->GetIdOffset();
                    
                    for(j = 0; j < bndgraphs[i]->GetNverts(); ++j)
                    {
                        // find edge in graph vert list
                        if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
                        {
                            edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];
                        

                            if(EdgeIdToElmts[edgeId][0] != -1)
                            {
                                ElmtInBndry[EdgeIdToElmts[edgeId][0]] = i;
                            }
                            if(EdgeIdToElmts[edgeId][1] != -1)
                            {
                                ElmtInBndry[EdgeIdToElmts[edgeId][1]] = i;
                            }
                        }
                    }
                }

                // Now search interior patches in this level for edges
                // that share the same element as a boundary edge and
                // assign this elmt that boundary edge
                vector<MultiRegions::SubGraphSharedPtr> intgraphs = bottomUpGraph->GetInteriorBlocks(n);
                for(i = 0; i < intgraphs.size(); ++i)
                {
                    int GlobIdOffset = intgraphs[i]->GetIdOffset();
                    bool SetEdge = false; 
                    int elmtid;
                    for(j = 0; j < intgraphs[i]->GetNverts(); ++j)
                    {
                        // Check to see if graph vert is an edge 
                        if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
                        {
                            edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];

                            for(k = 0; k < 2; ++k)
                            {
                                // relevant edge id
                                elmtid = EdgeIdToElmts[edgeId][k];

                                if(elmtid != -1)
                                {
                                    mapIt = ElmtInBndry.find(elmtid);
                                    
                                    if(mapIt != ElmtInBndry.end())
                                    {
                                        // now find a edge in the next level boundary graph 
                                        int GlobIdOffset1 = bndgraphs[mapIt->second]->GetIdOffset();
                                        for(int l = 0; l < bndgraphs[mapIt->second]->GetNverts(); ++l)
                                        {
                                            // find edge in graph vert list
                                            if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset1+l) != 0)
                                            {
                                                AddMeanPressureToEdgeId[elmtid] = HomGraphEdgeIdToEdgeId[GlobIdOffset1+l];
                        
                                                SetEdge = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // if we have failed to find matching edge in next
                    // level patch boundary then set last found elmt
                    // associated to this interior patch to the
                    // default edget value
                    if(SetEdge == false)
                    {
                        if(elmtid == -1) // find a elmtid in patch 
                        {
                            for(j = 0; j < intgraphs[i]->GetNverts(); ++j)
                            {
                                if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
                                {
                                    edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];
                                    for(k = 0; k < 2; ++k)
                                    {
                                        // relevant edge id
                                        elmtid = EdgeIdToElmts[edgeId][k];
                                        if(elmtid != -1)
                                        {
                                            break;
                                        }
                                    }
                                }
                                if(elmtid != -1)
                                {
                                    break;
                                }
                            }
                        }

                        AddMeanPressureToEdgeId[elmtid] = defedge;
                    }
                }
            }
        }

        // Set unset elmts to non-Dirichlet edges. 
        // special case of singular problem - need to fix one
        // pressure dof to and dirichlet edge 
        for(i = 0; i < locExpVector.size(); ++i)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            for(j = 0; j < locExpVector[eid]->GetNverts(); ++j)
            {
                edgeId = (locExpVector[eid]->GetGeom2D())->GetEid(j);

                if(IsDirEdgeDof.count(edgeId) == 0) // interior edge
                {
                    // setup AddMeanPressureToEdgeId to decide where to
                    // put pressure
                    if(AddMeanPressureToEdgeId[eid] == -1)
                    {
                        AddMeanPressureToEdgeId[eid] = edgeId;
                    }
                }
            }
            ASSERTL0((AddMeanPressureToEdgeId[eid] != -1),"Did not determine "
                     "a  edge to attach mean pressure dof");
            // Add the mean pressure degree of freedom to this edge
            Dofs[1][AddMeanPressureToEdgeId[eid]] += 1;
        }
        
        map<int,int> pressureEdgeOffset;
                
        /**
         * STEP 2: Count out the number of Dirichlet vertices and edges first
         */
        for(i = 0; i < bndCondExp.num_elements(); i++)
        {
            for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
            {
                bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));
                for(k = 0; k < nvel; ++k)
                {
                    if(bndConditionsVec[k][i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {
                        nLocDirBndCondDofs += bndSegExp->GetNcoeffs();
                    }
                    nLocBndCondDofs += bndSegExp->GetNcoeffs();
                }
            }
        }

        if(IsSingular)
        {
            locToGloMap->SetNumLocalDirBndCoeffs(nLocDirBndCondDofs+1);
        }
        else
        {
            locToGloMap->SetNumLocalDirBndCoeffs(nLocDirBndCondDofs);
        }

        /**
         * STEP 3: Set up an array which contains the offset information of
         * the different graph vertices.
         *
         * This basically means to identify how many global degrees of
         * freedom the individual graph vertices correspond. Obviously,
         * the graph vertices corresponding to the mesh-vertices account
         * for a single global DOF. However, the graph vertices
         * corresponding to the element edges correspond to 2*(N-2) global DOF
         * where N is equal to the number of boundary modes on this edge.
         */
        Array<OneD, int> graphVertOffset(nvel*(ReorderedGraphVertId[0].size()+
                                               ReorderedGraphVertId[1].size()),0);
        graphVertOffset[0] = 0;
        
        for(i = 0; i < nel; ++i)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[eid]);

            for(j = 0; j < locExpansion->GetNedges(); ++j)
            {
                nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);
                meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                meshVertId = (locExpansion->GetGeom2D())->GetVid(j);

                for(k = 0; k < nvel; ++k)
                {
                    graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nvel+k] = (nEdgeCoeffs-2);
                    graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nvel+k] = 1;
                }

                bType = locExpansion->GetEdgeBasisType(j);
                // need a sign vector for modal expansions if nEdgeCoeffs >=4
                if( (nEdgeCoeffs >= 4)&&
                    ( (bType == LibUtilities::eModified_A)||
                      (bType == LibUtilities::eModified_B) ) )
                {
                    signChange = true;
                }
            }
        }

        // Add mean pressure modes; 
        for(i = 0; i < nel; ++i)
        {
            graphVertOffset[(ReorderedGraphVertId[1][AddMeanPressureToEdgeId[i]]+1)*nvel-1] += 1;
        }
        
        locToGloMap->SetSignChange(signChange);

        // Negate the vertices and edges with only a partial
        // Dirichlet conditon. Essentially we check to see if an edge
        // has a mixed Dirichlet with Neumann/Robin Condition and if
        // so negate the offset associated with this vertex. 

        
        map<int,int> DirVertChk;

        for(i = 0; i < bndConditionsVec[0].num_elements(); ++i)
        {
            cnt = 0;
            for(j = 0; j < nvel; ++j)
            {
                if(bndConditionsVec[j][i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    cnt ++;
                }
            }
            
            // Case where partial Dirichlet boundary condition
            if((cnt > 0)&&(cnt < nvel)) 
            {
                for(j  = 0; j < nvel; ++j)
                {
                    if(bndConditionsVec[j][i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {

                        //negate graph offsets which should be
                        //Dirichlet conditions
                        for(k = 0; k < bndCondExp[i]->GetExpSize(); ++k)
                        {
                            // vertices with mix condition;
                            id = bndCondExp[i]->GetExp(k)->GetGeom1D()->GetVid(0);
                            if(DirVertChk.count(id) == 0)
                            {
                                DirVertChk[id] = 1;
                                graphVertOffset[ReorderedGraphVertId[0][id]*nvel+j] *= -1; 
                            }
                            
                            id = bndCondExp[i]->GetExp(k)->GetGeom1D()->GetVid(1);
                            if(DirVertChk.count(id) == 0)
                            {
                                DirVertChk[id] = 1;
                                graphVertOffset[ReorderedGraphVertId[0][id]*nvel+j] *= -1; 
                            }
                            
                            // edges with mixed id; 
                            id = bndCondExp[i]->GetExp(k)->GetGeom1D()->GetEid();
                            graphVertOffset[ReorderedGraphVertId[1][id]*nvel+j] *= -1; 
                        }
                    }
                }
            }
        }
        
        cnt = 0;
        // assemble accumulative list of full Dirichlet values. 
        for(i = 0; i < firstNonDirGraphVertId*nvel; ++i)
        {
            diff = abs(graphVertOffset[i]); 
            graphVertOffset[i] = cnt; 
            cnt += diff;
        }
        
        // set Dirichlet values with negative values to Dirichlet value
        for(i = firstNonDirGraphVertId*nvel; i <  graphVertOffset.num_elements(); ++i)
        {
            if(graphVertOffset[i] < 0)
            {
                diff = -graphVertOffset[i]; 
                graphVertOffset[i] = -cnt; 
                cnt += diff;
            }
        }
    
        // cnt currently contains the number of global dir degrees of freedom. 
        locToGloMap->SetNumGlobalDirBndCoeffs(cnt);

        // Accumulate all interior degrees of freedom with positive
        // offset values
        for(i = firstNonDirGraphVertId*nvel; i < graphVertOffset.num_elements(); ++i)
        {
            if(graphVertOffset[i] >= 0)
            {
                diff = graphVertOffset[i]; 
                graphVertOffset[i] = cnt; 
                cnt += diff; 
            }
        }

        // Finally set negative entries (corresponding to Dirichlet
        // values ) to be positive
        for(i = firstNonDirGraphVertId*nvel; i < graphVertOffset.num_elements(); ++i)
        {
            if(graphVertOffset[i] < 0)
            {
                graphVertOffset[i] = -graphVertOffset[i]; 
            }
        }

        // Allocate the proper amount of space for the class-data and fill
        // information that is already known
        cnt = 0; 
        int numLocalBndCoeffs = 0;
        int numLocalCoeffs = 0;
        for(i = 0; i < nel; ++i)
        {
            numLocalBndCoeffs += nvel*locExpVector[i]->NumBndryCoeffs() + 1;
        }
        numLocalCoeffs = m_pressure->GetNcoeffs()-nel + numLocalBndCoeffs;
        
        locToGloMap->SetNumLocalCoeffs(numLocalCoeffs);
        locToGloMap->SetNumLocalBndCoeffs(numLocalBndCoeffs);

        Array<OneD, int> localToGlobalMap(numLocalCoeffs,-1);
        Array<OneD, int> localToGlobalBndMap(numLocalBndCoeffs,-1);
        Array<OneD, int> bndCondCoeffsToGlobalCoeffsMap(nLocBndCondDofs,-1);
    
        Array<OneD, NekDouble> localToGlobalSign;
        Array<OneD, NekDouble> localToGlobalBndSign;

        // If required, set up the sign-vector
        if(signChange)
        {
            localToGlobalSign = Array<OneD, NekDouble>(numLocalCoeffs,1.0);
            localToGlobalBndSign = Array<OneD, NekDouble>(numLocalBndCoeffs,1.0);
        }

        locToGloMap->SetGlobalSysSolnType(solnType);
        locToGloMap->SetStaticCondLevel(staticCondLevel);
        locToGloMap->SetNumPatches(nel);
        Array<OneD, unsigned int> numLocalBndCoeffsPerPatch(nel);
        Array<OneD, unsigned int> numLocalIntCoeffsPerPatch(nel);

        for(i = 0; i < nel; ++i)
        {
            numLocalBndCoeffsPerPatch[i] = (unsigned int) nvel*locExpVector[m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i)]->NumBndryCoeffs() + 1;
            numLocalIntCoeffsPerPatch[i] = (unsigned int) m_pressure->GetExp(m_pressure->GetOffset_Elmt_Id(i))->GetNcoeffs()-1;
        }
        
        /**
         * STEP 4: Now, all ingredients are ready to set up the actual
         * local to global mapping.
         *
         * The remainder of the map consists of the element-interior
         * degrees of freedom. This leads to the block-diagonal submatrix
         * as each element-interior mode is globally orthogonal to modes
         * in all other elements.
         */
        cnt = 0;
        int nv,velnbndry;
        Array<OneD, unsigned int> bmap;
        
        // Loop over all the elements in the domain in shuffled
        // ordering (element type consistency)
        for(i = 0; i < nel; ++i)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[eid]);

            velnbndry = locExpansion->NumBndryCoeffs();

            // require an inverse ordering of the bmap system to store
            // local numbering system which takes matrix these
            // matrices. Therefore get hold of elemental bmap and set
            // up an inverse map
            map<int,int> inv_bmap;
            locExpansion->GetBoundaryMap(bmap);
            for(j = 0; j < bmap.num_elements(); ++j)
            {
                inv_bmap[bmap[j]] = j;
            }
            
            // Loop over all edges (and vertices) of element i
            for(j = 0; j < locExpansion->GetNedges(); ++j)
            {
                nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                edgeOrient       = (locExpansion->GetGeom2D())->GetEorient(j);
                meshEdgeId       = (locExpansion->GetGeom2D())->GetEid(j);
                meshVertId       = (locExpansion->GetGeom2D())->GetVid(j);
                
                locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);
                // Set the global DOF for vertex j of element i
                for(nv = 0; nv < nvel; ++nv)
                {
                    localToGlobalMap[cnt+nv*velnbndry+inv_bmap[locExpansion->GetVertexMap(j)]] = graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nvel+ nv];
                    
                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        localToGlobalMap[cnt+nv*velnbndry+inv_bmap[edgeInteriorMap[k]]] =  graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nvel+nv]+k;
                    }
                }
                
                // Fill the sign vector if required
                if(signChange)
                {
                    for(nv = 0; nv < nvel; ++nv)
                    {
                        for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                        {
                            localToGlobalSign[cnt+nv*velnbndry + inv_bmap[edgeInteriorMap[k]]] = (NekDouble) edgeInteriorSign[k];
                        }
                    }
                }
            }

            // sort out mean pressure dof. 
            nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(AddMeanPressureToEdgeId[eid])-2;
            localToGlobalMap[cnt + nvel*velnbndry] = graphVertOffset[(ReorderedGraphVertId[1][AddMeanPressureToEdgeId[eid]]+1)*nvel-1]+nEdgeInteriorCoeffs + pressureEdgeOffset[AddMeanPressureToEdgeId[eid]];
            
            pressureEdgeOffset[AddMeanPressureToEdgeId[eid]] += 1;
            
            cnt += velnbndry*nvel+ m_pressure->GetExp(eid)->GetNcoeffs();
        }
       
        // Set up the mapping for the boundary conditions
        offset = cnt = 0;
        for(nv = 0; nv < nvel; ++nv)
        {
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndSegExp  = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));

                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                    for(k = 0; k < 2; k++)
                    {
                        meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                        bndCondCoeffsToGlobalCoeffsMap[cnt+bndSegExp->GetVertexMap(k)] = graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nvel+nv];
                    }
                    
                    meshEdgeId = (bndSegExp->GetGeom1D())->GetEid();
                    bndEdgeCnt = 0;
                    nEdgeCoeffs = bndSegExp->GetNcoeffs();
                    for(k = 0; k < nEdgeCoeffs; k++)
                    {
                        if(bndCondCoeffsToGlobalCoeffsMap[cnt+k] == -1)
                        {
                            bndCondCoeffsToGlobalCoeffsMap[cnt+k] =
                                graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nvel+nv]+bndEdgeCnt;
                            bndEdgeCnt++;
                        }
                    }
                }
                offset += bndCondExp[i]->GetNcoeffs();
            }
        }
        
        globalId = Vmath::Vmax(numLocalCoeffs,&localToGlobalMap[0],1)+1;
        locToGloMap->SetNumGlobalBndCoeffs(globalId);
        
        /**
         * STEP 5: The boundary condition mapping is generated from the
         * same vertex renumbering and fill in a unique interior map.
         */
        cnt=0;
        for(i = 0; i < numLocalCoeffs; ++i)
        {
            if(localToGlobalMap[i] == -1)
            {
                localToGlobalMap[i] = globalId++;
            }
            else
            {
                if(signChange)
                {
                    localToGlobalBndSign[cnt]=localToGlobalSign[i];
                }
                localToGlobalBndMap[cnt++]=localToGlobalMap[i];
            }
        }
        locToGloMap->SetNumGlobalCoeffs(globalId);
        locToGloMap->SetLocalToGlobalMap(localToGlobalMap);
        locToGloMap->SetLocalToGlobalBndMap(localToGlobalBndMap);
        locToGloMap->SetBndCondCoeffsToGlobalCoeffsMap(bndCondCoeffsToGlobalCoeffsMap);
        locToGloMap->SetNumLocalBndCoeffsPerPatch(numLocalBndCoeffsPerPatch);
        locToGloMap->SetNumLocalIntCoeffsPerPatch(numLocalIntCoeffsPerPatch);
        if (signChange)
        {
            locToGloMap->SetLocalToGlobalSign(localToGlobalSign);
            locToGloMap->SetLocalToGlobalBndSign(localToGlobalBndSign);
        }

        // Set up the local to global map for the next level when using
        // multi-level static condensation
        if( (solnType == MultiRegions::eDirectMultiLevelStaticCond) )
        {
            if(staticCondLevel < (bottomUpGraph->GetNlevels()-1))
            {
                MultiRegions::LocalToGlobalBaseMap *locToGloPointer = m_locToGloMap.get();

                Array<OneD, int> vwgts_perm(Dofs[0].size()+Dofs[1].size()-firstNonDirGraphVertId);

                for(i = 0; i < ReorderedGraphVertId[0].size(); ++i)
                {
                    if(ReorderedGraphVertId[0][i] >= firstNonDirGraphVertId)
                    {
                        vwgts_perm[ReorderedGraphVertId[0][i]-firstNonDirGraphVertId] = Dofs[0][i];
                    }
                }

                for(i = 0; i < ReorderedGraphVertId[1].size(); ++i)
                {
                    if(ReorderedGraphVertId[1][i] >= firstNonDirGraphVertId)
                    {
                        vwgts_perm[ReorderedGraphVertId[1][i]-firstNonDirGraphVertId] = Dofs[1][i];
                    }
                }
        
                bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                //bottomUpGraph->Dump();

                locToGloMap->SetNextLevelLocalToGlobalMap(MemoryManager<MultiRegions::LocalToGlobalBaseMap>::AllocateSharedPtr(locToGloPointer,bottomUpGraph));
            }
        }
    }

    void CoupledLinearNS::v_Output(void)
    {
        
        Array<OneD, Array<OneD, NekDouble > > fieldcoeffs(m_fields.num_elements()+1);
        Array<OneD, std::string> variables(m_fields.num_elements()+1);
        int i;
        
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
        }
        
        fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());
        
        // project pressure field to velocity space
        m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());
        
        m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);
        variables[i] = "p"; 
  
        ADRBase::Output(m_fields[0],fieldcoeffs,variables);
    }
}

/**
* $Log: CoupledLinearNS.cpp,v $
**/
