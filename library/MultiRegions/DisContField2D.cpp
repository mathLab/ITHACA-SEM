//////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.cpp
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
// Description: Field definition for 1D domain with boundary
// conditions using LDG-H
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {

        DisContField2D::DisContField2D(void)
        {
        }

        DisContField2D::DisContField2D(const DisContField2D &In):
            ExpList2D(In),
            m_bndConstraint(In.m_bndConstraint),
            m_bndTypes(In.m_bndTypes),
            m_globalBndMat(In.m_globalBndMat),
            m_trace(In.m_trace),  
            m_elmtToTrace(In.m_elmtToTrace),
            m_bndEidToTraceEid(In.m_bndEidToTraceEid),
            m_bndExpAdjacentOrient(In.m_bndExpAdjacentOrient),
            m_elmtTraceMap(In.m_elmtTraceMap),
            m_elmtTraceSign(In.m_elmtTraceSign)
        {
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const int bc_loc):
            ExpList2D(graph2D)     
        {
            GenerateFieldBnd2D(graph2D,bcs,bcs.GetVariable(bc_loc));

            // Set up Trace space
            m_trace = MemoryManager<GenExpList1D>::AllocateSharedPtr(m_bndConstraint, m_bndTypes,*m_exp,graph2D);
            SetUpTraceMappings(graph2D);
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const std::string variable):
            ExpList2D(graph2D)  
        {
            GenerateFieldBnd2D(graph2D,bcs,variable);

            // Set up Trace space
            m_trace = MemoryManager<GenExpList1D>::AllocateSharedPtr(m_bndConstraint,m_bndTypes,*m_exp,graph2D);
            SetUpTraceMappings(graph2D);
        }



        void DisContField2D::GenerateFieldBnd2D(SpatialDomains::MeshGraph2D &graph2D,
                                                SpatialDomains::BoundaryConditions &bcs, 
                                                const std::string variable)
        {
            int mycnt = 0;
            int i,j,m;
            int cnt = 0;
            int coeffcnt = 0;
            int physcnt = 0;
            int nbnd;
            int nummodes;
            int npoints;
            SpatialDomains::BoundaryConditionShPtr locBCond;
            SpatialDomains::BoundaryRegionShPtr    locBregion;
            LocalRegions  ::SegExpSharedPtr        collSeg;
            StdRegions    ::StdSegExpSharedPtr     collStdSeg;
            SpatialDomains::SegGeomSharedPtr       SegGeom;

            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();    

            MultiRegions  ::ExpList1DSharedPtr  locExpList;   
            LocalRegions  ::SegExpSharedPtr     locSegExp;
            SpatialDomains::Geometry2DSharedPtr geom2D;
            StdRegions    ::StdExpMap           vmap;
            
            nbnd = bregions.size();

            m_bndConstraint = Array<OneD,MultiRegions::ExpList1DSharedPtr>(nbnd);
            m_bndTypes = Array<OneD,SpatialDomains::BoundaryConditionType>(nbnd);
            
            // Set up matrix map
            m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {  
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    locBregion = bregions[i];

                    locExpList = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*locBregion,graph2D);

                    coeffcnt=0;
                    physcnt=0;
                    for(j = 0; j < locExpList->GetExpSize(); j++)
                    {          
                        locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(locExpList->GetExp(j));

                        nummodes = locSegExp->GetNcoeffs();
                        npoints  = locSegExp->GetTotPoints();
                                
                        // Create a new BasisKey for projecting the
                        // dirichlet boundary conditions onto the
                        // boundary.  If the original BasisKey has N
                        // modes, the new BasisKey employs N GLL
                        // quadrature points such that the FwdTrans
                        // using this new basis is in fact a
                        // collocation projection trough this
                        // GLL-points. As a result, the expansion will
                        // be C0 continuous on the boundary.
                        
                        // The PointsKey used for the (collocation)
                        // projection
                        LibUtilities::PointsKey collPointsKey(nummodes,LibUtilities::eGaussLobattoLegendre);  
                        // The BasisKey used for the (collocation)
                        // projection
                        LibUtilities::BasisKey collBasisKey(locSegExp->GetBasisType(0),nummodes,collPointsKey);
                        
                        // Create a segment based on the new BasisKey
                        // in order to perfrom the projection
                        collSeg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(collBasisKey, locSegExp->GetGeom1D());
                        collStdSeg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(collBasisKey);
                                
                        // Calculate the coordinates of the N GLL
                        // quadrature points
                        Array<OneD,NekDouble> x0(nummodes,0.0);
                        Array<OneD,NekDouble> x1(nummodes,0.0);
                        Array<OneD,NekDouble> x2(nummodes,0.0);                         
                        collSeg->GetCoords(x0,x1,x2);
                                
                        // Evaluate the Dirichlet boundary condition
                        // at the N GLL quadrature points
                        for(m = 0; m < nummodes; ++m)
                        {
                            (collSeg->UpdatePhys())[m] = boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(locBCond)->
                                m_DirichletCondition.Evaluate(x0[m],x1[m],x2[m]);
                        }
                        // Perform a FwdTrans() to calculate the
                        // expansion coefficients.  As both the
                        // original as the new expansion (used for the
                        // projection) are of the same order, the
                        // coefficients will be identical and hence,
                        // can directly be stored at the right place
                        Array<OneD,NekDouble> outarray;
                        Array<OneD,NekDouble> inarray;
                        collStdSeg->FwdTrans(collSeg->GetPhys(),outarray = (locExpList->UpdateCoeffs()) + coeffcnt);   
                        // also update physical points 
                        locSegExp->BwdTrans( inarray = (locExpList->GetCoeffs()) + coeffcnt, outarray = locExpList->UpdatePhys() + physcnt);
                    
                        coeffcnt += nummodes;
                        physcnt += npoints;
                    }
                                
                    m_bndConstraint[cnt] = locExpList;
                    m_bndTypes[cnt++] = SpatialDomains::eDirichlet;                     
                } // end if Dirichlet
            }
            // list other boundaries
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];

                if(locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                    locBregion = bregions[i];
                    
                    locExpList = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*locBregion,graph2D);

                    npoints = locExpList->GetPointsTot();
                    Array<OneD,NekDouble> x0(npoints,0.0);
                    Array<OneD,NekDouble> x1(npoints,0.0);
                    Array<OneD,NekDouble> x2(npoints,0.0);                         
                    locExpList->GetCoords(x0,x1,x2);     

                    if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        for(m = 0; m < npoints; m++)
                        {
                            (locExpList->UpdatePhys())[m] = boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(locBCond)->
                                m_NeumannCondition.Evaluate(x0[m],x1[m],x2[m]);
                        }
                        locExpList->IProductWRTBase(*locExpList); 
                        m_bndConstraint[cnt] = locExpList;
                        m_bndTypes[cnt++] = SpatialDomains::eNeumann;      
                    }
                    else if(locBCond->GetBoundaryConditionType() == SpatialDomains::eRobin)
                    {        
                        boost::shared_ptr<SpatialDomains::RobinBoundaryCondition> robinBC = 
                            boost::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(locBCond);
                        for(m = 0; m < npoints; m++)
                        {
                            (locExpList->UpdatePhys())[m] = -robinBC->m_a.Evaluate(x0[m],x1[m],x2[m])/
                                robinBC->m_b.Evaluate(x0[m],x1[m],x2[m]);
                        }
                        locExpList->IProductWRTBase(*locExpList); 
                        m_bndConstraint[cnt] = locExpList;
                        m_bndTypes[cnt++] = SpatialDomains::eRobin;      
                    }
                    else
                    {
                        ASSERTL0(false,"This type of BC not implemented yet");
                    }                    
                } 
            }
        }

        void DisContField2D::SetUpTraceMappings(SpatialDomains::MeshGraph2D &graph2D)
        {
            int i,j,k,cnt,id, order_e,gid;
            int ntrace_exp = m_trace->GetExpSize();
            int nel = m_exp->size();
            LocalRegions::SegExpSharedPtr  locSegExp;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr  locTriExp;
            //SpatialDomains::SegGeomSharedPtr SegGeom;

            SpatialDomains::Geometry1DSharedPtr SegGeom;
            
            Array<OneD,int> MeshEdgeNo(graph2D.GetNseggeoms(),-1);

            // determine mapping from edge expansion to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_trace->GetExp(i)))
                {
                    MeshEdgeNo[(locSegExp->GetGeom1D())->GetEid()] = i;
                }
                else
                {
                    ASSERTL0(false,"Dynamics cast to segment expansion failed");
                }
            }

            // Count edges
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += (*m_exp)[i]->GetNedges();
            }
            

            Array<OneD, LocalRegions::SegExpSharedPtr> edgemap(cnt);
            m_elmtToTrace = Array<OneD, Array<OneD,LocalRegions::SegExpSharedPtr> >(nel);

            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = edgemap + cnt; 

                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*m_exp)[i]))
                {
                    for(j = 0; j < locQuadExp->GetNedges(); ++j)
                    {   
                        SegGeom = (locQuadExp->GetGeom2D())->GetEdge(j);
                        
                        id = SegGeom->GetEid();

                        if(MeshEdgeNo[id] != -1)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::SegExp> ((*m_trace).GetExp(MeshEdgeNo[id]));
                        }
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>((*m_exp)[i]))
                {
                    for(j = 0; j < locTriExp->GetNedges(); ++j)
                    {    
                        SegGeom = (locTriExp->GetGeom2D())->GetEdge(j);

                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeNo[id] != -1)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::SegExp> ((*m_trace).GetExp(MeshEdgeNo[id]));
                        }
                    }
                
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }                    
                cnt += (*m_exp)[i]->GetNedges();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < m_bndConstraint.num_elements(); ++i)
            {
                cnt += m_bndConstraint[i]->GetExpSize();
            }

            Array<OneD, int> bedgemap(cnt);
            m_bndEidToTraceEid = Array<OneD, Array<OneD, int> >(m_bndConstraint.num_elements());
            
            Array<OneD, AdjacentEdgeOrientation> bedgeorient(cnt);
            m_bndExpAdjacentOrient = Array<OneD, Array<OneD, AdjacentEdgeOrientation> > (m_bndConstraint.num_elements());
            

            cnt = 0;
            m_numTraceDirichletBCs = 0;
            m_numTraceDirichletPhysBCs = 0;
            for(i = 0; i < m_bndConstraint.num_elements(); ++i)
            {
                m_bndEidToTraceEid[i] = bedgemap + cnt; 
                m_bndExpAdjacentOrient[i] = bedgeorient + cnt;

                for(j = 0; j < m_bndConstraint[i]->GetExpSize(); ++j)
                {

                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_bndConstraint[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        
                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeNo[id] != -1)
                        {
                            (m_bndEidToTraceEid[i])[j] = MeshEdgeNo[id];
                        }

                        // Check to see which way boundar edge is
                        // orientated with respect to connecting
                        // element counter-clockwise convention.

                        SpatialDomains::SegGeomSharedPtr Sgeom;
                        if(!(Sgeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(SegGeom)))
                        {
                            ASSERTL0(false,"dynamic cast to a SegGeom failed"); 
                        }
                        SpatialDomains::ElementEdgeVectorSharedPtr con_elmt
                            = graph2D.GetElementsFromEdge(Sgeom);
                        
                        if((boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*con_elmt)[0]->m_Element))->GetEorient((*con_elmt)[0]->m_EdgeIndx) == StdRegions::eForwards)
                        {
                            m_bndExpAdjacentOrient[i][j] = eAdjacentEdgeIsForwards;
                        }
                        else
                        {
                            m_bndExpAdjacentOrient[i][j] = eAdjacentEdgeIsBackwards;
                        }
                        
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local Segment expansion failed");
                    }
                    
                    m_numTraceDirichletBCs     += locSegExp->GetNcoeffs();
                    m_numTraceDirichletPhysBCs += locSegExp->GetTotPoints();
                }
                cnt += j;
            }
            
            // Set up integer mapping array and sign change for each
            // degree of freedom

            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                nbndry += (*m_exp)[i]->NumDGBndryCoeffs();
            }

            m_elmtTraceMap  = Array<OneD, int > (nbndry);
            m_elmtTraceSign = Array<OneD, int > (nbndry,1);

            nbndry = cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                nbndry += (*m_exp)[i]->NumDGBndryCoeffs();

                for(j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                {   
                    locSegExp = m_elmtToTrace[i][j];
                    SegGeom = locSegExp->GetGeom1D();
                    
                    id  = SegGeom->GetEid();
                    gid = m_trace->GetCoeff_Offset(MeshEdgeNo[id]);
                    
                    order_e = locSegExp->GetNcoeffs();
                    
                    if((*m_exp)[i]->GetEorient(j) == StdRegions::eForwards)
                    {
                        for(k = 0; k < order_e; ++k)
                        {
                            m_elmtTraceMap[k+cnt] = gid + k;
                        }
                    }
                    else // backwards orientated
                    {
                        switch(locSegExp->GetBasisType(0))
                        {
                        case LibUtilities::eModified_A:
                            // reverse vertex order
                            m_elmtTraceMap[cnt] = gid + 1;
                            m_elmtTraceMap[cnt+1] = gid;
                            for(k = 2; k < order_e; ++k)
                            {
                                m_elmtTraceMap[k+cnt] = gid + k;
                            }

                            // negate odd modes
                            for(k = 3; k < order_e; k+=2)
                            {
                                m_elmtTraceSign[cnt+k] = -1;
                            }
                            
                            
                            break;
                        case LibUtilities::eGLL_Lagrange:
                            // reverse  order
                            for(k = 0; k < order_e; ++k)
                            {
                                m_elmtTraceMap[cnt+order_e-k-1] = gid + k;
                            }
                            break;
                        default:
                            ASSERTL0(false,"Boundary type not permitted");
                            
                        }
                    }
                    cnt += order_e;
                }
            }


        }
        
        
        DisContField2D::~DisContField2D()
        {
        }

        GlobalLinSysSharedPtr DisContField2D::GenGlobalBndLinSys(const GlobalLinSysKey &mkey)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt;
            StdRegions::StdExpansionVectorIter def;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;

            int totDofs       = m_trace->GetNcoeffs();
            int NumDirBCs     = m_numTraceDirichletBCs;
            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero    = 0.0,sign1,sign2; 
            NekDouble factor1 = mkey.GetFactor1();
            NekDouble factor2 = mkey.GetFactor2();
            StdRegions::MatrixType linsystype = mkey.GetLinSysType();
            
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);            
            ASSERTL0(linsystype == StdRegions::eUnifiedDGHelmBndLam,
                     "Routine currently only tested for UnifiedDGHelmholtz");
            
            // fill global matrix 
            cnt = 0;
            for(n = 0; n < (*m_exp).size(); ++n)
            {
                // Matrix to Bnd Sys
                LocalRegions::MatrixKey Umatkey(linsystype, (*m_exp)[n]->DetExpansionType(),*((*m_exp)[n]), factor1,factor2);
                DNekScalMat &BndSys = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 
                
                loc_lda = BndSys.GetColumns();
                
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = m_elmtTraceMap[cnt+i]; 
                    sign1 = m_elmtTraceSign[cnt+i];

                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = m_elmtTraceMap [cnt+j]; 
                            sign2 = m_elmtTraceSign[cnt+j];

                            if(gid2 >= NumDirBCs)
                            {

                                (*Gmat)(gid1 - NumDirBCs,gid2 - NumDirBCs) 
                                    += sign1*sign2*(BndSys)(i,j);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            
            // Believe that we need a call of the type:
            //linsys=MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,eWrapper);
            linsys       = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }


        GlobalLinSysSharedPtr DisContField2D::GetGlobalBndLinSys(const GlobalLinSysKey &mkey) 
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalBndMat->find(mkey);

            if(matrixIter == m_globalBndMat->end())
            {
                glo_matrix = GenGlobalBndLinSys(mkey);
                (*m_globalBndMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

        void DisContField2D::HelmSolve(DisContField2D &Fce, NekDouble lambda)
        {
            int e,i,j,n,cnt,cnt1,nbndry, order_e;
            int nexp = GetExpSize();
            static DNekScalBlkMatSharedPtr InvUDGHelm;
            LocalRegions::SegExpSharedPtr SegExp;
            StdRegions::StdExpansionSharedPtr BndExp;

            Array<OneD,NekDouble> f(m_ncoeffs);
            Array<OneD,NekDouble> e_f, e_bndsol;
            Array<OneD,NekDouble> e_lambda;
            NekDouble tau = 10;

            //----------------------------------
            // Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(Fce.GetPhys(),f);
            Vmath::Neg(m_ncoeffs,f,1);

            //----------------------------------
            // Solve continuous flux System
            //----------------------------------
            int GloBndDofs = m_trace->GetNcoeffs();
            int NumDirichlet = m_numTraceDirichletBCs;
            int e_ncoeffs, loc,id,offset;
            NekDouble sign;

            Array<OneD,NekDouble> BndSol = m_trace->UpdateCoeffs(); // linked data
            DNekVec Floc;

            // Zero trace space
            Vmath::Zero(GloBndDofs,BndSol,1);

            cnt = 0;
            for(n = 0; n < nexp; ++n)
            {
                cnt = max(cnt,(*m_exp)[n]->NumDGBndryCoeffs());
            }
            e_lambda = Array<OneD,NekDouble>(cnt);

#if 1       //Forcing based on analytical derivation From Cockburn
            for(cnt1 = cnt = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs(); 		    

                LocalRegions::MatrixKey Umatkey(StdRegions::eUnifiedDGLamToU, 
                   (*m_exp)[n]->DetExpansionType(),*((*m_exp)[n]), lambda, tau);

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f+cnt;
                
                DNekVec     ElmtFce(e_ncoeffs, e_f ,eWrapper);
                DNekScalMat &LamToU = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 

                Floc = Transpose(LamToU)*ElmtFce;

                // Put value into trace space expansion
                for(i = 0; i < nbndry; ++i)
                {
                    id = m_elmtTraceMap[cnt1+i];
                    sign = m_elmtTraceSign[cnt1+i];
                    m_trace->SetCoeffs(id,sign*Floc[i] + m_trace->GetCoeffs(id));
                }
                cnt  += e_ncoeffs;
                cnt1 += nbndry;
            }

#else      // Forcing based on flux derivative. 
            Array<OneD,NekDouble> ftmp;

            // Forcing  using  Direct form <\phi,[[qf ]]>
            offset = 0;
            for(cnt = e = 0; e < nexp; ++e)
            {
                nbndry = (*m_exp)[e]->NumDGBndryCoeffs();
                e_ncoeffs = (*m_exp)[e]->GetNcoeffs();

                // Get BndSysForce Matrix
                LocalRegions::MatrixKey Bmatkey(StdRegions::eUnifiedDGHelmBndFce,  (*m_exp)[e]->DetExpansionType(),*((*m_exp)[e]), lambda, tau);
                
                DNekScalMat &BndSys = *((*m_exp)[e]->GetLocMatrix(Bmatkey));                               
                DNekVec vin (e_ncoeffs,ftmp = f+offset,eWrapper);
                DNekVec vout(nbndry);
                
                vout = BndSys*vin;

                // Subtract vout from forcing terms 
                for(i = 0; i < nbndry; ++i)
                {
                    id = m_elmtTraceMap[cnt];
                    sign = m_elmtTraceSign[cnt];
                    m_trace->SetCoeff(id, m_trace->GetCoeff(id) - sign*vout[i]);
                    cnt++;
                }

                offset += e_ncoeffs;
            }
#endif

            // Copy Dirichlet boundary conditions into trace space        
            for(i = 0; i < m_bndTypes.num_elements(); ++i)
            {
                if(m_bndTypes[i] == SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndConstraint[i]->GetExpSize(); ++e)
                    {
                        id     = m_bndEidToTraceEid[i][e];
                        offset = m_trace->GetCoeff_Offset(id);
                        m_bndConstraint[i]->PutCoeffsInToElmtExp();
                        BndExp = m_bndConstraint[i]->GetExp(e);
                        for(j = 0; j < BndExp->GetNcoeffs(); ++j)
                        {
                            m_trace->SetCoeff(offset + j,BndExp->GetCoeff(j));
                        }
                    }
                }
            }


            // Dirichlet boundary forcing 
            for(cnt = e = 0; e < nexp; ++e)
            {
                nbndry = (*m_exp)[e]->NumDGBndryCoeffs();
                // check to see if element has Dirichlet boundary
                // Probably could use a quicker check here
                if(Vmath::Vmin(nbndry,&m_elmtTraceMap[cnt],1) < NumDirichlet)
                {
                    // Get BndSys Matrix

                    LocalRegions::MatrixKey Bmatkey(StdRegions::eUnifiedDGHelmBndLam,  (*m_exp)[e]->DetExpansionType(),*((*m_exp)[e]), lambda, tau);
                    
                    DNekScalMat &BndSys = *((*m_exp)[e]->GetLocMatrix(Bmatkey));               
                    DNekVec vin (nbndry);
                    DNekVec vout(nbndry);


                    // Set up Edge Dirichlet Values
                    for(cnt1 = i = 0; i < (*m_exp)[e]->GetNedges(); ++i)
                    {
                        e_ncoeffs = (*m_exp)[e]->GetEdgeNcoeffs(i); 

                        id = m_elmtTraceMap[cnt+cnt1];
                        
                        if(id < NumDirichlet)
                        {
                            for(j = 0; j < e_ncoeffs; ++j)
                            {
                                id = m_elmtTraceMap[cnt+cnt1];
                                sign = m_elmtTraceSign[cnt+cnt1];
                                vin[cnt1] = sign*m_trace->GetCoeffs(id); 
                                cnt1++;
                            } 
                        }
                        else
                        {
                            for(j = 0; j < e_ncoeffs; ++j)
                            {
                                vin[cnt1] = 0.0;
                                cnt1++;
                            }
                        }
                    }

                    vout = BndSys*vin;
                    
                    // Subtract vout from forcing terms 
                    for(i = 0; i < nbndry; ++i)
                    {
                        id = m_elmtTraceMap[cnt+i];
                        sign = m_elmtTraceSign[cnt+i];

                        if(id >= NumDirichlet)
                        {
                            m_trace->SetCoeff(id, m_trace->GetCoeff(id)
                                              - sign*vout[i]);
                        }
                    }
                }
                cnt += nbndry;
            }
            
            if(GloBndDofs - NumDirichlet > 0)
            {
                GlobalLinSysKey       key(StdRegions::eUnifiedDGHelmBndLam,
                                          lambda,tau,eDirectFullMatrix);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                
                Array<OneD,NekDouble> sln = BndSol+NumDirichlet;

                LinSys->Solve(sln,sln);
            }
            
            //----------------------------------
            // Setup forcing for local interior solves
            //----------------------------------
            Vmath::Zero(m_ncoeffs,m_coeffs,1);
            
            for(cnt = cnt1 = i = 0; i < nexp; ++i)
            {
                nbndry = (*m_exp)[i]->NumDGBndryCoeffs();
                e_f      = f + cnt;
                // put elemental solutions into local space; 
                
                for(j = 0; j < nbndry; ++j)
                {
                    e_lambda[j] = m_elmtTraceSign[cnt1+j]*m_trace->GetCoeffs(m_elmtTraceMap[cnt1+j]);
                }

                (*m_exp)[i]->AddUDGHelmholtzTraceTerms(tau, e_lambda, 
                                                       m_elmtToTrace[i], 
                                                       e_f);
                cnt  += (*m_exp)[i]->GetNcoeffs();
                cnt1 += nbndry;
            }
            
            // Inverse block diagonal interior solve
            if(!InvUDGHelm.get())
            {
                InvUDGHelm = SetupBlockMatrix(StdRegions::eInvUnifiedDGHelmholtz, 
                                              lambda, tau);
            }
            
            DNekVec in (m_ncoeffs,f,eWrapper);
            DNekVec out(m_ncoeffs,m_coeffs,eWrapper);            
            out = (*InvUDGHelm)*in;            
        }

        // Construct the two trace vectors of the inner and outer
        // trace solution from the field contained in m_phys, where
        // the Weak dirichlet boundary conditions are listed in the
        // outer part of the vecotr

        void DisContField2D::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                                    Array<OneD,NekDouble> &Bwd)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            StdRegions::EdgeOrientation edgedir;
            int nquad_e,cnt,n,e,npts,offset, phys_offset;
            Array<OneD,NekDouble> e_tmp, e_tmp1;

            // zero vectors; 
            Vmath::Zero(Fwd.num_elements(),Fwd,1);
            Vmath::Zero(Bwd.num_elements(),Bwd,1);

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eForwards)
                    {
                        offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e, e_tmp = m_phys + phys_offset, 
                                                     e_tmp1 = Fwd + offset);
                    }
                }
            }
            
            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eBackwards)
                    {
                        offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e,e_tmp = m_phys + phys_offset, 
                                                     e_tmp1 = Bwd+offset);
                    }
                }
            }
            

            // fill boundary conditions into missing elements
            
            int id1,id2;
            for(n = 0; n < m_bndConstraint.num_elements(); ++n)
            {
                
                for(e = 0; e < m_bndConstraint[n]->GetExpSize(); ++e)
                {
                    npts = m_bndConstraint[n]->GetExp(e)->GetNumPoints(0);

                    if(m_bndExpAdjacentOrient[n][e] == eAdjacentEdgeIsForwards)
                    {
                        id1 = m_bndConstraint[n]->GetPhys_Offset(e) ;
                        id2 = m_trace->GetPhys_Offset(m_bndEidToTraceEid[n][e]);
                        Vmath::Vcopy(npts,&(m_bndConstraint[n]->GetPhys())[id1],1,
                                     &Fwd[id2],1);
                            
                    }
                    else
                    {
                        id1 = m_bndConstraint[n]->GetPhys_Offset(e) ;
                        id2 = m_trace->GetPhys_Offset(m_bndEidToTraceEid[n][e]);
                        Vmath::Vcopy(npts,&(m_bndConstraint[n]->GetPhys())[id1],1,
                                     &Bwd[id2],1);
                    }
                }
            }
        }
        
        void DisContField2D::ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int nquad_e,cnt,n,e,npts,offset,phys_offset;
            Array<OneD,NekDouble> e_tmp,e_tmp1;

            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            ASSERTL1(outarray.num_elements() >= m_trace->GetNpoints(),
                     "input array is of insufficient length");

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->GetEdgePhysVals(e,  m_elmtToTrace[n][e], 
                                                 e_tmp = m_phys + phys_offset,
                                                 e_tmp1 = outarray + offset);
                }
            }
        }
            

        /// Note this routine changes m_trace->m_coeffs space; 
        void DisContField2D::AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                                              Array<OneD, const NekDouble> &Fy, 
                                              Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray,e_Fx,e_Fy;

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryInt(e,m_elmtToTrace[n][e],
                                                        e_Fx = Fx + t_offset,
                                                        e_Fy = Fy + t_offset,
                                                        e_outarray = outarray+offset);
                }    
            }
        }

    } // end of namespace
} //end of namespace
