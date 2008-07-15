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
            m_elmtTraceMap(In.m_elmtTraceMap),
            m_elmtTraceSign(In.m_elmtTraceSign),
            m_traceInnerNormals(In.m_traceInnerNormals)
        {
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const int bc_loc):
            ExpList2D(graph2D)     
        {
            GenerateFieldBnd2D(graph2D,bcs,bcs.GetVariable(bc_loc));

            // Set up Trace space
            m_trace = MemoryManager<ExpList1D>::AllocateSharedPtr(m_bndConstraint, m_bndTypes,*m_exp,graph2D);
            SetUpTraceMappings(graph2D);
            SetUpTraceInnerNormals();
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const std::string variable):
            ExpList2D(graph2D)  
        {
            GenerateFieldBnd2D(graph2D,bcs,variable);

            // Set up Trace space
            m_trace = MemoryManager<ExpList1D>::AllocateSharedPtr(m_bndConstraint,m_bndTypes,*m_exp,graph2D);
            SetUpTraceMappings(graph2D);
            SetUpTraceInnerNormals();
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
                        collSeg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(collBasisKey, locSegExp->GetGeom());
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

            m_elmtTraceMap  = Array<OneD, Array<OneD,int > > (ntrace_exp);
            m_elmtTraceSign = Array<OneD, Array<OneD,int > > (ntrace_exp);
            
            // determine mapping from edge expansion to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_trace->GetExp(i)))
                {
                    MeshEdgeNo[(locSegExp->GetGeom())->GetEid()] = i;
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
                        SegGeom = (locQuadExp->GetGeom())->GetEdge(j);
                        
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
                        SegGeom = (locTriExp->GetGeom())->GetEdge(j);

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
                                                                        
            cnt = 0;
            m_numTraceDirichletBCs = 0;
            m_numTraceDirichletPhysBCs = 0;
            for(i = 0; i < m_bndConstraint.num_elements(); ++i)
            {

                m_bndEidToTraceEid[i] = bedgemap + cnt; 

                for(j = 0; j < m_bndConstraint[i]->GetExpSize(); ++j)
                {

                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_bndConstraint[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom();
                        
                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeNo[id] != -1)
                        {
                            (m_bndEidToTraceEid[i])[j] = MeshEdgeNo[id];
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
            
            m_elmtTraceMap  = Array<OneD, Array<OneD,int > > (nel);
            m_elmtTraceSign = Array<OneD, Array<OneD,int > > (nel);
            
            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                nbndry += (*m_exp)[i]->NumDGBndryCoeffs();
            }
            Array<OneD, int>  Map(nbndry);
            Array<OneD, int>  Sign(nbndry,1);
            
            nbndry = cnt = 0;
            for(i = 0; i < nel; ++i)
            {
#if 1  
                m_elmtTraceMap [i] = Map  + nbndry;
                m_elmtTraceSign[i] = Sign + nbndry; 
#else
                m_elmtTraceMap [i] = Array<OneD, int> ((*m_exp)[i]->NumDGBndryCoeffs()); 
                m_elmtTraceSign[i] = Array<OneD, int> ((*m_exp)[i]->NumDGBndryCoeffs(),1); 
#endif

                nbndry += (*m_exp)[i]->NumDGBndryCoeffs();
                cnt = 0;
                for(j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                {   
                    locSegExp = m_elmtToTrace[i][j];
                    SegGeom = locSegExp->GetGeom();
                    
                    id  = SegGeom->GetEid();
                    gid = m_trace->GetCoeff_Offset(MeshEdgeNo[id]);
                    
                    order_e = locSegExp->GetNcoeffs();
                    
                    if((*m_exp)[i]->GetCartesianEorient(j) == StdRegions::eForwards)
                    {
                        for(k = 0; k < order_e; ++k)
                        {
                            m_elmtTraceMap[i][k+cnt] = gid + k;
                        }
                    }
                    else // backwards orientated
                    {
                        switch(locSegExp->GetBasisType(0))
                        {
                        case LibUtilities::eModified_A:
                            // reverse vertex order
                            m_elmtTraceMap[i][cnt] = gid + 1;
                            m_elmtTraceMap[i][cnt+1] = gid;
                            for(k = 2; k < order_e; ++k)
                            {
                                m_elmtTraceMap[i][k+cnt] = gid + k;
                            }
                            
                            // negate odd modes
                            for(k = 3; k < order_e; k+=2)
                            {
                                m_elmtTraceSign[i][k+cnt] = -1;
                            }
                            
                            break;
                        case LibUtilities::eGLL_Lagrange:
                            // reverse  order
                            for(k = 0; k < order_e; ++k)
                            {
                                m_elmtTraceMap[i][cnt+order_e-k-1] = gid + k;
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
        
        void DisContField2D::SetUpTraceInnerNormals(void)
        {
            int i,n,e,offset;
            int coordim = GetCoordim(0); // assume all element have the same coordim.
            int trace_npts = m_trace->GetNpoints();
            int npoints_e,nexp,cnt;
            Array<TwoD, const NekDouble> normals; 
            SpatialDomains::GeomType Gtype;
            Array<OneD, bool> normal_set;
            
            // Declare storage
            m_traceInnerNormals = MemoryManager<Array<OneD,Array<OneD,NekDouble> > >::AllocateSharedPtr(coordim);
            
            nexp = GetExpSize();

            // count mesh edge
            cnt = 0;
            for(n = 0; n < nexp; ++n)
            {
                cnt += (*m_exp)[n]->GetNedges();
            }
            normal_set = Array<OneD,bool> (cnt,false);
            
            for(i = 0; i < coordim; ++i)
            {
                (*m_traceInnerNormals)[i] = Array<OneD,NekDouble>(trace_npts);
                
            }

            //Fill normals
            cnt = 0;
            for(n = 0; n < GetExpSize(); ++n)
            {
                Gtype   = (*m_exp)[n]->GetMetricInfo()->GetGtype();
                normals = (*m_exp)[n]->GetMetricInfo()->GetNormals();

                if(Gtype == SpatialDomains::eDeformed)
                {
                    for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e, ++cnt)
                    {
                        npoints_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                        
                        if((*m_exp)[n]->GetEorient(e) == StdRegions::eForwards)
                        {
                            normal_set[cnt] = true;
                            offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                            for(i = 0 ; i < coordim; ++i)
                            {
                                Vmath::Vcopy(npoints_e,&normals[e][i],1,
                                            &((*m_traceInnerNormals)[i])[offset],1);
                            }

                        }
                            
                    }
                    
                }
                else
                {
                    for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e, ++cnt)
                    {
                        npoints_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                        
                        if((*m_exp)[n]->GetEorient(e) == StdRegions::eForwards)
                        { 
                            normal_set[cnt] = true;
                            offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                           for(i = 0 ; i < coordim; ++i)
                            {
                                Vmath::Fill(npoints_e,normals[e][i],
                                            &((*m_traceInnerNormals)[i])[offset],1);
                            }
                        }                            
                    }
                }
            }

            // finally loop through and set any unset boundaries that
            // are on boundary but orientated backwards
            cnt = 0;
            for(n = 0; n < GetExpSize(); ++n)
            {
                Gtype   = (*m_exp)[n]->GetMetricInfo()->GetGtype();
                normals = (*m_exp)[n]->GetMetricInfo()->GetNormals();

                if(Gtype == SpatialDomains::eDeformed)
                {
                    for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e, ++cnt)
                    {
                        npoints_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                        
                        if(((*m_exp)[n]->GetEorient(e) == StdRegions::eBackwards) && (normal_set[cnt] == false))
                        {
                            normal_set[cnt] = true;
                            offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                            for(i = 0 ; i < coordim; ++i)
                            {
                                Vmath::Vcopy(npoints_e,&normals[e][i],1,
                                            &((*m_traceInnerNormals)[i])[offset],1);
                            }

                        }
                            
                    }
                }
                else
                {
                    for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e, ++cnt)
                    {
                        npoints_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                        
                        if(((*m_exp)[n]->GetEorient(e) == StdRegions::eBackwards)&&(normal_set[cnt] == false))
                        { 
                            normal_set[cnt] = true;
                            offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                           for(i = 0 ; i < coordim; ++i)
                            {
                                Vmath::Fill(npoints_e,normals[e][i],
                                            &((*m_traceInnerNormals)[i])[offset],1);
                            }
                        }                            
                    }
                }
            }
            
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            // make sure all normals are set
            for(i = 0; i < normal_set.num_elements(); ++i)
            {
                ASSERTL1(normal_set[i] == true,"Normal not set");
            }
#endif
            
        }
        
        DisContField2D::~DisContField2D()
        {
        }

        GlobalLinSysSharedPtr DisContField2D::GenGlobalBndLinSys(const GlobalLinSysKey &mkey)
	{
            int i,j,n,gid1,gid2,loc_lda;
            StdRegions::StdExpansionVectorIter def;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;

            int totDofs    = m_trace->GetNcoeffs();
            int NumDirBCs     = m_numTraceDirichletBCs;
            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero = 0.0,sign1,sign2; 
            NekDouble factor1 = mkey.GetFactor1();
            NekDouble factor2 = mkey.GetFactor2();
            StdRegions::MatrixType linsystype = mkey.GetLinSysType();
            
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);            
            ASSERTL0(linsystype == StdRegions::eUnifiedDGHelmBndSys,
                     "Routine currently only tested for UnifiedDGHelmholtz");
            
            // fill global matrix 
            for(n = 0; n < (*m_exp).size(); ++n)
            {
                // Matrix to Bnd Sys
                LocalRegions::MatrixKey Umatkey(linsystype, (*m_exp)[n]->DetExpansionType(),*((*m_exp)[n]), factor1,factor2);
                DNekScalMat &BndSys = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 
                
                cout << "BndSys" << endl;
                cout << BndSys << endl;
                
                loc_lda = BndSys.GetColumns();
                
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = m_elmtTraceMap [n][i]; 
                    sign1 = m_elmtTraceSign[n][i];
                    
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            
                            gid2  = m_elmtTraceMap [n][j]; 
                            sign2 = m_elmtTraceSign[n][j];
                            
                            if(gid2 >= NumDirBCs)
                            {

                                (*Gmat)(gid1 - NumDirBCs,gid2 - NumDirBCs) 
                                    += sign1*sign2*(BndSys)(i,j);
                            }
                        }
                    }
                }
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

#if 0       //Forcing based on analytical derivation  From Cockburn paper -- Slow convergence
            // Set up boundary space forcing
            for(cnt = n = 0; n < nexp; ++n)
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
                    id = m_elmtTraceMap[n][i];
                    sign = m_elmtTraceSign[n][i];
                    m_trace->SetCoeffs(id,sign*Floc[i] + m_trace->GetCoeffs(id));
                }
                cnt  += e_ncoeffs;
            }
#else      // Forcing based on flux derivative. 
            Array<OneD,NekDouble> ftmp;

            // Forcing  using  Direct form <\phi,[[qf ]]>
            offset = 0;
            for(e = 0; e < nexp; ++e)
            {
                nbndry = (*m_exp)[e]->NumDGBndryCoeffs();
                e_ncoeffs = (*m_exp)[e]->GetNcoeffs();

                // Get BndSysForce Matrix
                LocalRegions::MatrixKey Bmatkey(StdRegions::eUnifiedDGHelmBndSysForce,  (*m_exp)[e]->DetExpansionType(),*((*m_exp)[e]), lambda, tau);
                
                DNekScalMat &BndSys = *((*m_exp)[e]->GetLocMatrix(Bmatkey));                               
                DNekVec vin (e_ncoeffs,ftmp = f+offset,eWrapper);
                DNekVec vout(nbndry);
                
                vout = BndSys*vin;

                // Subtract vout from forcing terms 
                cnt = 0;
                for(i = 0; i < nbndry; ++i)
                {
                    id = m_elmtTraceMap[e][i];
                    m_trace->SetCoeff(id, m_trace->GetCoeff(id)
                                      - m_elmtTraceSign[e][cnt]*vout[i]);
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
            for(e = 0; e < nexp; ++e)
            {
                nbndry = (*m_exp)[e]->NumDGBndryCoeffs();
                // check to see if element has Dirichlet boundary
                // Probably could use a quicker check here
                if(Vmath::Vmin(nbndry,m_elmtTraceMap[e],1) < NumDirichlet)
                {
                    // Get BndSys Matrix

                    LocalRegions::MatrixKey Bmatkey(StdRegions::eUnifiedDGHelmBndSys,  (*m_exp)[e]->DetExpansionType(),*((*m_exp)[e]), lambda, tau);
                    
                    DNekScalMat &BndSys = *((*m_exp)[e]->GetLocMatrix(Bmatkey));               
                    DNekVec vin (nbndry);
                    DNekVec vout(nbndry);


                    // Set up Edge Dirichlet Values
                    cnt = 0;
                    for(i = 0; i < (*m_exp)[e]->GetNedges(); ++i)
                    {
                        e_ncoeffs = (*m_exp)[e]->GetEdgeNcoeffs(i); 
                        if(m_elmtTraceMap[e][cnt] < NumDirichlet)
                        {
                            for(j = 0; j < e_ncoeffs; ++j)
                            {
                                vin[cnt] = m_elmtTraceSign[e][cnt]*m_trace->GetCoeffs(m_elmtTraceMap[e][cnt]); 
                                cnt++;
                            } 
                        }
                        else
                        {
                            for(j = 0; j < e_ncoeffs; ++j)
                            {
                                vin[cnt] = 0.0;
                                cnt++;
                            }
                        }
                    }
                    
                    vout = BndSys*vin;
                    
                    // Subtract vout from forcing terms 
                    cnt = 0;
                    for(i = 0; i < nbndry; ++i)
                    {
                        if((id = m_elmtTraceMap[e][cnt]) >= NumDirichlet)
                        {
                            m_trace->SetCoeff(id, m_trace->GetCoeff(id)
                                           - m_elmtTraceSign[e][cnt]*vout[i]);
                        }
                        cnt++;
                    }
                }
            }
            
            if(GloBndDofs - NumDirichlet > 0)
            {
                GlobalLinSysKey       key(StdRegions::eUnifiedDGHelmBndSys,
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
                
#if 1
                for(j = 0; j < nbndry; ++j)
                {
                    e_lambda[j] = m_elmtTraceSign[i][j]*m_trace->GetCoeffs(m_elmtTraceMap[i][j]);
                }
#else
                int nedge,  cnt2 = 0;
                m_trace->PutCoeffsInToElmtExp();
                
                for(cnt2 = j = 0; j < (*m_exp)[i]->GetNedges(); ++j, cnt2 += nedge)
                {
                    nedge  = m_elmtToTrace[i][j]->GetNcoeffs();
                    Vmath::Vcopy(nedge, &(m_elmtToTrace[i][j]->GetCoeffs())[0],1,
                                 &(e_lambda)[0]+cnt2,1);
                }
#endif          
                (*m_exp)[i]->AddUDGHelmholtzTraceTerms(tau, e_lambda, 
                                                       m_elmtToTrace[i], 
                                                       e_f);
                cnt  += (*m_exp)[i]->GetNcoeffs();
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

        void DisContField2D::GetInnerOuterTracePhys(Array<OneD,NekDouble> &Inner, 
                                                    Array<OneD,NekDouble> &Outer)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            StdRegions::EdgeOrientation edgedir;
            int nquad_e,cnt,n,e,npts,offset;
            Array<OneD,NekDouble> e_tmp;

            // zero vectors; 
            Vmath::Zero(Inner.num_elements(),Inner,1);
            Vmath::Zero(Outer.num_elements(),Outer,1);

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eForwards)
                    {
                        offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e,(*m_exp)[n]->GetPhys(),
                                                     e_tmp = Inner + offset);
                    }
                }
            }
            
            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eBackwards)
                    {
                        offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e,(*m_exp)[n]->GetPhys(),
                                                     e_tmp = Outer+offset);
                    }
                }
            }
            
            // Since the Dirichlet boundaries are order first in the
            // trace space add these two vectors together and then
            // fill in Boundary conditions into the Outer array.
            Vmath::Vadd(m_numTraceDirichletPhysBCs,Inner,1,Outer,1,Inner,1);

            for(cnt = n = 0; n < m_bndConstraint.num_elements(); ++n)
            {
                npts = m_bndConstraint[n]->GetNpoints();
                Vmath::Vcopy(npts, m_bndConstraint[n]->GetPhys(),1,
                             e_tmp = Outer + cnt,1);
                cnt += npts;
            }
            
        }
        
        void DisContField2D::ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int nquad_e,cnt,n,e,npts,offset;
            Array<OneD,NekDouble> e_tmp;

            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            ASSERTL1(outarray.num_elements() >= m_trace->GetNpoints(),
                     "input array is of insufficient length");

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    offset = m_trace->GetPhys_Offset(m_elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->GetEdgePhysVals(e,(*m_exp)[n]->GetPhys(),
                                                 e_tmp = outarray + offset);
                }
            }
        }
            

        /// Note this routine changes m_trace->m_coeffs space; 
        void DisContField2D::AddTraceIntegral(Array<OneD, const NekDouble> &inarray, 
                                              Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset;
            Array<OneD, NekDouble> e_outarray;

            m_trace->PutPhysInToElmtExp(inarray);
            
            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetPhys_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    (*m_exp)[n]->AddEdgeBoundaryInt(e,m_elmtToTrace[n][e],
                                                    e_outarray = outarray+offset);
                }    
            }
        }

    } // end of namespace
} //end of namespace
