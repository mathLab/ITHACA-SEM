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
// Description: Field definition for 2D domain with boundary
// conditions using LDG flux
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {

        DisContField2D::DisContField2D(void):
            ExpList2D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField2D::DisContField2D(const DisContField2D &In):
            ExpList2D(In),
            m_bndCondExpansions   (In.m_bndCondExpansions),
            m_bndConditions       (In.m_bndConditions),
            m_globalBndMat        (In.m_globalBndMat),
            m_trace               (In.m_trace),
            m_traceMap            (In.m_traceMap)
        {
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       const GlobalSysSolnType solnType,
                                       bool SetUpJustDG):
            ExpList2D(graph2D),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            ApplyGeomInfo(graph2D);

            if(SetUpJustDG)
            {
                ASSERTL0(false, "Set up Trace space for no boundary conditions");
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>
                                                    ::AllocateSharedPtr();
//                map<int,int> periodicEdges;
//                vector<map<int,int> >periodicVertices;
//                GetPeriodicEdges(graph2D,bcs,bcs.GetVariable(bc_loc),
//                                 periodicVertices,periodicEdges);

                // Set up Trace space
/*                bool UseGenSegExp = true;
                m_trace = MemoryManager<ExpList1D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                *m_exp,graph2D, periodicEdges, UseGenSegExp);

                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(graph2D,m_trace,m_exp,solnType,
                                      m_bndCondExpansions,m_bndConditions,
                                      periodicEdges);
*/            }
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const int bc_loc,
                                       const GlobalSysSolnType solnType,
                                       bool SetUpJustDG):
            ExpList2D(graph2D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,
                                               bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph2D);

            if(SetUpJustDG)
            {
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>
                                                    ::AllocateSharedPtr();
                map<int,int> periodicEdges;
                vector<map<int,int> >periodicVertices;
                GetPeriodicEdges(graph2D,bcs,bcs.GetVariable(bc_loc),
                                 periodicVertices,periodicEdges);

                // Set up Trace space
                m_trace = MemoryManager<ExpList1D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                *m_exp,graph2D, periodicEdges);

                // Scatter trace segments to 2D elements. For each element,
                // we find the trace segment associated to each edge. The
                // element then retains a pointer to the trace space segments,
                // to ensure uniqueness of normals when retrieving from two
                // adjoining elements which do not lie in a plane.
                SpatialDomains::Geometry1DSharedPtr ElmtSegGeom;
                SpatialDomains::Geometry1DSharedPtr TraceSegGeom;
                for (int i = 0; i < m_exp->size(); ++i)
                {
                    for (int j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                    {
                        ElmtSegGeom  = ((*m_exp)[i]->GetGeom2D())->GetEdge(j);
                        for (int k = 0; k < m_trace->GetExpSize(); ++k)
                        {
                            TraceSegGeom = m_trace->GetExp(k)->GetGeom1D();
                            if (TraceSegGeom == ElmtSegGeom)
                            {
                                LocalRegions::Expansion2DSharedPtr exp2d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[i]);
                                StdRegions::StdExpansion1DSharedPtr exp1d
                                    = boost::dynamic_pointer_cast<StdRegions::StdExpansion1D>(m_trace->GetExp(k));

                                exp2d->SetEdgeExp(j,exp1d);
                                break;
                            }
                        }
                    }
                }

                // Finally set up the trace map between element edges and
                // trace segment expansions.
                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(graph2D,m_trace,*this,solnType,
                                      m_bndCondExpansions,m_bndConditions,
                                      periodicEdges);
                
            }
            else
            { 
                // set elmt edges to point to robin bc edges if required. 
                int i,cnt;
                Array<OneD, int> ElmtID,EdgeID;
                GetBoundaryToElmtMap(ElmtID,EdgeID);
                
                for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                {
                    MultiRegions::ExpList1DSharedPtr locExpList;
                    
                    if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
                    {
                        int e;                    
                        locExpList = m_bndCondExpansions[i];
                        
                        for(e = 0; e < locExpList->GetExpSize(); ++e)
                        {
                            LocalRegions::Expansion2DSharedPtr exp2d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[ElmtID[cnt+e]]);
                            StdRegions::StdExpansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast<StdRegions::StdExpansion1D>(locExpList->GetExp(e));
                            
                            exp2d->SetEdgeExp(EdgeID[cnt+e],exp1d);
                        }
                    }
                    cnt += m_bndCondExpansions[i]->GetExpSize();
                }
            }
        }

        DisContField2D::DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const std::string variable,
                                       const GlobalSysSolnType solnType,
                                       bool SetUpJustDG):

            ExpList2D(graph2D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph2D);

            if(SetUpJustDG)
            {
                // Set up matrix map
                m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

                map<int,int> periodicEdges;
                vector<map<int,int> > periodicVertices;
                GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);

                // Set up Trace space
                m_trace = MemoryManager<ExpList1D>::AllocateSharedPtr(m_bndCondExpansions,m_bndConditions,*m_exp,graph2D,periodicEdges);

                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(graph2D,m_trace,*this,solnType,
                                      m_bndCondExpansions,m_bndConditions, periodicEdges);
            }
        }


        void DisContField2D::GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                                                SpatialDomains::BoundaryConditions &bcs,
                                                                const std::string variable)
        {
            int i,cnt  = 0;
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            int cnt2 = 0;
            for(i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    cnt++;
                    if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() == SpatialDomains::eDirichlet )
                    {
                        cnt2++;
                    }
                }
            }
            bool UseGenSegExp = true;
            m_numDirBndCondExpansions = cnt2;
            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpList1DSharedPtr>(cnt);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);

            SetBoundaryConditionExpansion(graph2D,bcs,variable,m_bndCondExpansions,m_bndConditions);

            // Set up normals on non-Dirichlet boundary conditions
            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType()
                   != SpatialDomains::eDirichlet)
                {
                    m_bndCondExpansions[i]->SetUpPhysNormals(*m_exp);
                }
            }
        }

        DisContField2D::~DisContField2D()
        {
        }

        GlobalLinSysSharedPtr DisContField2D::GetGlobalBndLinSys(const GlobalLinSysKey &mkey)
        {
            ASSERTL0(mkey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");
            ASSERTL1(mkey.GetGlobalSysSolnType()!=eDirectFullMatrix,
                     "Full matrix global systems are not supported for HDG expansions");
            ASSERTL1(mkey.GetGlobalSysSolnType()==m_traceMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested solution type");

            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalBndMat->find(mkey);

            if(matrixIter == m_globalBndMat->end())
            {
                glo_matrix = GenGlobalBndLinSys(mkey,m_traceMap);
                (*m_globalBndMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

        // Construct the two trace vectors of the inner and outer
        // trace solution from the field contained in m_phys, where
        // the Weak dirichlet boundary conditions are listed in the
        // outer part of the vecotr
        void DisContField2D::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd,
                                                Array<OneD,NekDouble> &Bwd)
        {
            GetFwdBwdTracePhys(m_phys,Fwd,Bwd);
        }

        void DisContField2D::GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field,
                                                Array<OneD,NekDouble> &Fwd,
                                                Array<OneD,NekDouble> &Bwd)
        {
            // Loop over elements and collect forward expansion
            int nexp = GetExpSize();
            StdRegions::EdgeOrientation edgedir;
            int nquad_e,cnt,n,e,npts,offset, phys_offset;
            Array<OneD,NekDouble> e_tmp;

            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            // zero vectors;
            Vmath::Zero(Fwd.num_elements(),Fwd,1);
            Vmath::Zero(Bwd.num_elements(),Bwd,1);

            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eForwards)
                    {
                        offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Fwd + offset);
                    }
                }
            }

            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eBackwards)
                    {
                        offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Bwd + offset);
                    }
                }
            }

            // fill boundary conditions into missing elements
            int id1,id2 = 0;
            cnt = 0;
            for(n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                if(m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {

                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0);

                        if(m_traceMap->GetBndExpAdjacentOrient(cnt+e) == eAdjacentEdgeIsForwards)
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e) ;
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&(m_bndCondExpansions[n]->GetPhys())[id1],1,&Bwd[id2],1);
                        }
                        else
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e) ;
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&(m_bndCondExpansions[n]->GetPhys())[id1],1,&Fwd[id2],1);
                        }
                    }

                    cnt +=e;
                }
                else if((m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eNeumann)||(m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eRobin))
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0);

                        if(m_traceMap->GetBndExpAdjacentOrient(cnt+e) == eAdjacentEdgeIsForwards)
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e);
                            ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1] == 0.0,"method not set up for non-zero Neumann boundary condition");
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&Fwd[id2],1,&Bwd[id2],1);
                        }
                        else
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e);
                            ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1] == 0.0,"method not set up for non-zero Neumann boundary condition");
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&Bwd[id2],1,&Fwd[id2],1);
                        }
                    }

                    cnt +=e;
                }
                else
                {
                    ASSERTL0(false,"method not set up for non-Dirichlet conditions");
                }
            }

        }

        void DisContField2D::ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {

            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            ExtractTracePhys(m_phys, outarray);
        }

        void DisContField2D::ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int nquad_e,n,e,offset,phys_offset;
            Array<OneD,NekDouble> e_tmp;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            ASSERTL1(outarray.num_elements() >= m_trace->GetNpoints(),
                     "input array is of insufficient length");

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->GetEdgePhysVals(e,  elmtToTrace[n][e],
                                                 inarray + phys_offset,
                                                 e_tmp = outarray + offset);
                }
            }
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField2D::AddTraceIntegral(const Array<OneD, const NekDouble> &Fx,
                                              const Array<OneD, const NekDouble> &Fy,
                                              Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fx + t_offset,
                                                        Fy + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField2D::AddTraceIntegral(const Array<OneD, const NekDouble> &Fn, Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fn + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }

        // Set up a list of element ids and edge ids that link to the
        // boundary conditions
        void DisContField2D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, Array<OneD,int> &EdgeID)
        {
            map<int, int> EdgeGID;
            int i,n,id;
            int bid,cnt,Eid;
            int nbcs = 0;

            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {
                nbcs += m_bndCondExpansions[i]->GetExpSize();
            }

            // make sure arrays are of sufficient length
            if(ElmtID.num_elements() != nbcs)
            {
                ElmtID = Array<OneD, int>(nbcs,-1);
            }
            else
            {
                fill(ElmtID.get(), ElmtID.get()+nbcs, -1);
            }

            if(EdgeID.num_elements() != nbcs)
            {
                EdgeID = Array<OneD, int>(nbcs);
            }

            // setup map of all global ids along boundary
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                for(i = 0; i < m_bndCondExpansions[n]->GetExpSize(); ++i)
                {
                    Eid =  m_bndCondExpansions[n]->GetExp(i)->GetGeom1D()->GetEid();
                    EdgeGID[Eid] = cnt++;
                }
            }


            // Loop over elements and find edges that match;
            for(cnt = n = 0; n < GetExpSize(); ++n)
            {
                for(i = 0; i < (*m_exp)[n]->GetNedges(); ++i)
                {
                    id = (*m_exp)[n]->GetGeom2D()->GetEid(i);

                    if(EdgeGID.count(id) > 0)
                    {
                        bid = EdgeGID.find(id)->second;
                        ASSERTL1(ElmtID[bid] == -1,"Edge already set");
                        ElmtID[bid] = n;
                        EdgeID[bid] = i;
                        cnt ++;
                    }
                }
            }

            ASSERTL1(cnt == nbcs,"Failed to visit all boundary condtiions");
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField2D::AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd,
                                                const Array<OneD, const NekDouble> &Bwd,
                                                Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryBiInt(e,elmtToTrace[n][e],
                                                        Fwd + t_offset,
                                                        Bwd + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }

        /** Calculate the L2 error of the Q_dir derivative using the
            consistent DG evaluation of Q_dir. The soln provided is of
            the primative variation at the quadrature points and the
            derivative is compared to the discrete derivative at these
            points which is likely to be undesireable unless using a
            much higher number of quadrature points than the
            polynomial order used to evaluate Q_dir
         */

       NekDouble DisContField2D::L2_DGDeriv(const int dir,
                                           const Array<OneD, const NekDouble> &soln)
      {

          int    i,e,ncoeff_edge;
          Array<OneD, const NekDouble> tmp_coeffs;
          Array<OneD, NekDouble> out_d(m_ncoeffs), out_tmp;

          Array<OneD, Array< OneD, StdRegions::StdExpansion1DSharedPtr> > elmtToTrace = m_traceMap->GetElmtToTrace();

          StdRegions::EdgeOrientation edgedir;

          int     eid,cnt;
          int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
          Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;
          m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

          edge_lambda = loc_lambda;
          // Calculate Q using standard DG formulation.
          for(i =cnt = 0; i < GetExpSize(); ++i)
          {
              eid = m_offset_elmt_id[i];
              // Probably a better way of setting up lambda than this.
              // Note cannot use PutCoeffsInToElmts since lambda space
              // is mapped during the solve.
              for(e = 0; e < (*m_exp)[eid]->GetNedges(); ++e)
              {
                  edgedir = (*m_exp)[eid]->GetEorient(e);

                  ncoeff_edge = elmtToTrace[eid][e]->GetNcoeffs();
                  elmtToTrace[eid][e]->SetCoeffsToOrientation(edgedir,edge_lambda,edge_lambda);
                  Vmath::Vcopy(ncoeff_edge,edge_lambda,1,elmtToTrace[eid][e]->UpdateCoeffs(),1);
                  edge_lambda = edge_lambda + ncoeff_edge;
              }

              (*m_exp)[eid]->DGDeriv(dir,tmp_coeffs = m_coeffs+cnt,
                                     elmtToTrace[eid],
                                     out_tmp = out_d+cnt);
              cnt  += (*m_exp)[eid]->GetNcoeffs();
          }
          BwdTrans(out_d,m_phys);
          Vmath::Vsub(m_npoints,m_phys,1,soln,1,m_phys,1);

          return L2();
      }

        void DisContField2D::v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveDG(inarray, outarray, lambda, varLambda, varCoeff, 1);
        }

        void DisContField2D::v_HelmSolveDG(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                          NekDouble tau)
        {
            int i,j,n,cnt,cnt1,nbndry;
            int nexp = GetExpSize();
            StdRegions::StdExpansionSharedPtr BndExp;

            Array<OneD,NekDouble> f(m_ncoeffs);
            DNekVec F(m_ncoeffs,f,eWrapper);
            Array<OneD,NekDouble> e_f, e_l;

            //----------------------------------
            //  Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(inarray,f);
            Vmath::Neg(m_ncoeffs,f,1);

            //----------------------------------
            //  Solve continuous flux System
            //----------------------------------
            int GloBndDofs   = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirichlet = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs,id;

            // linked data
            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,lambda,tau,varCoeff);
            const DNekScalBlkMatSharedPtr &HDGLamToU = GetBlockMatrix(HDGLamToUKey);

            Array<OneD,NekDouble> BndSol = m_trace->UpdateCoeffs();
            Array<OneD,NekDouble> BndRhs(GloBndDofs,0.0);
            // Zero trace space
            Vmath::Zero(GloBndDofs,BndSol,1);

            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);

            //----------------------------------
            // Evaluate Trace Forcing
            //----------------------------------

            // Determing <u_lam,f> terms using HDGLamToU matrix
            for(cnt = cnt1 = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[m_offset_elmt_id[n]]->NumDGBndryCoeffs();

                e_ncoeffs = (*m_exp)[m_offset_elmt_id[n]]->GetNcoeffs();
                e_f       = f + cnt;
                e_l       = loc_lambda + cnt1;

                // use outarray as tmp space
                DNekVec     Floc    (nbndry, e_l, eWrapper);
                DNekVec     ElmtFce (e_ncoeffs, e_f, eWrapper);
                Floc = Transpose(*(HDGLamToU->GetBlock(n,n)))*ElmtFce;

                cnt   += e_ncoeffs;
                cnt1  += nbndry;
            }

            // Assemble into global operator
            m_traceMap->AssembleBnd(loc_lambda,BndRhs);

            cnt = 0;
            // Copy Dirichlet boundary conditions into trace space
            for(i = 0; i < m_numDirBndCondExpansions; ++i)
            {
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                {
                    id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(cnt++);
                    BndSol[id] = m_bndCondExpansions[i]->GetCoeffs()[j];
                }
            }

            //Add weak boundary condition to trace forcing
            for(i = m_numDirBndCondExpansions; i < m_bndCondExpansions.num_elements(); ++i)
            {
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                {
                    id   = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(cnt++);
                    BndRhs[id] += m_bndCondExpansions[i]->GetCoeffs()[j];
                }
            }

            //----------------------------------
            // Solve trace problem
            //----------------------------------
            if(GloBndDofs - NumDirichlet > 0)
            {
                GlobalLinSysKey       key(StdRegions::eHybridDGHelmBndLam,
                                          m_traceMap,lambda,tau,
                                          varCoeff);                           
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                LinSys->Solve(BndRhs,BndSol,m_traceMap);
            }

            //----------------------------------
            // Internal element solves
            //----------------------------------
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,lambda,tau,varCoeff);
            const DNekScalBlkMatSharedPtr& InvHDGHelm = GetBlockMatrix(invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

            // get local trace solution from BndSol
            m_traceMap->GlobalToLocalBnd(BndSol,loc_lambda);

            //  out =  u_f + u_lam = (*InvHDGHelm)*f + (LamtoU)*Lam
            out = (*InvHDGHelm)*F + (*HDGLamToU)*LocLambda;
        }
        
        /** 
         * Search through the edge expansions and identify which ones
         * have Robin/Mixed type boundary conditions. If find a Robin
         * boundary then store the edge id of the boundary condition
         * and the array of points of the physical space boundary
         * condition which are hold the boundary condition primitive
         * variable coefficient at the quatrature points
         *
         * \return std map containing the robin boundary condition
         * info using a key of the element id 
         *
         * There is a next member to allow for more than one Robin
         * boundary condition per element
         */

        map<int, RobinBCInfoSharedPtr> DisContField2D::GetRobinBCInfo(void)
        {
            int i,cnt;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);
            
            for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                MultiRegions::ExpList1DSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
                {
                    int e,elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];
                    
                    for(e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        RobinBCInfoSharedPtr rInfo = MemoryManager<RobinBCInfo>::AllocateSharedPtr(EdgeID[cnt+e],Array_tmp = locExpList->GetPhys() + locExpList->GetPhys_Offset(e));
                        elmtid = ElmtID[cnt+e];
                        // make link list if necessary
                        if(returnval.count(elmtid) != 0) 
                        {
                            rInfo->next = returnval.find(elmtid)->second;
                        }
                        returnval[elmtid] = rInfo;
                    }
                }
                cnt += m_bndCondExpansions[i]->GetExpSize();
            }
            
            return returnval; 
        }
    } // end of namespace
} //end of namespace
