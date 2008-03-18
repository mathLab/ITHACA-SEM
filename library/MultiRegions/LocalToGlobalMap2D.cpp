///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap2D.cpp
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
// Description: Local to Global mapping routines in 2D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalMap2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        LocalToGlobalMap2D::LocalToGlobalMap2D()
        {
        }

        LocalToGlobalMap2D::LocalToGlobalMap2D(const int loclen, 
                                               const StdRegions::StdExpansionVector &locexp, 
                                               const SpatialDomains::MeshGraph2D &graph2D)
        {
            int i,j,k;
            int gid = 0;
            int vid;
            int eid;
            int cnt = 0;
            int vcnt = 0;
            int ecnt = 0;
            int nedge, nedge_coeffs;
            StdRegions::EdgeOrientation eorient;
            LibUtilities::BasisType Btype;

            // set up Local to Continuous mapping 
            StdRegions::StdExpMap vmap;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr locTriExp;
           
            m_totLocDofs = loclen; 
            m_locToContMap =  Array<OneD, int>(m_totLocDofs,-1);

            // Reserve storage for the re-ordering (give a new vertex and edge id) of the vertices and edges.
            // This is needed because the set-up of the mapping is based on the vertex and edge id's, but
            // because the domain does not necessarily encompassed the entire meshgraph, the original id's
            // might be unsuitable.
            Array<OneD, int> renumbVerts(graph2D.GetNvertices(),-1);
            Array<OneD, int> renumbEdges(graph2D.GetNseggeoms(),-1);

            m_totLocBndDofs = 0;
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < locQuadExp->GetNverts(); ++j)
                    {   
                        vid = (locQuadExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = vcnt++;
                        }
                        
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }  
                    }
                    m_totLocBndDofs += locQuadExp->NumBndryCoeffs();
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < locTriExp->GetNverts(); ++j)
                    {    
                        vid = (locTriExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = vcnt++;
                        }
                        
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }  
                    }
                    m_totLocBndDofs += locTriExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
            }
            m_locToContBndMap = Array<OneD, int>(m_totLocBndDofs,-1);

            // Calculate the number of DOFs for every edge           
            Array<OneD, int> edge_offset(ecnt+1);
            m_sign_change = false;

            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < locQuadExp->GetNedges(); ++j)
                    {   
                        //set up nedge coefficients for each edge     
                        nedge_coeffs = locQuadExp->GetEdgeNcoeffs(j);
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs-2;
                        
                        // need a sign vector if edge_nceoff >=4 
                        if((nedge_coeffs >= 4)&&(locQuadExp->GetEdgeBasisType(0) ==  LibUtilities::eModified_A))
                        {
                            m_sign_change = true;
                        }
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < locTriExp->GetNedges(); ++j)
                    {                    
                        //set up nedge coefficients for each edge     
                        nedge_coeffs = locTriExp->GetEdgeNcoeffs(j);
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs-2;
                        
                        // need a sign vector if edge_nceoff >=4 
                        if((nedge_coeffs >= 4)&&(locTriExp->GetEdgeBasisType(0) ==  LibUtilities::eModified_A))
                        {
                            m_sign_change = true;
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
            }
            
            // set up sign vector
            if(m_sign_change)
            {
                m_sign = Array<OneD, NekDouble>(m_totLocDofs,1.0);
                m_bndSign = Array<OneD, NekDouble>(m_totLocBndDofs,1.0);
            }
            
            edge_offset[0] = vcnt; 
            
            // set up consecutive list for edge entries
            for(i = 1; i < ecnt; ++i)
            {
                edge_offset[i] += edge_offset[i-1];
            }            
            
            // set up simple map;
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < (nedge=locQuadExp->GetNedges()); ++j)
                    {
                        locQuadExp->MapTo_ModalFormat(nedge_coeffs = locQuadExp->GetEdgeNcoeffs(j),
                                                      Btype = locQuadExp->GetEdgeBasisType(j), j,
                                                      eorient = (locQuadExp->GetGeom())->GetEorient(j),
                                                      vmap);
                        
                        // set edge ids before setting vertices 
                        // so that can be reverse for nodal expansion when j>2
                        
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        
                        for(k = 2; k < nedge_coeffs; ++k)
                        {
                            m_locToContMap[cnt+vmap[k]] =  edge_offset[renumbEdges[eid]]+(k-2);
                        }
                        
                        // vmap is set up according to cartesian
                        // coordinates so have to change setting
                        // depending on which edge we are considering
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    renumbVerts[(locQuadExp->GetGeom())->GetVid(j)];
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    renumbVerts[(locQuadExp->GetGeom())->GetVid(j)];
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }   
                            }
                        }
                        else
                        {
                            // set edge global ids
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    renumbVerts[(locQuadExp->GetGeom())->GetVid(j)];
                                
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    renumbVerts[(locQuadExp->GetGeom())->GetVid(j)];
                                
                                if(Btype == LibUtilities::eGLL_Lagrange)
                                {
                                    for(k = 2; k < nedge_coeffs; ++k)
                                    {
                                        m_locToContMap[cnt+vmap[nedge_coeffs+1-k]] = edge_offset[renumbEdges[eid]]+(k-2);
                                    }
                                }
                            }
                        }
                        
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < (nedge=locTriExp->GetNedges()); ++j)
                    {
                        locTriExp->MapTo_ModalFormat(nedge_coeffs = locTriExp->GetEdgeNcoeffs(j),
                                                     Btype = locTriExp->GetEdgeBasisType(j), j,
                                                     eorient = (locTriExp->GetGeom())->GetEorient(j),
                                                     vmap);
                        
                        // set edge ids before setting vertices 
                        // so that can be reverse for nodal expansion when j>2
                        
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        
                        for(k = 2; k < nedge_coeffs; ++k)
                        {
                            m_locToContMap[cnt+vmap[k]] =  edge_offset[renumbEdges[eid]]+(k-2);
                        }
                        
                        // vmap is set up according to cartesian
                        // coordinates so have to change setting
                        // depending on which edge we are considering
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    renumbVerts[(locTriExp->GetGeom())->GetVid(j)];
                            }
                        else
                        {
                            m_locToContMap[cnt+vmap[1]] = 
                                renumbVerts[(locTriExp->GetGeom())->GetVid(j)];
                            if(m_sign_change)
                            {
                                for(k = 3; k < nedge_coeffs; k+=2)
                                {
                                    m_sign[cnt+vmap[k]] = -1;
                                }
                            }   
                        }
                        }
                        else
                        {
                            // set edge global ids
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    renumbVerts[(locTriExp->GetGeom())->GetVid(j)];
                                
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    renumbVerts[(locTriExp->GetGeom())->GetVid(j)];
                                
                                if(Btype == LibUtilities::eGLL_Lagrange)
                                {
                                    for(k = 2; k < nedge_coeffs; ++k)
                                    {
                                        m_locToContMap[cnt+vmap[nedge_coeffs+1-k]] = edge_offset[renumbEdges[eid]]+(k-2);
                                    }
                                }
                            }
                        }
                        
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
                cnt += locexp[i]->GetNcoeffs();            
            }
                
            gid = Vmath::Vmax(m_totLocDofs,&m_locToContMap[0],1)+1;
            m_totGloBndDofs = gid;
            
            cnt = 0;
            // setup interior mapping 
            for(i = 0; i < m_totLocDofs; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = gid++;
                }
                else
                {
                    if(m_sign_change)
                    {
                        m_bndSign[cnt]=m_sign[i];
                    }
                    m_locToContBndMap[cnt++]=m_locToContMap[i];
                }
            }
            m_totGloDofs = gid;
        }    
        
        LocalToGlobalMap2D::LocalToGlobalMap2D(const int loclen, 
                                               const StdRegions::StdExpansionVector &locexp, 
                                               const SpatialDomains::MeshGraph2D &graph2D,
                                               const ConstArray<OneD,MultiRegions::ExpList1DSharedPtr> &bndCondExp,
                                               const ConstArray<OneD,SpatialDomains::BoundaryConditionType> &bndCondTypes)
        {
            int i,j,k;
            int vid;
            int vid2;
            int eid;
            int j2;
            int bndID;
            int gid = 0;
            int cnt = 0;
            int cnt2;
            int vcnt = 0;
            int ecnt = 0;
            int nedge, nedge_coeffs;
            int nDirVerts = 0;
            int nDirEdges = 0;
            int nonDirOffset;
            bool bndCondEdge = false;
            int nLocBndCondDofs = 0;
            int nBndCondEdges = 0;
            StdRegions::EdgeOrientation eorient;
            LibUtilities::BasisType Btype;

            StdRegions::StdExpMap vmap;
            StdRegions::StdExpMap vmapEdge;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr locTriExp;
            LocalRegions::NodalTriExpSharedPtr locNodalTriExp;
            LocalRegions::SegExpSharedPtr bndSegExp;

            MultiRegions::ExpList1DSharedPtr bndExpList;
           
            m_totLocDofs = loclen; 
            m_locToContMap =  Array<OneD, int>(m_totLocDofs,-1);

            // Reserve storage for the re-ordering (give a new vertex and edge id) of the vertices and edges.
            // This is needed because the set-up of the mapping is based on the vertex and edge id's, but
            // because the domain does not necessarily encompassed the entire meshgraph, the original id's
            // might be unsuitable.
            // We will also use these array's to number the dirichlet vertices and edges first.
            Array<OneD, int> renumbVerts(graph2D.GetNvertices(),-1);
            Array<OneD, int> renumbEdges(graph2D.GetNseggeoms(),-1);

            // This array is used to indicate whether an edge is part of the boundary of the domain.
            Array<OneD, int> bndCondEdgeID(graph2D.GetNseggeoms(),-1);
            
            // Calculate the number of boundary condition edges
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                nBndCondEdges += (bndCondExp[i])->GetExpSize();
            }

            // This array should not be needed once the routine to retrieve the element en local edge-id of a
            // certain edge
            Array<OneD,LocalRegions::SegExpSharedPtr> bndCondEdges(nBndCondEdges);

            // Order the Dirichlet vertices and edges first.
            Array<OneD, int> bndCond_offset(nBndCondEdges);
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                bndExpList = bndCondExp[i];
                for(j = 0; j < bndExpList->GetExpSize(); j++)
                {
                    if(!(bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndExpList->GetExp(j))))
                    {
                        ASSERTL0(false,"dynamic cast to a SegExp failed");
                    }
                    eid = (bndSegExp->GetGeom())->GetEid();
                    bndCondEdgeID[eid] = cnt;
                    bndCond_offset[cnt] = nLocBndCondDofs;
                    bndCondEdges[cnt] = bndSegExp;
                    nLocBndCondDofs += bndSegExp->GetNcoeffs();
                    
                    if(bndCondTypes[i]==SpatialDomains::eDirichlet)
                    {
                        renumbEdges[eid] = ecnt++;  
                        nDirEdges++;
                        for(k = 0; k < bndSegExp->GetNverts(); k++)
                        {
                            vid = (bndSegExp->GetGeom())->GetVid(k);
                            if(renumbVerts[vid]==-1)
                            {
                                renumbVerts[vid] = vcnt++;
                                nDirVerts++;
                            }
                        }
                    }
                    cnt++;
                }
            }
            m_locToContBndCondMap = Array<OneD,int>(nLocBndCondDofs,-1);

            // Now order all other vertices and edges in renumbVerts and renumbEdges
            m_totLocBndDofs = 0;
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < locQuadExp->GetNverts(); ++j) // number verts = number edges for 2D geom
                    {   
                        vid = (locQuadExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = vcnt++;
                        }
                        
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }                       
                    }
                    m_totLocBndDofs += locQuadExp->NumBndryCoeffs();
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < locTriExp->GetNverts(); ++j) // number verts = number edges for 2D geom
                    {                    
                        vid = (locTriExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = vcnt++;
                        }
                        
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }   
                    }
                    m_totLocBndDofs += locTriExp->NumBndryCoeffs();
                }
                else if(locNodalTriExp = boost::dynamic_pointer_cast<LocalRegions::NodalTriExp>(locexp[i]))
                {
                    for(j = 0; j < locNodalTriExp->GetNverts(); ++j) // number verts = number edges for 2D geom
                    {                    
                        vid = (locNodalTriExp->GetGeom())->GetVid(j);
                        if(renumbVerts[vid]==-1)
                        {
                            renumbVerts[vid] = vcnt++;
                        }
                        
                        eid = (locNodalTriExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }   
                    }
                    m_totLocBndDofs += locNodalTriExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
            }
            m_locToContBndMap = Array<OneD, int>(m_totLocBndDofs,-1);
            
            // Calculate the number of DOFs for every edge
            Array<OneD, int> edge_offset(ecnt+1);
            m_sign_change = false;

            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < locQuadExp->GetNedges(); ++j)
                    {   
                        //set up nedge coefficients for each edge     
                        nedge_coeffs = locQuadExp->GetEdgeNcoeffs(j);
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs-2;
                        
                        // need a sign vector if edge_nceoff >=4 
                        if((nedge_coeffs >= 4)&&(locQuadExp->GetEdgeBasisType(0) == 
                                                 LibUtilities::eModified_A))
                        {
                            m_sign_change = true;
                        }
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < locTriExp->GetNedges(); ++j)
                    {                    
                        //set up nedge coefficients for each edge     
                        nedge_coeffs = locTriExp->GetEdgeNcoeffs(j);
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs-2;
                        
                        // need a sign vector if edge_nceoff >=4 
                        if((nedge_coeffs >= 4)&&(locTriExp->GetEdgeBasisType(0) == 
                                                 LibUtilities::eModified_A))
                        {
                            m_sign_change = true;
                        }
                    }
                }
                else if(locNodalTriExp = boost::dynamic_pointer_cast<LocalRegions::NodalTriExp>(locexp[i]))
                {
                    for(j = 0; j < locNodalTriExp->GetNedges(); ++j)
                    {                    
                        //set up nedge coefficients for each edge     
                        nedge_coeffs = locNodalTriExp->GetEdgeNcoeffs(j);
                        eid = (locNodalTriExp->GetGeom())->GetEid(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs-2;
                        
//                         // need a sign vector if edge_nceoff >=4 
//                         if((nedge_coeffs >= 4)&&(locNodalTriExp->GetEdgeBasisType(0) == 
//                                                  LibUtilities::eModified_A))
//                         {
//                             m_sign_change = true;
//                         }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
            }
            
            // set up sign vector
            if(m_sign_change)
            {
                m_sign = Array<OneD, NekDouble>(m_totLocDofs,1.0);
                m_bndSign = Array<OneD, NekDouble>(m_totLocBndDofs,1.0);
            }

            edge_offset[0] = nDirVerts;

            if(ecnt>nDirEdges)
            {
                edge_offset[nDirEdges] += (vcnt-nDirVerts);
            }
            
            // set up consecutive list for edge entries
            for(i = 1; i < ecnt+1; ++i)
            {
                edge_offset[i] += edge_offset[i-1];
            }

            m_numDirichletDofs = edge_offset[nDirEdges]-(vcnt-nDirVerts);  
            nonDirOffset = m_numDirichletDofs - nDirVerts;
            
            cnt = 0;
            // Set up the mapping
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < (nedge=locQuadExp->GetNedges()); ++j)
                    {
                        locQuadExp->MapTo_ModalFormat(nedge_coeffs = locQuadExp->GetEdgeNcoeffs(j),
                                                      Btype = locQuadExp->GetEdgeBasisType(j), j,
                                                      eorient = (locQuadExp->GetGeom())->GetEorient(j),
                                                      vmap);
                        
                        // set edge ids before setting vertices 
                        // so that can be reverse for nodal expansion when j>2                        
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        vid = renumbVerts[(locQuadExp->GetGeom())->GetVid(j)];

                        bndCondEdge = false;
                        if((bndID = bndCondEdgeID[eid])!=-1)
                        {
                            bndCondEdge = true;
                            cnt2 = bndCond_offset[bndID];
                            bndCondEdges[bndID]->MapTo(eorient,vmapEdge);

                            m_locToContBndCondMap[cnt2+vmapEdge[0]] = 
                                (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                            j2 = (j==3) ? 0 : (j+1);                     
                            vid2 = renumbVerts[(locQuadExp->GetGeom())->GetVid(j2)];
                            m_locToContBndCondMap[cnt2+vmapEdge[1]] = 
                                (vid2 < nDirVerts) ? (vid2) : (nonDirOffset + vid2);

                            int bndedgecnt = 0;
                            for(k = 0; k < nedge_coeffs; ++k)
                            {
                                if(m_locToContBndCondMap[cnt2+k] == -1)
                                {
                                    m_locToContBndCondMap[cnt2+k] = edge_offset[renumbEdges[eid]]+bndedgecnt;  
                                    bndedgecnt++;
                                }
                            }                            
                        }
                        for(k = 2; k < nedge_coeffs; ++k)
                        {
                            m_locToContMap[cnt+vmap[k]] =  edge_offset[renumbEdges[eid]]+(k-2);
                        }  
                                
                                
                        // vmap is set up according to cartesian
                        // coordinates so have to change setting
                        // depending on which edge we are considering
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }   
                            }
                        }
                        else
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);
                            }
                            if(Btype == LibUtilities::eGLL_Lagrange)
                            {
                                for(k = 2; k < nedge_coeffs; ++k)
                                {
                                    m_locToContMap[cnt+vmap[nedge_coeffs+1-k]] = edge_offset[renumbEdges[eid]]+(k-2);
                                }
                            }
                        }                      
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < (nedge=locTriExp->GetNedges()); ++j)
                    {
                        locTriExp->MapTo_ModalFormat(nedge_coeffs = locTriExp->GetEdgeNcoeffs(j),
                                                     Btype = locTriExp->GetEdgeBasisType(j), j,
                                                     eorient = (locTriExp->GetGeom())->GetEorient(j),
                                                     vmap);
                        
                        // set edge ids before setting vertices 
                        // so that can be reverse for nodal expansion when j>2
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        vid = renumbVerts[(locTriExp->GetGeom())->GetVid(j)];

                        bndCondEdge = false;
                        if((bndID = bndCondEdgeID[eid])!=-1)
                        {
                            bndCondEdge = true;
                            cnt2 = bndCond_offset[bndID];
                            bndCondEdges[bndID]->MapTo(eorient,vmapEdge);

                            m_locToContBndCondMap[cnt2+vmapEdge[0]] = 
                                (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                            j2 = (j==2) ? 0 : (j+1);                     
                            vid2 = renumbVerts[(locTriExp->GetGeom())->GetVid(j2)];
                            m_locToContBndCondMap[cnt2+vmapEdge[1]] = 
                                (vid2 < nDirVerts) ? (vid2) : (nonDirOffset + vid2);
                        }
                        
                        for(k = 2; k < nedge_coeffs; ++k)
                        {
                            m_locToContMap[cnt+vmap[k]] =  edge_offset[renumbEdges[eid]]+(k-2);
                            if(bndCondEdge)
                            {
                                m_locToContBndCondMap[cnt2+k] = edge_offset[renumbEdges[eid]]+(k-2);
                            }
                        }
                        
                        // vmap is set up according to cartesian
                        // coordinates so have to change setting
                        // depending on which edge we are considering
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);
                            }
                            else
                            {

                                m_locToContMap[cnt+vmap[1]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }   
                            }
                        }
                        else
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);
                                
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);                          
                            }
                        }
                        
                    }
                }
                else if(locNodalTriExp = boost::dynamic_pointer_cast<LocalRegions::NodalTriExp>(locexp[i]))
                {
                    for(j = 0; j < (nedge=locNodalTriExp->GetNedges()); ++j)
                    {
                        locNodalTriExp->MapTo_ModalFormat(nedge_coeffs = locNodalTriExp->GetEdgeNcoeffs(j),
                                                     Btype = locNodalTriExp->GetEdgeBasisType(j), j,
                                                     eorient = (locNodalTriExp->GetGeom())->GetEorient(j),
                                                     vmap);
                        
                        // set edge ids before setting vertices 
                        // so that can be reverse for nodal expansion when j>2
                        eid = (locNodalTriExp->GetGeom())->GetEid(j);
                        vid = renumbVerts[(locNodalTriExp->GetGeom())->GetVid(j)];

                        bndCondEdge = false;
                        if((bndID = bndCondEdgeID[eid])!=-1)
                        {
                            bndCondEdge = true;
                            cnt2 = bndCond_offset[bndID];
                            bndCondEdges[bndID]->MapTo(eorient,vmapEdge);

                            m_locToContBndCondMap[cnt2+vmapEdge[0]] = 
                                (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                            j2 = (j==2) ? 0 : (j+1);                     
                            vid2 = renumbVerts[(locNodalTriExp->GetGeom())->GetVid(j2)];
                            m_locToContBndCondMap[cnt2+vmapEdge[1]] = 
                                (vid2 < nDirVerts) ? (vid2) : (nonDirOffset + vid2);
                                
                            int bndedgecnt = 0;
                            for(k = 0; k < nedge_coeffs; ++k)
                            {
                                if(m_locToContBndCondMap[cnt2+k] == -1)
                                {
                                    m_locToContBndCondMap[cnt2+k] = edge_offset[renumbEdges[eid]]+bndedgecnt;  
                                    bndedgecnt++;
                                }
                            }  
                        }
                        
                        for(k = 2; k < nedge_coeffs; ++k)
                        {
                            m_locToContMap[cnt+vmap[k]] =  edge_offset[renumbEdges[eid]]+(k-2);
                        }
                        
                        // vmap is set up according to cartesian
                        // coordinates so have to change setting
                        // depending on which edge we are considering
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);
                            }
                            else
                            {

                                m_locToContMap[cnt+vmap[1]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);

                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }   
                            }
                        }
                        else
                        {
                            if(eorient == StdRegions::eForwards)
                            {
                                m_locToContMap[cnt+vmap[1]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);
                                
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_sign[cnt+vmap[k]] = -1;
                                    }
                                }
                            }
                            else
                            {
                                m_locToContMap[cnt+vmap[0]] = 
                                    (vid < nDirVerts) ? (vid) : (nonDirOffset + vid);                          
                            }
//                             if(Btype == LibUtilities::eGLL_Lagrange)
//                             {
                                 for(k = 2; k < nedge_coeffs; ++k)
                                 {
                                     m_locToContMap[cnt+vmap[nedge_coeffs+1-k]] = edge_offset[renumbEdges[eid]]+(k-2);
                                 }
                                //                           }
                        }
                        
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
                cnt += locexp[i]->GetNcoeffs();            
            }
                
            gid = Vmath::Vmax(m_totLocDofs,&m_locToContMap[0],1)+1;
            m_totGloBndDofs = gid;
            
            cnt=0;
            // setup interior mapping 
            for(i = 0; i < m_totLocDofs; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = gid++;
                }
                else
                {
                    if(m_sign_change)
                    {
                        m_bndSign[cnt]=m_sign[i];
                    }
                    m_locToContBndMap[cnt++]=m_locToContMap[i];
                }
            }
            m_totGloDofs = gid;           
        }    
        
        LocalToGlobalMap2D::~LocalToGlobalMap2D()
        {
        } 
    }
}

/**
* $Log: LocalToGlobalMap2D.cpp,v $
* Revision 1.9  2008/01/23 21:50:52  sherwin
* Update from EdgeComponents to SegGeoms
*
* Revision 1.8  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.7  2007/10/03 11:37:51  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.6  2007/07/20 02:04:13  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.5  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
