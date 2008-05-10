///////////////////////////////////////////////////////////////////////////////
//
// File LocaltoGlobalBndryMap2D.cpp
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

#include <MultiRegions/LocalToGlobalBndryMap2D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        LocalToGlobalBndryMap2D::LocalToGlobalBndryMap2D(const StdRegions::StdExpansionVector &locexp,  
                                                         const SpatialDomains::MeshGraph1D &graph2D,
                                                         const ConstArray<OneD,LocalRegions::ExpList1DSharedPtr> &bndCondExp,
                                                         const ConstArray<OneD,SpatialDomains::BoundaryConditionType> &bndCondTypes)
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
            
            // set up mapping from continuous edge (lambda) space to elements 
            StdRegions::StdExpMap vmap;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr  locTriExp;
              
            // Reserve storage for the re-ordering (give a new vertex
            // and edge id) of the vertices and edges.  This is needed
            // because the set-up of the mapping is based on the
            // vertex and edge id's, but because the domain does not
            // necessarily encompassed the entire meshgraph, the
            // original id's might be unsuitable.
            Array<OneD, int> renumbEdges(graph2D.GetNseggeoms(),-1);

            // Calculate the number of DOFs for every edge           
            Array<OneD, int> edge_offset(ecnt+1);
            m_sign_change = false;

            // Number Dirichlet edges first 
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                bndExpList = bndCondExp[i];
                for(j = 0; j < bndExpList->GetExpSize(); j++)
                {
                    if(!(bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndExpList->GetExp(j))))
                    {
                        ASSERTL0(false,"dynamic cast to a SegExp failed");
                    }
                    if(bndCondTypes[i] == SpatialDomains::eDirichlet)
                    {
                        eid = (bndSegExp->GetGeom())->GetEid();
                        renumbEdges[eid] = ecnt++;  
                    }
                }
            }
            
            m_totEdgeDofs = 0;
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < locQuadExp->GetNverts(); ++j)
                    {   
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }  

                        //set up nedge coefficients for each edge     
                        nedge_coeffs =  locQuadExp->GetEdgeNcoeffs(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs;
                        
                        // need a sign vector if edge_nceoff >=4 
                        if((nedge_coeffs >= 4)&&(locQuadExp->GetEdgeBasisType(0) ==  LibUtilities::eModified_A))
                        {
                            m_sign_change = true;
                        }

                        m_totEdgeDofs += nedge_coeffs;
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locexp[i]))
                {
                    for(j = 0; j < locTriExp->GetNverts(); ++j)
                    {    
                        eid = (locTriExp->GetGeom())->GetEid(j);
                        if(renumbEdges[eid]==-1)
                        {
                            renumbEdges[eid] = ecnt++;
                        }  

                        //set up nedge coefficients for each edge     
                        nedge_coeffs = locTriExp->GetEdgeNcoeffs(j);
                        edge_offset[renumbEdges[eid]+1] = nedge_coeffs;
                        
                        // need a sign vector if edge_nceoff >=4 
                        if((nedge_coeffs >= 4)&&(locTriExp->GetEdgeBasisType(0) ==  LibUtilities::eModified_A))
                        {
                            m_sign_change = true;
                        }
                        m_totLocBndDofs += nedge_coeffs; 
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
            }

            m_edgeToElmtMap = Array<OneD, int>(m_totEdgeDofs,-1);
            
            // set up sign vector
            if(m_sign_change)
            {
                m_bndSign = Array<OneD, NekDouble>(m_totEdgeDofs,1.0);
            }
            
            edge_offset[0] = 0; 
            
            // set up consecutive list for edge entries
            for(i = 1; i < ecnt; ++i)
            {
                edge_offset[i] += edge_offset[i-1];
            }            

           
            // set up simple map;
            vcnt = 0;
            for(i = 0; i < locexp.size(); ++i)
            {
                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locexp[i]))
                {
                    for(j = 0; j < (nedge=locQuadExp->GetNedges()); ++j)
                    {

                        // set up vertices continuously as we go through elements

                        locQuadExp->MapTo_ModalFormat(nedge_coeffs = locQuadExp->GetEdgeNcoeffs(j),
                                                      Btype = locQuadExp->GetEdgeBasisType(j), j,
                                                      eorient = (locQuadExp->GetGeom())->GetEorient(j),
                                                      vmap);
                        
                        
                        eid = (locQuadExp->GetGeom())->GetEid(j);
                        offset = edge_offset[renumEdge[eid]];
                        
                        for(k = 0; k < nedge_coeffs; ++k)
                        {
                            m_edgeToElmtMap[offset+k] = cnt+vmap[k]; 
                        }
                        

                        // Sort out sign or ordering issues.
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eBackwards)
                            {
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_bndSign[offset+k] = -1;
                                    }
                                }   
                            }
                        }
                        else
                        {
                            // set edge global ids
                            if(eorient == StdRegions::eForwards)
                            {
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    { 
                                        m_bndSign[offset + k] = -1;
                                    }
                                }
                            }
                            else
                            {
                                if(Btype == LibUtilities::eGLL_Lagrange)
                                {
                                    for(k = 0; k < nedge_coeffs; ++k)
                                    {
                                        m_edgeToElmtMap[offset + k] = cnt+vmap[nedge_coeffs-1-k];
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
                        offset = edge_offset[renumEdge[eid]];
                        
                        for(k = 0; k < nedge_coeffs; ++k)
                        {
                            m_edgeToElmtMap[offset+k] = cnt+vmap[k]; 
                        }
                        
                        // Sort out sign issues
                        if(j < 2)
                        {
                            if(eorient == StdRegions::eBackwards)
                            {
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_bndSign[cnt+vmap[k]] = -1;
                                    }
                                }   
                            }
                        }
                        else
                        {
                            // set edge global ids
                            if(eorient == StdRegions::eForwards)
                            {
                                if(m_sign_change)
                                {
                                    for(k = 3; k < nedge_coeffs; k+=2)
                                    {
                                        m_bndSign[cnt+vmap[k]] = -1;
                                    }
                                }
                            }
                            else
                            {
                                if(Btype == LibUtilities::eGLL_Lagrange)
                                {
                                    for(k = 0; k < nedge_coeffs; ++k)
                                    {
                                        m_locToContBndMap[cnt+vmap[nedge_coeffs-1-k]] = edge_offset[renumbEdges[eid]]+k;
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
                cnt += locexp[i]->NumBndryCoeffs(); 
            }
            
            gid = Vmath::Vmax(m_totLocBndDofs,&m_locToContBndMap[0],1)+1;
            m_totGloBndDofs = gid;
            
        }
 
        LocalToGlobalBndryMap2D::~LocalToGlobalBndryMap2D()
        {
        }
    }
}

/**
* $Log: LocalToGlobalBndryMap1D.cpp,v $
*
**/
