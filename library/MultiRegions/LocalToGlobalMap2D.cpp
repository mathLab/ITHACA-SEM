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
        LocalToGlobalMap2D::LocalToGlobalMap2D(const int loclen, 
                                 const StdRegions::StdExpansionVector &locexp, 
                                 const SpatialDomains::MeshGraph2D &graph2D)
        {
            int i,j,k,gid,cnt, nedge, nedge_coeffs,edgid;
            int nGloVerts = 0;
            StdRegions::EdgeOrientation eorient;
            LibUtilities::BasisType Btype;
           
            m_totLocLen = loclen; 
            m_locToContMap =  Array<OneD, int>(m_totLocLen,-1);
            
            // set up mapping based 
            StdRegions::StdExpMap vmap;

            // need to put into LocalToGloabl initialiser. 
            Array<OneD, int> edge_offset(graph2D.GetNecomps()+1);
            m_sign_change = false;

            // determine global vertex ids
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < locexp[i]->GetNverts(); ++j)
                {                    
                    nGloVerts = max(nGloVerts,graph2D.GetVidFromElmt(locexp[i]->DetShapeType(),j,i));
                    //set up nedge coefficients for each edge     
                    nedge_coeffs = locexp[i]->GetEdgeNcoeffs(j);
                    edge_offset[graph2D.GetEidFromElmt(locexp[i]->DetShapeType(),j,i)+1] = nedge_coeffs-2;
                        
                    // need a sign vector if edge_nceoff >=4 
                    if((nedge_coeffs >= 4)&&(locexp[i]->GetEdgeBasisType(0) == 
                         LibUtilities::eModified_A))
                    {
                        m_sign_change = true;
                    }
                }
            }
            nGloVerts++; // set actual value to one plus maximum id
            
            // set up sign vector
            if(m_sign_change)
            {
                m_sign = Array<OneD, NekDouble>(m_totLocLen,1.0);
            }

            edge_offset[0] = nGloVerts; 

            // set up consecutive list for edge entries
            for(i = 1; i < graph2D.GetNecomps(); ++i)
            {
                edge_offset[i] += edge_offset[i-1];
            }
            

            // set up simple map;
            cnt = 0;
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < (nedge=locexp[i]->GetNedges()); ++j)
                {
                    locexp[i]->MapTo_ModalFormat(nedge_coeffs = locexp[i]->GetEdgeNcoeffs(j),
                                                 Btype = locexp[i]->GetEdgeBasisType(j), j,
                                                 eorient = graph2D.GetEorientFromElmt(locexp[i]->DetShapeType(),j,i),
                                                 vmap);
                    
                    // set edge ids before setting vertices 
                    // so that can be reverse for nodal expansion when j>2
                    
                    edgid = graph2D.GetEidFromElmt(locexp[i]->DetShapeType(),j,i);
                    
                    for(k = 2; k < nedge_coeffs; ++k)
                    {
                        m_locToContMap[cnt+vmap[k]] =  edge_offset[edgid]+(k-2);
                    }
                    
                    // vmap is set up according to cartesian
                    // coordinates so have to change setting
                    // depending on which edge we are considering
                    if(j < 2)
                    {
                        if(eorient == StdRegions::eForwards)
                        {
                            m_locToContMap[cnt+vmap[0]] = 
                                graph2D.GetVidFromElmt(locexp[i]->DetShapeType(),j,i);
                        }
                        else
                        {
                            m_locToContMap[cnt+vmap[1]] = 
                                graph2D.GetVidFromElmt(locexp[i]->DetShapeType(),j,i);
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
                                graph2D.GetVidFromElmt(locexp[i]->DetShapeType(),j,i);
                            
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
                                graph2D.GetVidFromElmt(locexp[i]->DetShapeType(),j,i);
                            
                            if(Btype == LibUtilities::eGLL_Lagrange)
                            {
                                for(k = 2; k < nedge_coeffs; ++k)
                                {
                                    m_locToContMap[cnt+vmap[nedge_coeffs+1-k]] = edge_offset[edgid]+(k-2);
                                }
                            }
                        }
                    }
                    
                }
                cnt += locexp[i]->GetNcoeffs();
            
            }
                
            gid = Vmath::Vmax(m_totLocLen,&m_locToContMap[0],1)+1;
            
            // setup interior mapping 
            for(i = 0; i < m_totLocLen; ++i)
            {
                if(m_locToContMap[i] == -1)
                {
                    m_locToContMap[i] = gid++;
                }
            }
            m_totGloLen = gid;
        }
    
        
        LocalToGlobalMap2D::~LocalToGlobalMap2D()
        {
        }
        
    }
}

/**
* $Log: LocalToGlobalMap2D.cpp,v $
* Revision 1.5  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
