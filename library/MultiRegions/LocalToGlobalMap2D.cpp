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
	LocalToGlobalMap2D::LocalToGlobalMap2D(int loclen, 
		       std::vector<StdRegions::StdExpansionVector> &exp_shapes, 
		       SpatialDomains::MeshGraph2D &graph2D)
	{
	    int i,j,k,gid,cnt, nedge, nedge_coeffs,edgid;
	    int nGloVerts = 0;
	    int *locToContMap;
	    std::vector<StdRegions::StdExpansionVector>::iterator  def;
	    StdRegions::EdgeOrientation eorient;
	    StdRegions::BasisType Btype;
	   
	    m_totLocLen = loclen; 

	    locToContMap = new int [m_totLocLen];
	    Vmath::Fill(m_totLocLen,-1,locToContMap,1);
	    
	    // set up mapping based 
	    StdRegions::StdExpMap vmap;
    
	    // need to put into LocalToGloabl initialiser. 
	    int *edge_offset  = new int [graph2D.GetNecomps()+1];
	    double *sign;
	    m_sign_change = false;

	    // determine global vertex ids
	    for(def = exp_shapes.begin(); def != exp_shapes.end(); ++def)
	    {
		for(i = 0; i < (*def).size(); ++i)
		{
		    for(j = 0; j < (*(*def)[i]).GetNverts(); ++j)
		    {
			
			nGloVerts = std::max(nGloVerts,
					     graph2D.GetVidFromElmt((*(*def)[i]).DetShapeType(),j,i));
			//set up nedge coefficients for each edge 	
			nedge_coeffs = (*(*def)[i]).GetEdgeNcoeffs(j);
			edge_offset[graph2D.GetEidFromElmt((*(*def)[i]).DetShapeType(),
							   j,i)+1] =   nedge_coeffs -2;
			
			// need a sign vector if edge_nceoff >=4 
			if((nedge_coeffs >= 4)&&((*(*def)[i]).GetEdgeBasisType(0) == 
						 StdRegions::eModified_A))
			{
			    m_sign_change = true;
			}
		    }
		}
	    }
	    nGloVerts++; // set actual value to one plus maximum id
	    
	    // set up sign vector
	    if(m_sign_change){
		sign = new double[m_totLocLen];
		Vmath::Fill(m_totLocLen,1.0,sign,1);
		m_sign.reset(sign);
	    }

	    edge_offset[0] = nGloVerts; 

	    // set up consecutive list for edge entries
	    for(i = 1; i < graph2D.GetNecomps(); ++i)
	    {
		edge_offset[i] += edge_offset[i-1];
	    }
	    

	    // set up simple map;
	    cnt = 0;
	    for(def = exp_shapes.begin(); def != exp_shapes.end(); ++def)
	    {
		for(i = 0; i < (*def).size(); ++i)
		{
		    for(j = 0; j < (nedge = (*(*def)[i]).GetNedges()); ++j)
		    {
			(*(*def)[i]).MapTo_ModalFormat(
			 nedge_coeffs = (*(*def)[i]).GetEdgeNcoeffs(j),
			 Btype = (*(*def)[i]).GetEdgeBasisType(j), j,
			 eorient = graph2D.GetEorientFromElmt((*(*def)[i]).DetShapeType(),j,i), vmap);


			// set edge ids before setting vertices 
			// so that can be reverse for nodal expansion when j>2
			
			edgid = graph2D.GetEidFromElmt((*(*def)[i]).DetShapeType(),j,i);
			
			for(k = 2; k < nedge_coeffs; ++k)
			{
			    locToContMap[cnt+vmap[k]] =  edge_offset[edgid]+(k-2);
			}
			
			// vmap is set up according to cartesian
			// coordinates so have to change setting
			// depending on which edge we are considering
			if(j < 2)
			{
			    if(eorient == StdRegions::eForwards)
			    {
				locToContMap[cnt+vmap[0]] = 
				graph2D.GetVidFromElmt((*(*def)[i]).DetShapeType(),j,i);
			    }
			    else
			    {
				locToContMap[cnt+vmap[1]] = 
				graph2D.GetVidFromElmt((*(*def)[i]).DetShapeType(),j,i);
				if(m_sign_change)
				{
				    for(k = 3; k < nedge_coeffs; k+=2)
				    {
					sign[cnt+vmap[k]] = -1;
				    }
				}
				
			    }
			    
			}
			else
			{
			    // set edge global ids
			    if(eorient == StdRegions::eForwards)
			    {
				locToContMap[cnt+vmap[1]] = 
				graph2D.GetVidFromElmt((*(*def)[i]).DetShapeType(),j,i);

				if(m_sign_change)
				{
				    for(k = 3; k < nedge_coeffs; k+=2)
				    {
					sign[cnt+vmap[k]] = -1;
				    }
				}
			    }
			    else
			    {
				locToContMap[cnt+vmap[0]] = 
				graph2D.GetVidFromElmt((*(*def)[i]).DetShapeType(),j,i);

				if(Btype == StdRegions::eGLL_Lagrange)
				{
				    for(k = 2; k < nedge_coeffs; ++k)
				    {
					locToContMap[cnt+vmap[nedge_coeffs+1-k]] = edge_offset[edgid]+(k-2);
				    }
				}
			    }
			}

		    }
		    cnt += (*(*def)[i]).GetNcoeffs();
		}
	    }

	    gid = Vmath::Vmax(m_totLocLen,locToContMap,1)+1;
	    
	    // setup interior mapping 
	    for(i = 0; i < m_totLocLen; ++i)
	    {
		
		if(locToContMap[i] == -1)
		{
		    locToContMap[i] = gid++;
		}
	    }
	    
	  m_totGloLen = gid;
	  
	  m_locToContMap.reset(locToContMap);
	  
	  delete[] edge_offset;
	  


	}
	

	LocalToGlobalMap2D::~LocalToGlobalMap2D()
	{
	}
      
  }
}
