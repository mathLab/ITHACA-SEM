///////////////////////////////////////////////////////////////////////////////
//
// File LocaltoGlobalMap1D.cpp
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
// Description: Local to Global mapping routines in 1D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalMap1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	LocalToGlobalMap1D::LocalToGlobalMap1D(int loclen, std::vector<StdRegions::StdExpansionVector> &exp_shapes, SpatialDomains::MeshGraph1D &graph1D)
	{
	    int i,j,gid,cnt;
	    int *locToContMap;
	    
	    m_totLocLen = loclen; 
	    
	    // set up Local to Continuous mapping 
	    locToContMap = new int [m_totLocLen];
	    Vmath::Fill(m_totLocLen,-1,locToContMap,1);
	  
	    StdRegions::StdExpMap vmap;
	    
	    // set up simple map based on vertex and edge id's
	    cnt = 0;
	    for(i = 0; i < exp_shapes[0].size(); ++i)
	    {
		exp_shapes[0][i]->MapTo(StdRegions::eForwards,vmap);
		for(j = 0; j < 2; ++j)
		{
		    locToContMap[cnt+vmap[j]] = graph1D.GetVidFromElmt(j,i);
		}
		cnt += exp_shapes[0][i]->GetNcoeffs();
	    }
	    
	    gid = Vmath::Vmax(m_totLocLen,locToContMap,1)+1;
	    for(i = 0; i < m_totLocLen; ++i)
	    {
		if(locToContMap[i] == -1)
		{
		    locToContMap[i] = gid++;
		}
	    }
	    
	    m_totGloLen = gid;
	    
	    m_locToContMap.reset(locToContMap);
	}
	
	
	LocalToGlobalMap1D::~LocalToGlobalMap1D()
	{
	}
	
    }
}
