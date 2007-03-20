///////////////////////////////////////////////////////////////////////////////
//
// File ContExpList1D.cpp
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
// Description: Continuous Expansion list definition in 1D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContExpList1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	ContExpList1D::ContExpList1D()
	{
	}
	
	ContExpList1D::~ContExpList1D()
	{
	}

	ContExpList1D::ContExpList1D(const LibUtilities::BasisKey &Ba, 
				     SpatialDomains::MeshGraph1D &graph1D):
	    ExpList1D(Ba,graph1D)
	{
	    
	    ASSERTL1((Ba.GetBasisType() == LibUtilities::Modified_A)
		     ||(Ba.GetBasisType() == LibUtilities::GLL_Lagrange),
		     "Expansion not of an boundary-interior type");
	    
	    // setup mapping array 
	    m_locToGloMap.reset(new LocalToGlobalMap1D(m_exp,graph1D));
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloLen();
	    m_contCoeffs  = MemoryManager::AllocateSharedArray<double> 
		(m_contNcoeffs);
	}
              
	void ContExpList1D::IProductWRTBase(const ExpList &In)
	{
	    ExpList1D::IProductWRTBase(In);
	    Assemble();
	    m_transState = eLocalCont;
	}
	
	void ContExpList1D::FwdTrans(const ExpList &In)
	{
	    IProductWRTBase(In);

	    if(!(m_mass.get()))
	    {
		GenMassMatrix();
	    }

	    DNekVec v(m_ncoeffs,m_contCoeffs,eWrapper);
	    m_mass->Solve(v,v);
	    m_transState = eContinuous;
	    m_physState = false;
	}
	
	void ContExpList1D::BwdTrans(const ExpList &In)
	{
	    
	    if(m_transState == eContinuous)
	    {
		ContToLocal();
	    }
	    
	    ExpList1D::BwdTrans(In);
	}
	
	
	void ContExpList1D::GenMassMatrix(void)
	{
	    if(!(m_mass.get()))
	    {
		int   i,j,n,gid1,gid2,loc_lda;
		double *loc_mat;
		DNekMatSharedPtr loc_mass;
		StdRegions::StdExpansionVectorIter def;
		
		DNekLinSys Gmass(m_contNcoeffs,m_contNcoeffs);
		m_mass.reset(Gmass);
		
		// fill global matrix 
		for(n = 0; n < m_exp.size(); ++n)
		{
		    loc_mass = GetExp(n)->GetLocMass(StdRegions::eMassMatrix);
		    loc_lda = loc_mass->GetColumns();
		    loc_mat = loc_mass->GetMatrix();
		    
		    for(i = 0; i < loc_lda; ++i)
		    {
			gid1 = m_locToGloMap->GetMap(n,i);
			
			for(j = 0; j < loc_lda; ++j)
			{
			    gid2 = m_locToGloMap->GetMap(n,j);
			    Gmass(gid1*m_contNcoeffs,gid2) += loc_mat[i*loc_lda+j];
			}
		    }
		}
	    }
	}
    } //end of namespace
} //end of namespace



