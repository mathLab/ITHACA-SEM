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

        ContExpList1D::ContExpList1D(const ContExpList1D &In):
            ExpList1D(In),
            m_contNcoeffs(In.m_contNcoeffs),
            m_locToGloMap(In.m_locToGloMap),
            m_mass(In.m_mass)

        {
            m_contCoeffs = Array<OneD,NekDouble>(m_contNcoeffs);
        }
        
	ContExpList1D::ContExpList1D(const LibUtilities::BasisKey &Ba, 
				     const SpatialDomains::MeshGraph1D &graph1D):
	    ExpList1D(Ba,graph1D)
	{
	    
	    ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)
		     ||(Ba.GetBasisType() == LibUtilities::eGLL_Lagrange),
		     "Expansion not of an boundary-interior type");
	    
	    // setup mapping array 
	    m_locToGloMap = MemoryManager<LocalToGlobalMap1D>::AllocateSharedPtr(m_ncoeffs,*m_exp,graph1D);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloLen();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs);
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
		GenMassMatrixLinSys();
	    }

	    DNekVec v(m_contNcoeffs,m_contCoeffs,eWrapper);
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
	
	
	void ContExpList1D::GenMassMatrixLinSys(void)
	{
	    if(!(m_mass.get()))
	    {
		int   i,j,n,gid1,gid2,loc_lda,cnt;
		DNekMatSharedPtr loc_mass;
		StdRegions::StdExpansionVectorIter def;
		
		DNekMatSharedPtr Gmass = MemoryManager<DNekMat>::AllocateSharedPtr(m_contNcoeffs,m_contNcoeffs);
		
		// fill global matrix 
		for(n = cnt = 0; n < (*m_exp).size(); ++n)
		{
		    loc_mass = (*m_exp)[n]->GetLocMatrix(StdRegions::eMassMatrix);
		    loc_lda = loc_mass->GetColumns();
		    
		    for(i = 0; i < loc_lda; ++i)
		    {
			gid1 = m_locToGloMap->GetMap(cnt + i);
			for(j = 0; j < loc_lda; ++j)
			{
			    gid2 = m_locToGloMap->GetMap(cnt + j);
			    (*Gmass)(gid1,gid2) += (*loc_mass)(i,j);
			}		
                    }
                    cnt += (*m_exp)[n]->GetNcoeffs();
		}
                m_mass = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmass);
	    }
	}
    } //end of namespace
} //end of namespace

/**
* $Log: ContExpList1D.cpp,v $
**/

