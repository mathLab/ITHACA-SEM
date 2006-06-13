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
      m_mass = (StdRegions::StdMatContainer *)NULL;
      m_locToContMap = (int *)NULL;
  }
  
    ContExpList1D::~ContExpList1D()
    {
      if(m_contCoeffs)
      {
	delete[] m_contCoeffs;
      }
      
      if(m_locToContMap)
      {
	delete[] m_locToContMap;
      }
    }

    ContExpList1D::ContExpList1D(const StdRegions::BasisKey &Ba, 
				 SpatialDomains::MeshGraph1D &graph1D):
	ExpList1D(Ba,graph1D)
    {
      int gid,i,j;
      int order = Ba.GetBasisOrder();
      
      ASSERTL1((Ba.GetBasisType() == StdRegions::Modified_A)
	       ||(Ba.GetBasisType() == StdRegions::GLL_Lagrange),
	       "Expansion not of an boundary-interior type");

      m_mass = (StdRegions::StdMatContainer *)NULL;
      
      m_locToContMap = new int [m_ncoeffs];
      Vmath::Fill(m_ncoeffs,-1,m_locToContMap,1);
      
      // set up mapping based 
      StdRegions::StdExpMap vmap;
    
      // assume all elements have the same mapping and expasion order
      static_cast<StdRegions::StdSegExp*>(m_exp_shapes[0][0].get())->MapTo(StdRegions::eForwards, vmap);
   
      // set up simple map;
      for(gid = i = 0; i < m_exp_shapes[0].size(); ++i,++gid)
      {
	for(j = 0; j < 2; ++j)
	{
	  m_locToContMap[order*i+vmap[j]] = gid+j;
	}
      }
      ++gid;

      for(i = 0; i < m_exp_shapes[0].size(); ++i)
      {
	for(j = 0; j < order; ++j)
	{
	  if(m_locToContMap[order*i+j] == -1)
	  {
	    m_locToContMap[order*i+j] = gid++;
	  }
	}
      }

      m_contNcoeffs = gid;
      m_contCoeffs = new double [m_contNcoeffs];
    }
    
    
    void ContExpList1D::IProductWRTBase(const double *inarray, double *outarray)
    {
      ExpList1D::IProductWRTBase(inarray,m_coeffs);
      Assemble(m_coeffs,outarray);
      m_transState = eLocalCont;
    }

    void ContExpList1D::FwdTrans(const double *inarray)
    {
      IProductWRTBase(inarray,m_contCoeffs);
      if(!m_mass)
      {
	GenMassMatrix();
      }
      m_mass->Solve(m_contCoeffs,1);
      m_transState = eContinuous;
      m_physState = false;
    }

    void ContExpList1D::BwdTrans(double *outarray)
    {

      if(m_transState == eContinuous)
      {
	ContToLocal();
      }

      ExpList1D::BwdTrans(outarray);
    }
    

    void ContExpList1D::GenMassMatrix(void)
    {
      if(!m_mass)
      {
	int   i,j,cnt,gid1,gid2,loc_lda;
	double *loc_mat;
	StdRegions::StdMatContainer *loc_mass;
	StdRegions::StdExpansionVectorIter def;
	
	double *mmat = new double [m_contNcoeffs*m_contNcoeffs];
	Vmath::Zero(m_contNcoeffs*m_contNcoeffs,mmat,1);
      
	m_mass = new StdRegions::StdMatContainer(mmat);
	m_mass->SetLda     (m_contNcoeffs);
	m_mass->SetMatForm (StdRegions::eSymmetric_Positive);
	
      // fill global matrix 
	for(cnt = 0, def = m_exp_shapes[0].begin(); 
	    def != m_exp_shapes[0].end(); ++def)
	{
	  loc_mass = (*def)->GetMassMatrix();
	  loc_lda = loc_mass->GetLda();
	  loc_mat = loc_mass->GetMatrix();
	
	  for(i = 0; i < loc_lda; ++i)
	  {
	    gid1 = m_locToContMap[i+cnt];

	    for(j = 0; j < loc_lda; ++j)
	    {
	      gid2 = m_locToContMap[j+cnt];
	      mmat[gid1*m_contNcoeffs + gid2] += loc_mat[i*loc_lda + j];
	    }
	  }
	  cnt+=(*def)->GetNcoeffs();
	}
      }
    }
  } //end of namespace
} //end of namespace

