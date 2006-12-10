///////////////////////////////////////////////////////////////////////////////
//
// File ContExpList2D.cpp
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
// Description: Continuous Expansion list definition in 2D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContExpList2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	ContExpList2D::ContExpList2D()
	{
	    m_mass = (StdRegions::StdMatContainer *)NULL;
	}
	
	ContExpList2D::~ContExpList2D()
	{
	    if(m_contCoeffs)
	    {
		delete[] m_contCoeffs;
	    }
	}
	
	ContExpList2D::ContExpList2D(const StdRegions::BasisKey &TriBa, 
				     const StdRegions::BasisKey &TriBb, 
				     const StdRegions::BasisKey &QuadBa, 
				     const StdRegions::BasisKey &QuadBb, 
				     SpatialDomains::MeshGraph2D &graph2D):
	    ExpList2D(TriBa,TriBb,QuadBa,QuadBb,graph2D)
	{
	    
	    ASSERTL1((TriBa.GetBasisType() == StdRegions::eModified_A),
		     "Expansion not of an boundary-interior type");
	    
	    ASSERTL1((TriBb.GetBasisType() == StdRegions::eModified_B),
		     "Expansion not of an boundary-interior type");

	    ASSERTL1((QuadBa.GetBasisType() == StdRegions::eModified_A)
		     ||(QuadBa.GetBasisType() == StdRegions::eGLL_Lagrange),
		     "Expansion not of an boundary-interior type");

	    ASSERTL1((QuadBb.GetBasisType() == StdRegions::eModified_A)
		     ||(QuadBb.GetBasisType() == StdRegions::eGLL_Lagrange),
		     "Expansion not of an boundary-interior type");
	    
	    if(TriBa.GetBasisType() == StdRegions::eModified_A)
	    {
		ASSERTL1((QuadBb.GetBasisType() == StdRegions::eModified_A),
			 "Quad and Tri Expansions are not of the same type");
	    }
	    
	    m_mass = (StdRegions::StdMatContainer *)NULL;
	    
	    // setup mapping array 
	    m_locToGloMap.reset(new LocalToGlobalMap2D(m_ncoeffs, m_exp_shapes,
							graph2D));
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloLen();
	    m_contCoeffs = new double [m_contNcoeffs];
	}


	ContExpList2D::ContExpList2D(const StdRegions::BasisKey &TriBa, 
				     const StdRegions::BasisKey &TriBb, 
				     const StdRegions::NodalBasisType &TriNb,
				     const StdRegions::BasisKey &QuadBa, 
				     const StdRegions::BasisKey &QuadBb, 
				     SpatialDomains::MeshGraph2D &graph2D):
	    ExpList2D(TriBa,TriBb,TriNb,QuadBa,QuadBb,graph2D)
	{
	    
	    ASSERTL1((TriBa.GetBasisType() == StdRegions::eModified_A)&&
		     (TriNb == NULL),
		     "Expansion not of an boundary-interior type");
	    
	    ASSERTL1((TriBb.GetBasisType() == StdRegions::eModified_B)&&
		     (TriNb == NULL),
		     "Expansion not of an boundary-interior type");

	    ASSERTL1((QuadBa.GetBasisType() == StdRegions::eModified_A)
		     ||(QuadBa.GetBasisType() == StdRegions::eGLL_Lagrange),
		     "Expansion not of an boundary-interior type");

	    ASSERTL1((QuadBb.GetBasisType() == StdRegions::eModified_A)
		     ||(QuadBb.GetBasisType() == StdRegions::eGLL_Lagrange),
		     "Expansion not of an boundary-interior type");
	    
	    if(TriBa.GetBasisType() == StdRegions::eModified_A)
	    {
		ASSERTL1((QuadBb.GetBasisType() == StdRegions::eModified_A),
			 "Quad and Tri Expansions are not of the same type");
	    }
	    
	    m_mass = (StdRegions::StdMatContainer *)NULL;
	    
	    // setup mapping array 
	    m_locToGloMap.reset(new LocalToGlobalMap2D(m_ncoeffs, m_exp_shapes,
							graph2D));
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloLen();
	    m_contCoeffs = new double [m_contNcoeffs];
	}
	
	void ContExpList2D::IProductWRTBase(const double *inarray, 
					    double *outarray)
	{
	    ExpList2D::IProductWRTBase(inarray,m_coeffs);
	    Assemble(m_coeffs,outarray);
	    m_transState = eLocalCont;
	}

	void ContExpList2D::FwdTrans(const double *inarray)
	{
	    IProductWRTBase(inarray,m_contCoeffs);
	    if(!m_mass)
	    {
		GenMassMatrix();
	    }
	    
	    m_mass->ShowMatrixStructure(stdout);
	    m_mass->Solve(m_contCoeffs,1);
	    m_transState = eContinuous;
	    m_physState = false;
	}

	void ContExpList2D::BwdTrans(double *outarray)
	{
	    
	    if(m_transState == eContinuous)
	    {
		ContToLocal();
	    }
	    
	    ExpList2D::BwdTrans(outarray);
	}
	

	void ContExpList2D::GenMassMatrix(void)
	{
	    if(!m_mass)
	    {
		int   i,j,cnt,gid1,gid2,loc_lda;
		double *loc_mat,sign1,sign2;
		StdRegions::StdMatContainer *loc_mass;
		std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
		StdRegions::StdExpansionVectorIter def;
	
		double *mmat = new double [m_contNcoeffs*m_contNcoeffs];
		Vmath::Zero(m_contNcoeffs*m_contNcoeffs,mmat,1);
		
		m_mass = new StdRegions::StdMatContainer(mmat);
		m_mass->SetLda     (m_contNcoeffs);
		m_mass->SetMatForm (StdRegions::eSymmetric_Positive);
		
		// fill global matrix 
		for(cnt = 0, sdef = m_exp_shapes.begin(); sdef != m_exp_shapes.end();
		    ++sdef)
		{
		    
		    for(def = (*sdef).begin(); def != (*sdef).end(); ++def)
		    {

			loc_mass = (*def)->GetMassMatrix();
			loc_lda = loc_mass->GetLda();
			loc_mat = loc_mass->GetMatrix();
		    
			for(i = 0; i < loc_lda; ++i)
			{
			    gid1  = m_locToGloMap->GetMap(i+cnt);
			    sign1 = m_locToGloMap->GetSign(i+cnt);
			    
			    for(j = 0; j < loc_lda; ++j)
			    {
				gid2 = m_locToGloMap->GetMap(j+cnt);	
				sign2 = m_locToGloMap->GetSign(j+cnt);
				
				mmat[gid1*m_contNcoeffs + gid2] 
				+= sign1*sign2*loc_mat[i*loc_lda + j];
			    }
			}
			cnt+=(*def)->GetNcoeffs();
		    }
		}
	    }
	}
    } //end of namespace
} //end of namespace

