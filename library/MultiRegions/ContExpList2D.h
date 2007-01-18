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

#ifndef MULTIREGIONS_CONTEXPLIST1D_H
#define MULTIREGIONS_CONTEXPLIST1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/LocalToGlobalMap2D.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.h>
#include <StdRegions/StdMatrix.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	class ContExpList2D: 
	public ExpList2D 
	{
	public:
	    ContExpList2D();
	    ContExpList2D(const StdRegions::BasisKey &TriBa, 
		const StdRegions::BasisKey &TriBb, 
		const StdRegions::BasisKey &QuadBa, 
		const StdRegions::BasisKey &QuadBb, 
		SpatialDomains::MeshGraph2D &graph2D);

	    ContExpList2D(const StdRegions::BasisKey &TriBa, 
			  const StdRegions::BasisKey &TriBb, 
			  const StdRegions::NodalBasisType TriNb,
			  const StdRegions::BasisKey &QuadBa, 
			  const StdRegions::BasisKey &QuadBb, 
			  SpatialDomains::MeshGraph2D &graph2D);
	    
	    ~ContExpList2D();
	    
	    inline int getContNcoeffs()
	    {
		return m_contNcoeffs;
	    }

	    inline double *getContCoeffs()
	    {
		return m_contCoeffs;
	    }
	    
	    inline void ContToLocal()
	    {
		m_locToGloMap->ContToLocal(m_contCoeffs,m_coeffs);
	    }
	    
	    inline void LocalToCont()
	    {
		m_locToGloMap->LocalToCont(m_coeffs,m_contCoeffs);
	    }
	    
	    
	    inline void Assemble()
	    {
		Assemble(m_coeffs,m_contCoeffs);
	    }
	    
	    inline void Assemble(double *loc, double *cont)
	    {
		m_locToGloMap->Assemble(loc,cont);
	    }
	    
	    void IProductWRTBase(const double *inarray, double *outarray);
      
	    void FwdTrans(const double *inarray);
	    
	    void BwdTrans(double *outarray);
	    
	    void GenMassMatrix(void);
	    
	protected:
    
	    
	private:
	    int     m_contNcoeffs;
	    double *m_contCoeffs;
	    
	    boost::shared_ptr<LocalToGlobalMap2D> m_locToGloMap;

#ifdef NEKMAT
	    boost::shared_ptr< DNekMat > m_mass;
#else
	    StdRegions::StdMatContainer *m_mass;
#endif
	    
	    virtual void v_BwdTrans(double *outarray)
	    {
		BwdTrans(outarray);
	    }
	    
	};
    } //end of namespace
} //end of namespace

#endif // end of define
