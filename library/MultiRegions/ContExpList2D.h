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

#ifndef MULTIREGIONS_CONTEXPLIST2D_H
#define MULTIREGIONS_CONTEXPLIST2D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/LocalToGlobalMap2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	class ContExpList2D: 
	public ExpList2D 
	{
	public:
	    ContExpList2D();
            
	    ContExpList2D(const LibUtilities::BasisKey &TriBa, 
			  const LibUtilities::BasisKey &TriBb, 
			  const LibUtilities::BasisKey &QuadBa, 
			  const LibUtilities::BasisKey &QuadBb, 
			  const SpatialDomains::MeshGraph2D &graph2D,
                          const LibUtilities::PointsType 
                          TriNb = LibUtilities::SIZE_PointsType);

            ContExpList2D::ContExpList2D(const ContExpList2D &In);
	    
	    ~ContExpList2D();
	    
	    inline int getContNcoeffs()
	    {
		return m_contNcoeffs;
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
		m_locToGloMap->Assemble(m_coeffs,m_contCoeffs);
	    }
	   
	    void IProductWRTBase(const ExpList &In);
	    
	    void FwdTrans(const ExpList &In);
	    
	    void BwdTrans(const ExpList &In);
	    
	    void GenMassMatrixLinSys(void);
	    
	protected:
    
	    
	private:
	    int                    m_contNcoeffs;
	    Array<OneD, NekDouble> m_contCoeffs;
	    
	    boost::shared_ptr<LocalToGlobalMap2D> m_locToGloMap;
	    
	    DNekLinSysSharedPtr m_mass;
	    
	};

        typedef boost::shared_ptr<ContExpList2D>      ContExpList2DSharedPtr;
        typedef std::vector<ContExpList2DSharedPtr>   ContExpList2DVector;
        typedef std::vector<ContExpList2DSharedPtr>::iterator ContExpList2DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define

/**
* $Log: ContExpList2D.h,v $
**/
