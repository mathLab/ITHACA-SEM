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
// Description: Continusou Expansion list definition in 1D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_CONTEXPLIST1D_H
#define NEKTAR_LIB_MULTIREGIONS_CONTEXPLIST1D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/LocalToGlobalMap1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	class ContExpList1D: 
	public ExpList1D 
	{
	public:
	    ContExpList1D();
	    ContExpList1D(const LibUtilities::BasisKey &Ba, 
			  const SpatialDomains::MeshGraph1D &graph1D);
            ContExpList1D::ContExpList1D(const ContExpList1D &In);
	    ~ContExpList1D();
	    
	    inline int GetContNcoeffs()
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
	    
	    boost::shared_ptr<LocalToGlobalMap1D> m_locToGloMap;
	    
	    DNekLinSysSharedPtr m_mass;
	    
	};
	
        typedef boost::shared_ptr<ContExpList1D>      ContExpList1DSharedPtr;
        typedef std::vector<ContExpList1DSharedPtr>   ContExpList1DVector;
        typedef std::vector<ContExpList1DSharedPtr>::iterator ContExpList1DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define
