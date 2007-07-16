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
            //ContExpList1D(const LibUtilities::BasisKey &Ba, 
            //const SpatialDomains::MeshGraph1D &graph1D);
            ContExpList1D(const LibUtilities::BasisKey &Ba, 
                          const SpatialDomains::Composite &cmps);
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
	    
	    inline void ContToLocal(const ConstArray<OneD,NekDouble> &inarray,
                                    Array<OneD,NekDouble> &outarray)
	    {
		m_locToGloMap->ContToLocal(inarray,outarray);
	    }
	    
	    inline void LocalToCont()
	    {
		m_locToGloMap->LocalToCont(m_coeffs,m_contCoeffs);
	    }
	    	    
	    inline void Assemble()
	    {
		m_locToGloMap->Assemble(m_coeffs,m_contCoeffs);
	    }
	    
	    inline void Assemble(const ConstArray<OneD,NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray)
	    {
		m_locToGloMap->Assemble(inarray,outarray);
	    }
	    
	    void IProductWRTBase(const ExpList &In);
	    
	    void FwdTrans(const ExpList &In);

	    void HelmSolve(const ExpList &In, NekDouble lambda);
	    
	    void BwdTrans(const ExpList &In);

            void GeneralMatrixOp(const StdRegions::MatrixType     mtype,
                                 const ConstArray<OneD,NekDouble> &inarray,
                                 Array<OneD, NekDouble>          &outarray,
                                 NekDouble lambda);
	    
	    void GenMassMatrixLinSys(void);
	    void GenHelmholtzMatrixLinSys(NekDouble lambda);
	    
	protected:
      	    int                    m_contNcoeffs;
	    Array<OneD, NekDouble> m_contCoeffs;
	    
	    boost::shared_ptr<LocalToGlobalMap1D> m_locToGloMap;
	    
	    DNekLinSysSharedPtr m_mass;
	    DNekLinSysSharedPtr m_helm;
	    
	private:

	};
	
        typedef boost::shared_ptr<ContExpList1D>      ContExpList1DSharedPtr;
        typedef std::vector<ContExpList1DSharedPtr>   ContExpList1DVector;
        typedef std::vector<ContExpList1DSharedPtr>::iterator ContExpList1DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define

/**
* $Log: ContExpList1D.h,v $
* Revision 1.15  2007/07/13 15:22:12  sherwin
* Update for Helmholtz (working without bcs )
*
* Revision 1.14  2007/07/13 09:02:23  sherwin
* Mods for Helmholtz solver
*
* Revision 1.13  2007/07/10 08:54:29  pvos
* Updated ContField1D constructor
*
* Revision 1.12  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.11  2007/06/08 12:58:26  sherwin
* Added ContField1D and remove previous structure using Fields
*
* Revision 1.10  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
