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
            m_globalMat(In.m_globalMat)

        {
            m_contCoeffs = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }
        

        ContExpList1D::ContExpList1D(const LibUtilities::BasisKey &Ba, 
                                     const SpatialDomains::Composite &cmps):
	    ExpList1D(Ba,cmps)
	{
	    
	    ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)
		     ||(Ba.GetBasisType() == LibUtilities::eGLL_Lagrange),
		     "Expansion not of an boundary-interior type");
	    

            // setup Matrix Map
            m_globalMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

	    // setup mapping array 
	    m_locToGloMap = MemoryManager<LocalToGlobalMap1D>::AllocateSharedPtr(m_ncoeffs,*m_exp,cmps);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloLen();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
	}
              
	void ContExpList1D::IProductWRTBase(const ExpList &In)
	{
            if(m_transState == eContinuous)
            {
                ContToLocal();
            }
	    ExpList1D::IProductWRTBase(In);
	    Assemble();
	    m_transState = eLocalCont;
	}

	void ContExpList1D::GeneralMatrixOp(const StdRegions::MatrixType mtype,
                                            const ConstArray<OneD, NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray, 
                                            NekDouble lambda = 1.0)

	{
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(m_ncoeffs);
            ContToLocal(inarray,tmp);
	    ExpList1D::GeneralMatrixOp(mtype,tmp,tmp,lambda);
	    Assemble(tmp,outarray);
	}
	
	void ContExpList1D::FwdTrans(const ExpList &In)
	{
            IProductWRTBase(In);

            GlobalLinSysSharedPtr mass_matrix;
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalLinSysMap::iterator matrixIter = m_globalMat->find(key);
           
            if(matrixIter == m_globalMat->end())
            {
                mass_matrix = GenGlobalLinSys(key,0);
                (*m_globalMat)[key] = mass_matrix;
            }
            else
            {
                mass_matrix = matrixIter->second;
            }

	    DNekVec v(m_contNcoeffs,m_contCoeffs,eWrapper);
	    mass_matrix->GetLinSys()->Solve(v,v);
	    m_transState = eContinuous;
	    m_physState = false;
	}

	void ContExpList1D::HelmSolve(const ExpList &In, NekDouble lambda)
	{
            IProductWRTBase(In);
            Vmath::Neg(m_contNcoeffs,&m_contCoeffs[0],1);

            GlobalLinSysSharedPtr helm_matrix;
            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            GlobalLinSysMap::iterator matrixIter = m_globalMat->find(key);
           
            if(matrixIter == m_globalMat->end())
            {
                helm_matrix = GenGlobalLinSys(key,0);
                (*m_globalMat)[key] = helm_matrix;
            }
            else
            {
                helm_matrix = matrixIter->second;
            }


	    DNekVec v(m_contNcoeffs,m_contCoeffs,eWrapper);
	    helm_matrix->GetLinSys()->Solve(v,v);
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
	
	
	GlobalLinSysSharedPtr ContExpList1D::GenGlobalLinSys(const GlobalLinSysKey &mkey, const int NumDirBCs)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt;
            DNekScalMatSharedPtr loc_mat;
            StdRegions::StdExpansionVectorIter def;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;
            
            unsigned int rows = m_contNcoeffs - NumDirBCs;
            unsigned int cols = m_contNcoeffs - NumDirBCs;
            NekDouble zero = 0.0;
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);
            
            // fill global matrix 
            for(n = cnt = 0; n < (*m_exp).size(); ++n)
            {
                LocalRegions::MatrixKey matkey(mkey.GetLinSysType(),
                                          (*m_exp)[n]->DetShapeType(),
                                         *(*m_exp)[n],mkey.GetScaleFactor());
                
                loc_mat = (*m_exp)[n]->GetLocMatrix(matkey);
                loc_lda = loc_mat->GetColumns();
		    
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = m_locToGloMap->GetMap(cnt + i);
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = m_locToGloMap->GetMap(cnt + j);
                            if(gid2 >= NumDirBCs)
                            {
                                (*Gmat)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                    += (*loc_mat)(i,j);
                            }
                        }		
                    }
                }
                cnt += (*m_exp)[n]->GetNcoeffs();
            }
            
            linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
            
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }

    } //end of namespace
} //end of namespace

/**
* $Log: ContExpList1D.cpp,v $
* Revision 1.20  2007/07/29 07:45:30  sherwin
* Updated for new memory manager call
*
* Revision 1.19  2007/07/27 03:10:48  bnelson
* Fixed g++ compile error.
*
* Revision 1.18  2007/07/23 16:06:30  sherwin
* Put a std::map to hold global matrix systems
*
* Revision 1.17  2007/07/19 20:02:24  sherwin
* Generalised global matrix solver
*
* Revision 1.16  2007/07/16 18:28:42  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.15  2007/07/13 15:22:12  sherwin
* Update for Helmholtz (working without bcs )
*
* Revision 1.14  2007/07/13 09:02:23  sherwin
* Mods for Helmholtz solver
*
* Revision 1.13  2007/07/10 08:54:29  pvos
* Updated ContField1D constructor
*
* Revision 1.12  2007/07/06 18:39:33  pvos
* ContField1D constructor updates
*
* Revision 1.11  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

