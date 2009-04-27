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

        ContExpList1D::ContExpList1D():
            ExpList1D(),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalLinSys()
        {
        }

        ContExpList1D::~ContExpList1D()
        {
        }

        ContExpList1D::ContExpList1D(const ContExpList1D &In):
            ExpList1D(In),
            m_locToGloMap(In.m_locToGloMap),
            m_contNcoeffs(In.m_contNcoeffs),
            m_contCoeffs(m_contNcoeffs,0.0),
            m_globalLinSys(In.m_globalLinSys)
        {
        }

        ContExpList1D::ContExpList1D(const LibUtilities::BasisKey &Ba,
                                     const SpatialDomains::MeshGraph1D &graph1D,
                                     const bool constructMap):
	    ExpList1D(Ba,graph1D),
            m_locToGloMap(),
            m_contNcoeffs(),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {   
            ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)
                     ||(Ba.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
            
	    // setup mapping array 
            if(constructMap)
            {
                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp);
                m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
                m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
            }
	}
        
        ContExpList1D::ContExpList1D(SpatialDomains::MeshGraph1D &graph1D,
                                     const bool constructMap):
	    ExpList1D(graph1D),
            m_locToGloMap(),
            m_contNcoeffs(),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            const SpatialDomains::ExpansionVector &expansions = graph1D.GetExpansions();
            
            for(int i = 0; i < expansions.size(); ++i)
            {	    
                ASSERTL1((expansions[i]->m_BasisKeyVector[0].GetBasisType() == LibUtilities::eModified_A)
                          ||(expansions[i]->m_BasisKeyVector[0].GetBasisType() == LibUtilities::eGLL_Lagrange),
                         "Expansion not of an boundary-interior type");
            }
                        
	    // setup mapping array 
            if(constructMap)
            {
                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp);
                m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
                m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
            }
	}
        
        void ContExpList1D::IProductWRTBase(const Array<OneD, const NekDouble> &inarray,  
                                                  Array<OneD,       NekDouble> &outarray,
                                            bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD, NekDouble> wsp(m_ncoeffs);
                IProductWRTBase_IterPerExp(inarray,wsp);
                Assemble(wsp,outarray);
            }
            else
            {
                IProductWRTBase_IterPerExp(inarray,outarray);
            }
        }

        void ContExpList1D::GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                            const Array<OneD,const NekDouble> &inarray, 
                                                  Array<OneD,      NekDouble> &outarray,
                                            bool  UseContCoeffs)

        {
            if(UseContCoeffs)
            {
                Array<OneD,NekDouble> tmp1(2*m_ncoeffs);
                Array<OneD,NekDouble> tmp2(tmp1+m_ncoeffs);
                GlobalToLocal(inarray,tmp1);
                GeneralMatrixOp_IterPerExp(gkey,tmp1,tmp2);
                Assemble(tmp2,outarray);
            }
            else
            {
                GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
            }
        }

        void ContExpList1D::FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD,       NekDouble> &outarray,
                                     bool  UseContCoeffs)
        {
            GlobalLinSysSharedPtr mass_matrix;
            GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
            GlobalLinSysMap::iterator matrixIter = m_globalLinSys->find(key);
            
            if(matrixIter == m_globalLinSys->end())
            {
                mass_matrix = GenGlobalLinSys(key,m_locToGloMap);
                (*m_globalLinSys)[key] = mass_matrix;
            }
            else
            {
                mass_matrix = matrixIter->second;
            }
            
            IProductWRTBase(inarray,outarray,true);

            if(UseContCoeffs)
            {                
                mass_matrix->Solve(outarray,outarray,*m_locToGloMap);
            }
            else
            {
                Array<OneD, NekDouble> wsp(m_contNcoeffs);
                mass_matrix->Solve(outarray,wsp,*m_locToGloMap);
                GlobalToLocal(wsp,outarray);
            }

        }
        
        void ContExpList1D::BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD,       NekDouble> &outarray,
                                     bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD, NekDouble> wsp(m_ncoeffs);
                GlobalToLocal(inarray,wsp);
                BwdTrans_IterPerExp(wsp,outarray);
            }
            else
            {
                BwdTrans_IterPerExp(inarray,outarray);
            }
        }

    } //end of namespace
} //end of namespace

/**
* $Log: ContExpList1D.cpp,v $
* Revision 1.37  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.36  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.35  2009/02/08 09:06:19  sherwin
* Added updated definition of GlobalLinSysKey to include localtoglobalbasemap
*
* Revision 1.34  2009/01/06 21:05:56  sherwin
* Added virtual function calls for BwdTrans, FwdTrans and IProductWRTBase from the class ExpList. Introduced _IterPerExp versions of these methods in ExpList.cpp§
*
* Revision 1.33  2008/09/16 13:36:05  pvos
* Restructured the LocalToGlobalMap classes
*
* Revision 1.32  2008/07/10 13:02:33  pvos
* Added periodic boundary conditions functionality
*
* Revision 1.31  2008/04/06 06:00:06  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.30  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.29  2007/12/17 13:05:03  sherwin
* Made files compatible with modifications in StdMatrixKey which now holds constants
*
* Revision 1.28  2007/12/06 22:52:29  pvos
* 2D Helmholtz solver updates
*
* Revision 1.27  2007/11/20 16:27:15  sherwin
* Zero Dirichlet version of UDG Helmholtz solver
*
* Revision 1.26  2007/11/07 20:29:48  jfrazier
* Modified to use new expansion list contained in meshgraph.
*
* Revision 1.25  2007/10/04 12:10:04  sherwin
* Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
*
* Revision 1.24  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.23  2007/09/25 14:25:29  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.22  2007/08/11 23:43:25  sherwin
* Expansion bases reader part for Helmholtz1D
*
* Revision 1.21  2007/08/02 12:36:18  sherwin
* ...
*
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
\*
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

