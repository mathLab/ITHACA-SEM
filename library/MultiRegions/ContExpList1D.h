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
#include <MultiRegions/GlobalLinSys.h>

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
                          const SpatialDomains::MeshGraph1D &graph1D,
                          const bool constructMap = true);
            ContExpList1D(SpatialDomains::MeshGraph1D &graph1D,
                          const bool constructMap = true);
            ContExpList1D(const ContExpList1D &In);
            ~ContExpList1D();

            inline Array<OneD, NekDouble> &UpdateContCoeffs()
            {
                m_transState = eContinuous;
                return m_contCoeffs;
            }

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

            void BwdTrans(const ExpList &In);

            void GeneralMatrixOp(const GlobalLinSysKey            &gkey,
                                 const ConstArray<OneD,NekDouble> &inarray,
                                 Array<OneD, NekDouble>           &outarray);
	    
	protected:
      	    int                       m_contNcoeffs;
	    Array<OneD, NekDouble>    m_contCoeffs;
	    LocalToGlobalMapSharedPtr m_locToGloMap;	
            GlobalLinSysMapShPtr      m_globalMat;
            
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
* Revision 1.27  2007/12/06 22:52:29  pvos
* 2D Helmholtz solver updates
*
* Revision 1.26  2007/11/07 20:29:52  jfrazier
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
* Revision 1.22  2007/09/03 19:58:30  jfrazier
* Formatting.
*
* Revision 1.21  2007/07/26 00:07:50  bnelson
* Fixed linux compiler errors.
*
* Revision 1.20  2007/07/23 16:06:30  sherwin
* Put a std::map to hold global matrix systems
*
* Revision 1.19  2007/07/22 23:04:19  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.18  2007/07/20 02:04:10  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.17  2007/07/19 20:02:25  sherwin
* Generalised global matrix solver
*
* Revision 1.16  2007/07/16 18:28:43  sherwin
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
