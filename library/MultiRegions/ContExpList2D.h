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
                          TriNb = LibUtilities::SIZE_PointsType,
                          const bool constructMap = true);

            ContExpList2D(SpatialDomains::MeshGraph2D &graph2D,
                          const bool constructMap = true);
            
            ContExpList2D(const ContExpList2D &In);
            
            ~ContExpList2D();

            inline Array<OneD, NekDouble> &UpdateContCoeffs()
            {
                m_transState = eContinuous;
                return m_contCoeffs;
            }
        
            inline int getContNcoeffs()
            {
                return m_contNcoeffs;
            }
            
            inline void ContToLocal()
            {
                m_locToGloMap->ContToLocal(m_contCoeffs,m_coeffs);
            }

            inline void ContToLocal(const Array<OneD, const NekDouble> &inarray,
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

            inline void Assemble(const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->Assemble(inarray,outarray);
            }

            inline const LocalToGlobalMapSharedPtr& GetLocalToGlobalMap() const
            {
                return  m_locToGloMap;
            }
       
            void IProductWRTBase(const ExpList &In);
        
            void FwdTrans(const ExpList &In);
        
            void BwdTrans(const ExpList &In);

            void GeneralMatrixOp(const GlobalLinSysKey            &gkey,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble>          &outarray);
            
        protected:
            LocalToGlobalMapSharedPtr m_locToGloMap;
      	    int                       m_contNcoeffs;
	    Array<OneD, NekDouble>    m_contCoeffs;
            GlobalLinSysMapShPtr      m_globalMat;
                        
        private:
            
        };
        
        typedef boost::shared_ptr<ContExpList2D>      ContExpList2DSharedPtr;
        typedef std::vector<ContExpList2DSharedPtr>   ContExpList2DVector;
        typedef std::vector<ContExpList2DSharedPtr>::iterator ContExpList2DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define

/**
* $Log: ContExpList2D.h,v $
* Revision 1.10  2008/04/06 06:00:07  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.9  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.8  2008/01/20 16:31:11  bnelson
* Fixed linux compile errors.
*
* Revision 1.7  2007/12/17 13:05:04  sherwin
* Made files compatible with modifications in StdMatrixKey which now holds constants
*
* Revision 1.6  2007/12/06 22:52:29  pvos
* 2D Helmholtz solver updates
*
* Revision 1.5  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.4  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
