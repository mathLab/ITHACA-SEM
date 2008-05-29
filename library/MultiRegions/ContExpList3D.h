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

#ifndef MULTIREGIONS_CONTEXPLIST3D_H
#define MULTIREGIONS_CONTEXPLIST3D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/LocalToGlobalMap3D.h>

namespace Nektar
{
    namespace MultiRegions
    {
    
    class ContExpList3D: public ExpList3D
        {
        public:
            ContExpList3D();

            ContExpList3D(const LibUtilities::BasisKey &Ba,
                          const LibUtilities::BasisKey &Bb,
                          const LibUtilities::BasisKey &Bc,
                          const SpatialDomains::MeshGraph3D &graph3D,
                          const LibUtilities::PointsType 
                          TetNb = LibUtilities::SIZE_PointsType,
                          const bool constructMap = true);

                          ContExpList3D(SpatialDomains::MeshGraph3D &graph3D,
                          const bool constructMap = true);
            
                          ContExpList3D(const ContExpList3D &In);
                        
            ~ContExpList3D();

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
                m_locToGloMap->ContToLocal(m_contCoeffs, m_coeffs);
            }
             inline void ContToLocal(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->ContToLocal(inarray, outarray);
            }
            
            inline void LocalToCont()
            {
                m_locToGloMap->LocalToCont(m_coeffs, m_contCoeffs);
            }        
        
            inline void Assemble()
            {
                m_locToGloMap->Assemble(m_coeffs, m_contCoeffs);
            }

            inline void Assemble(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->Assemble(inarray, outarray);
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
        
        typedef boost::shared_ptr<ContExpList3D>      ContExpList3DSharedPtr;
        typedef std::vector<ContExpList3DSharedPtr>   ContExpList3DVector;
        typedef std::vector<ContExpList3DSharedPtr>::iterator ContExpList3DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define

