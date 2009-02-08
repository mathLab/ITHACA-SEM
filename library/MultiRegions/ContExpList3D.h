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
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/GlobalLinSys.h>

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
        
            inline int GetContNcoeffs()
            {
                return m_contNcoeffs;
            }
            
            inline void GlobalToLocal()
            {
                m_locToGloMap->GlobalToLocal(m_contCoeffs, m_coeffs);
            }
            inline void GlobalToLocal(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->GlobalToLocal(inarray, outarray);
            }
            
            inline void LocalToGlobal()
            {
                m_locToGloMap->LocalToGlobal(m_coeffs, m_contCoeffs);
            }        
        
            inline void Assemble()
            {
                m_locToGloMap->Assemble(m_coeffs, m_contCoeffs);
            }

            inline void Assemble(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->Assemble(inarray, outarray);
            }

            inline const LocalToGlobalC0ContMapSharedPtr& GetLocalToGlobalMap() const
            {
                return  m_locToGloMap;
            }

            void IProductWRTBase(const ExpList &In);

            void IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
                                 Array<OneD, NekDouble> &outarray, 
                                 Array<OneD, NekDouble> &wksp = NullNekDouble1DArray);
        
            void FwdTrans(const ExpList &In);
        
            void FwdTrans(const Array<OneD, const NekDouble> &inarray,
                          Array<OneD, NekDouble> &outarray);
                             

            void BwdTrans(const ExpList &In);

            void BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray);

            void GeneralMatrixOp(const GlobalLinSysKey            &gkey,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble>          &outarray);
            

        protected:
            LocalToGlobalC0ContMapSharedPtr m_locToGloMap;
      	    int                             m_contNcoeffs;
            Array<OneD, NekDouble>          m_contCoeffs;
            GlobalLinSysMapShPtr            m_globalMat;
                        
        private:
            
            
            // virtual functions
            virtual  const Array<OneD, const NekDouble> &v_GetContCoeffs() const 
            {
                return m_contCoeffs;;
            }

            virtual void v_BwdTrans(const ExpList &Sin)
            {
                BwdTrans(Sin);
            }

            virtual void v_BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
            }

            
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }

            virtual void v_FwdTrans(const ExpList &Sin)
            {
                FwdTrans(Sin);
            }

            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual void v_IProductWRTBase(const ExpList &Sin)
            {
                IProductWRTBase(Sin);
            }

        };
        
        typedef boost::shared_ptr<ContExpList3D>      ContExpList3DSharedPtr;
        typedef std::vector<ContExpList3DSharedPtr>   ContExpList3DVector;
        typedef std::vector<ContExpList3DSharedPtr>::iterator ContExpList3DVectorIter;

    } //end of namespace
} //end of namespace

#endif // end of define

