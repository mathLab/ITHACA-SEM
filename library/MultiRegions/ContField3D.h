///////////////////////////////////////////////////////////////////////////////
//
// File ContField2D.h
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
// Description: Field definition in three-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField3D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/GlobalLinSys.h>

#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/BoundaryConditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class ContField3D: public DisContField3D
        {
        public:
            ContField3D();

            ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs, 
                        const int bc_loc = 0,
                        const GlobalSysSolnType solnType = eDirectStaticCond);

            ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs, 
                        const std::string variable,
                        const GlobalSysSolnType solnType = eDirectStaticCond);

            ContField3D(const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc,
                        SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const int bc_loc = 0,
                        const LibUtilities::PointsType 
                        TetNb = LibUtilities::SIZE_PointsType,
                        const GlobalSysSolnType solnType = eDirectStaticCond);

            ContField3D(const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc,
                        SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        const LibUtilities::PointsType 
                        TetNb = LibUtilities::SIZE_PointsType,
                        const GlobalSysSolnType solnType = eDirectStaticCond);

            ContField3D(const ContField3D &In);

            ~ContField3D();

            void FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,      NekDouble> &outarray,
                          bool  UseContCoeffs = false);
            
            void MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, 
                                               Array<OneD,       NekDouble> &outarray,
                                         bool  UseContCoeffs = false);

            void EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                ExpList3D::EvaluateBoundaryConditions(time,m_bndCondExpansions,m_bndConditions);
            }
            /**
             * \brief This function return the boundary conditions expansion.
             */ 
            inline const Array<OneD,const MultiRegions::ExpList2DSharedPtr>&GetBndCondExp()
            {
                return m_bndCondExpansions;
            }
            
            void GenerateDirBndCondForcing(const GlobalLinSysKey &key, 
                                                        Array<OneD, NekDouble> &inout, 
                                                        Array<OneD, NekDouble> &outarray);
                       inline void GlobalToLocal()
            {
                m_locToGloMap->GlobalToLocal(m_contCoeffs, m_coeffs);
            }
            inline void GlobalToLocal(const Array<OneD, const NekDouble>        &inarray, Array<OneD,NekDouble> &outarray)
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

            inline void Assemble(const Array<OneD, const NekDouble> &inarray,   Array<OneD,NekDouble> &outarray)
            {
                m_locToGloMap->Assemble(inarray, outarray);
            }

            inline const LocalToGlobalC0ContMapSharedPtr& GetLocalToGlobalMap() const
            {
                return  m_locToGloMap;
            } 

        protected:
            LocalToGlobalC0ContMapSharedPtr m_locToGloMap;
            int                             m_contNcoeffs;
            Array<OneD, NekDouble>          m_contCoeffs;
            GlobalLinSysMapShPtr            m_globalLinSys;

        private:
            /**
             * \brief The number of boundary segments on which
             * Dirichlet boundary conditions are imposed
             */ 
            int m_numDirBndCondExpansions;

            Array<OneD,MultiRegions::ExpList2DSharedPtr>           m_bndCondExpansions;
            Array<OneD,SpatialDomains::BoundaryConditionShPtr>     m_bndConditions;

            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            void GlobalSolve(const GlobalLinSysKey &key, 
                             const Array<OneD, const  NekDouble> &rhs, 
                             Array<OneD, NekDouble> &inout,
                             const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);

            void GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph3D &graph3D,
                                                    SpatialDomains::BoundaryConditions &bcs,
                                                    const std::string variable);
          
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,       NekDouble> &outarray,
                                    bool  UseContCoeffs)
            {
                FwdTrans(inarray,outarray,UseContCoeffs);
            }

            virtual void v_MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD,       NekDouble> &outarray,
                                                   bool  UseContCoeffs)
            {
                MultiplyByInvMassMatrix(inarray,outarray,UseContCoeffs);
            }

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff);
            
            virtual void v_HelmSolveCG(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                          bool UseContCoeffs,
                    const Array<OneD, const NekDouble> &dirForcing);

        };
        typedef boost::shared_ptr<ContField3D>      ContField3DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3D_H
