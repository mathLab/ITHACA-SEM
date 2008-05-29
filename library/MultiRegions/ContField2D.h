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
// Description: Field definition in tow-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>

//#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>


namespace Nektar
{
    namespace MultiRegions
    {
        class ContField2D: public ContExpList2D
        {
        public:
            ContField2D();

            ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                SpatialDomains::BoundaryConditions &bcs, 
                const int bc_loc = 0);

            ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                SpatialDomains::BoundaryConditions &bcs, 
                const std::string variable);

            ContField2D(const LibUtilities::BasisKey &TriBa, 
                        const LibUtilities::BasisKey &TriBb, 
                        const LibUtilities::BasisKey &QuadBa, 
                        const LibUtilities::BasisKey &QuadBb, 
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const int bc_loc = 0,
                        const LibUtilities::PointsType 
                        TriNb = LibUtilities::SIZE_PointsType);
            
            ContField2D(const LibUtilities::BasisKey &TriBa, 
                        const LibUtilities::BasisKey &TriBb, 
                        const LibUtilities::BasisKey &QuadBa, 
                        const LibUtilities::BasisKey &QuadBb, 
                        SpatialDomains::MeshGraph2D &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        const LibUtilities::PointsType 
                        TriNb = LibUtilities::SIZE_PointsType);

            ContField2D(const ContField2D &In);

            ~ContField2D();

/*             void SetBoundaryCondition(const int loc, const NekDouble value) */
/*             { */
/*                 m_bndConstraint[loc]->SetValue(value); */
/*             } */

            void FwdTrans (const ExpList &In);
            void HelmSolve(const ExpList &In, NekDouble lambda);

            void EvaluateBoundaryConditions(const NekDouble time = 0.0);

            inline const Array<OneD,const MultiRegions::ExpList1DSharedPtr>& 
                GetBndCondExp()
            {
                return m_bndConstraint;
            }

        protected:

        private:
            Array<OneD,MultiRegions::ExpList1DSharedPtr>       m_bndConstraint;
            Array<OneD,SpatialDomains::BoundaryConditionType>  m_bndTypes;
            Array<OneD,SpatialDomains::Equation>               m_bndCondEquations;

            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            void GlobalSolve(const GlobalLinSysKey &key, const ExpList &Rhs, 
                             NekDouble ScaleForcing=1.0);

            void GenerateField2D(SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable);

            void GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                                    SpatialDomains::BoundaryConditions &bcs, 
                                                    const std::string variable);
        };
        typedef boost::shared_ptr<ContField2D>      ContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD2D_H
