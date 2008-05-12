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
#include <MultiRegions/ContExpList3D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>

#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/BoundaryConditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class ContField3D: public ContExpList3D
        {
        public:
            ContField3D();
            ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs, 
                        const int bc_loc = 0);
                        
            ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable);

            ContField3D(const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc,
                        SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const int bc_loc = 0,
                        const LibUtilities::PointsType 
                        TriNb = LibUtilities::SIZE_PointsType);

            ContField3D(const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc,
                        SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        const LibUtilities::PointsType 
                        TriNb = LibUtilities::SIZE_PointsType);

            ContField3D(const ContField3D &In);

            ~ContField3D();

            void FwdTrans (const ExpList &In);
            void HelmSolve(const ExpList &In, NekDouble lambda);

        protected:

        private:
            Array<OneD,MultiRegions::ExpList2DSharedPtr>           m_bndConstraint;
            Array<OneD,SpatialDomains::BoundaryConditionType>      m_bndTypes;
            
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            void GlobalSolve(const GlobalLinSysKey &key, const ExpList &Rhs, NekDouble ScaleForcing=1.0);

            void GenerateField3D(SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable);

        };
        typedef boost::shared_ptr<ContField3D>      ContField3DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3D_H
