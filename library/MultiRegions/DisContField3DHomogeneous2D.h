///////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous2D.h
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
// Description: Field definition in three-dimensions for a discontinuous
// LDG-H expansion with 2 homogeneous directions in 2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO2D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/ExpList1DHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/DisContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField3DHomogeneous2D: public ExpList3DHomogeneous2D
        {
        public:
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D();

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(const LibUtilities::BasisKey &HomoBasis_y,
															 const LibUtilities::BasisKey &HomoBasis_z,
															 const NekDouble lhom_y,
															 const NekDouble lhom_z,
															 bool useFFT);

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(const LibUtilities::BasisKey &HomoBasis_y,
															 const LibUtilities::BasisKey &HomoBasis_z,
															 const NekDouble lhom_y,
															 const NekDouble lhom_z,
															 bool useFFT,
															 SpatialDomains::MeshGraph1D &graph1D,
															 SpatialDomains::BoundaryConditions &bcs, 
															 const int bc_loc = 0,
															 const GlobalSysSolnType solnType = eDirectMultiLevelStaticCond);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(const DisContField3DHomogeneous2D &In,
                                        bool DeclareLinesSetCoeffPhys = true);

            /// Destructor. 
            MULTI_REGIONS_EXPORT ~DisContField3DHomogeneous2D();
            
			
            MULTI_REGIONS_EXPORT void SetupBoundaryConditions(const LibUtilities::BasisKey &HomoBasis_y,
															  const LibUtilities::BasisKey &HomoBasis_z,
															  const NekDouble lhom_y,
															  const NekDouble lhom_z,
															  SpatialDomains::BoundaryConditions &bcs);
            
            MULTI_REGIONS_EXPORT void EvaluateBoundaryConditions(const NekDouble time = 0.0);
            
        protected:
           
            Array<OneD,MultiRegions::ExpList1DHomogeneous2DSharedPtr>   m_bndCondExpansions;

            Array<OneD,SpatialDomains::BoundaryConditionShPtr>          m_bndConditions;

        private:
            // virtual functions
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff);
            
            virtual void v_HelmSolveDG(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                    NekDouble tau);
			
			virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0,
													  const NekDouble x2_in = NekConstants::kNekUnsetDouble);
        };

        typedef boost::shared_ptr<DisContField3DHomogeneous2D>  
            DisContField3DHomogeneous2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3DHOMO2D_H
