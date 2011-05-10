///////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous1D.h
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
// LDG-H expansion with a homogeneous direction in 1D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO1D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO1D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/DisContField2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField3DHomogeneous1D: public ExpList3DHomogeneous1D
        {
        public:
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D();

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D(const LibUtilities::BasisKey &HomoBasis,
                                        const NekDouble lhom, bool useFFT);

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D(
                           const LibUtilities::BasisKey &HomoBasis,
                           const NekDouble lhom,
						   bool useFFT,
                           SpatialDomains::MeshGraph2D &graph2D,
                           SpatialDomains::BoundaryConditions &bcs, 
                           const int bc_loc = 0,
                           const GlobalSysSolnType solnType = eDirectMultiLevelStaticCond);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D(const DisContField3DHomogeneous1D &In,
                                        bool DeclarePlanesSetCoeffPhys = true);

            /// Destructor. 
            MULTI_REGIONS_EXPORT ~DisContField3DHomogeneous1D();
            
            MULTI_REGIONS_EXPORT void SetupBoundaryConditions(const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom, SpatialDomains::BoundaryConditions &bcs);
            /**
             * \brief This function evaluates the boundary conditions at a certain 
             * time-level.
             *
             * Based on the boundary condition \f$g(\boldsymbol{x},t)\f$ evaluated
             * at a given time-level \a t, this function transforms the boundary 
             * conditions onto the coefficients of the (one-dimensional) boundary 
             * expansion. Depending on the type of boundary conditions, these
             * expansion coefficients are calculated in different ways:
             * - <b>Dirichlet boundary conditions</b><BR>
             *   In order to ensure global \f$C^0\f$ continuity of the spectral/hp 
             *   approximation, the Dirichlet boundary conditions are projected onto 
             *   the boundary expansion by means of a modified \f$C^0\f$ continuous  
             *   Galerkin projection. This projection can be viewed as a collocation
             *   projection at the vertices, followed by an \f$L^2\f$ projection on 
             *   the interior modes of the edges. The resulting coefficients 
             *   \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ will be stored for the 
             *   boundary expansion.
             * - <b>Neumann boundary conditions</b>
             *   In the discrete Galerkin formulation of the problem to be solved, 
             *   the Neumann boundary conditions appear as the set of surface 
             *   integrals: \f[\boldsymbol{\hat{g}}=\int_{\Gamma}
             *   \phi^e_n(\boldsymbol{x})g(\boldsymbol{x})d(\boldsymbol{x})\quad
             *   \forall n \f]
             *   As a result, it are the coefficients \f$\boldsymbol{\hat{g}}\f$ 
             *   that will be stored in the boundary expansion
             *
             * \param time The time at which the boundary conditions should be 
             * evaluated
             */ 
            MULTI_REGIONS_EXPORT void EvaluateBoundaryConditions(const NekDouble time = 0.0);
            
        protected:
            /**
             * \brief An object which contains the discretised
             * boundary conditions.
             *
             * It is an array of size equal to the number of boundary
             * regions and consists of entries of the type
             * MultiRegions#ExpList1D. Every entry corresponds to the
             * one-dimensional spectral/hp expansion on a single
             * boundary region.  The values of the boundary conditions
             * are stored as the coefficients of the one-dimensional
             * expansion.
             */ 
            Array<OneD,MultiRegions::ExpList2DHomogeneous1DSharedPtr>   m_bndCondExpansions;

            /**
             * \brief An array which contains the information about
             * the boundary condition on the different boundary
             * regions.
             */ 
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

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

        typedef boost::shared_ptr<DisContField3DHomogeneous1D>  
            DisContField3DHomogeneous1DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3DHOMO1D_H
