///////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.h
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
// Description: Field definition in two-dimensions for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/GenExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/LocalToGlobalDGMap.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField2D: public ExpList2D
        {
        public:
            DisContField2D();

            DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                SpatialDomains::BoundaryConditions &bcs, 
                           const int bc_loc = 0);

            DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                           SpatialDomains::BoundaryConditions &bcs, 
                           const std::string variable);
            
            DisContField2D(const DisContField2D &In);

            ~DisContField2D();
            
            inline GenExpList1DSharedPtr &GetTrace(void)
            {
                return m_trace;
            }

            
            void HelmSolve(DisContField2D &Fce, NekDouble lambda);
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
            void EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                ExpList2D::EvaluateBoundaryConditions(time,m_bndCondExpansions,m_bndConditions);
            }

            GlobalLinSysSharedPtr GetGlobalBndLinSys(const GlobalLinSysKey &mkey);
        
            void GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd);

            void ExtractTracePhys()
            {
                ExtractTracePhys(m_trace->UpdatePhys());
            }

            void ExtractTracePhys(Array<OneD,NekDouble> &outarray);

            void AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                                  Array<OneD, const NekDouble> &Fy, 
                                  Array<OneD, NekDouble> &outarray);

        protected:

        private:
            Array<OneD,MultiRegions::ExpList1DSharedPtr>       m_bndCondExpansions;
            Array<OneD,SpatialDomains::BoundaryConditionShPtr>  m_bndConditions;
            GlobalLinSysMapShPtr                               m_globalBndMat;
            GenExpList1DSharedPtr                              m_trace;
            LocalToGlobalDGMapSharedPtr                        m_traceMap;


            /**
             * \brief This function discretises the boundary conditions by setting up
             * a list of one-dimensional boundary expansions.    
             *
             * According to their boundary region, the separate segmental boundary 
             * expansions are bundled together in an object of the class 
             * MultiRegions#ExpList1D. 
             * The list of expansions of the Dirichlet boundary regions are listed 
             * first in the array #m_bndCondExpansions.
             *
             * \param graph2D A mesh, containing information about the domain and 
             * the spectral/hp element expansion.
             * \param bcs An entity containing information about the boundaries and 
             * boundary conditions.
             * \param variable An optional parameter to indicate for which variable 
             * the boundary conditions should be discretised.
             */ 
            void GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                                    SpatialDomains::BoundaryConditions &bcs, 
                                                    const std::string variable);
        };

        typedef boost::shared_ptr<DisContField2D>   DisContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
