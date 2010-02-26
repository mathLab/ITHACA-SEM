///////////////////////////////////////////////////////////////////////////////
//
// File DisContField3D.h
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
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/LocalToGlobalDGMap.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/BoundaryConditions.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Abstraction of a global discontinuous three-dimensional spectral/hp
        /// element expansion which approximates the solution of a set of
        /// partial differential equations.
        class DisContField3D: public ExpList3D
        {
        public:
            /// Default constructor
            DisContField3D();
            
            /// Constructs a global discontinuous field based on an input mesh
            /// with boundary conditions.
            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
                           SpatialDomains::BoundaryConditions &bcs,
                           const int bc_loc,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);
                           
            /// Constructs a global discontinuous field based on an input mesh
            /// with boundary conditions.
            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
                           SpatialDomains::BoundaryConditions &bcs,
                           const std::string variable,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);
                           
            /// Constructs a global discontinuous field based on an input mesh.
            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
                          SpatialDomains::BoundaryConditions &bcs,
                          const GlobalSysSolnType solnType = eDirectStaticCond,
                          const bool constructMap = true);
                          
            /// Constructs a global discontinuous field based on another
            /// discontinuous field.
            DisContField3D(const DisContField3D &In);
            
            /// Destructor.
            ~DisContField3D();

            void EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                ExpList3D::EvaluateBoundaryConditions(time,m_bndCondExpansions,m_bndConditions);
            }

        protected:
            /**
             * \brief The number of boundary segments on which
             * Dirichlet boundary conditions are imposed
             */ 
            int m_numDirBndCondExpansions;

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
            Array<OneD,MultiRegions::ExpList2DSharedPtr>       m_bndCondExpansions;

            /**
             * \brief An array which contains the information about
             * the boundary condition on the different boundary
             * regions.
             */ 
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

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
            void GenerateBoundaryConditionExpansion(
                    SpatialDomains::MeshGraph3D &graph3D, 
                    SpatialDomains::BoundaryConditions &bcs, 
                    const std::string variable);
        
        private:
            GlobalLinSysMapShPtr                                m_globalBndMat;
            ExpList1DSharedPtr                                  m_trace;
            LocalToGlobalDGMapSharedPtr                         m_traceMap;        
        };

        typedef boost::shared_ptr<DisContField3D>   DisContField3DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3D_H
