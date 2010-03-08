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
#include <MultiRegions/ExpList1D.h>
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
            /// Default constructor
            DisContField2D();

            /// Construct a Discontinuous Galerkin field with boundary
            /// conditions.
            DisContField2D( SpatialDomains::MeshGraph2D &graph2D,
                            SpatialDomains::BoundaryConditions &bcs, 
                            const int bc_loc = 0,
                            const GlobalSysSolnType solnType 
                                                = eDirectMultiLevelStaticCond,
                            bool SetUpJustDG = true);

            /// Construct a Discontinuous Galerkin field with boundary
            /// conditions.
            DisContField2D( SpatialDomains::MeshGraph2D &graph2D,
                            SpatialDomains::BoundaryConditions &bcs, 
                            const std::string variable,
                            const GlobalSysSolnType solnType 
                                                = eDirectMultiLevelStaticCond,
                            bool SetUpJustDG = true);

            /// Copy constructor
            DisContField2D(const DisContField2D &In);

            /// Destructor
            ~DisContField2D();
            
            /// Retrieve the 1D trace space expansion list.
            inline ExpList1DSharedPtr &GetTrace();

            /// Retrieve the local to global mapping of the trace space.
	        inline LocalToGlobalDGMapSharedPtr &GetTraceMap();

            /// This function evaluates the boundary conditions at a certain 
            /// time-level.
            inline void EvaluateBoundaryConditions(const NekDouble time = 0.0);

            /// Retrieves the fields global boundary linear system.
            GlobalLinSysSharedPtr GetGlobalBndLinSys(
                                const GlobalLinSysKey &mkey);

            /// This method extracts the "forward" and "backward" trace data 
            /// from the array #m_phys and puts the data into output vectors 
            /// \a Fwd and \a Bwd.
            void GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd);

            /// This method extracts the "forward" and "backward" trace data 
            /// from the array \a field and puts the data into output vectors 
            /// \a Fwd and \a Bwd.
            void GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field, 
                                    Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd);

            /// Wrapper call to extract the trace space.
            inline void ExtractTracePhys();

            /// Extract trace space physical space values of edges in the 2D
            /// mesh and stores them in outarray.
            inline void ExtractTracePhys(Array<OneD,NekDouble> &outarray);
            
            /// This method extracts the trace (edges in 2D) from the field 
            /// \a inarray and puts the values in \a outarray. 
            void ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, 
                                  Array<OneD, NekDouble> &outarray);
            

            void AddTraceIntegral(const Array<OneD, const NekDouble> &Fx, 
                                  const Array<OneD, const NekDouble> &Fy, 
                                  Array<OneD, NekDouble> &outarray);

            
            void AddTraceIntegral(const Array<OneD, const NekDouble> &Fn, 
                                  Array<OneD, NekDouble> &outarray);

            void AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd, 
                                    const Array<OneD, const NekDouble> &Bwd, 
                                    Array<OneD, NekDouble> &outarray);

            /// Retrieves the array of boundary condition expansion lists for
            /// each of the boundary regions.
            inline const Array<OneD, const MultiRegions::ExpList1DSharedPtr>& 
                                                        GetBndCondExpansions();
            
            /// Retrieves the array of boundary conditions for each of the
            /// boundary regions.
            inline const Array<OneD, const SpatialDomains
                                ::BoundaryConditionShPtr>& GetBndConditions();
            
            /// Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, 
                                      Array<OneD,int> &EdgeID);

            /// Calculate the L2 error of the Q_dir derivative using the 
            /// consistent DG evaluation of Q_dir.
            NekDouble L2_DGDeriv(const int dir, 
                                 const Array<OneD, const NekDouble> &soln);

        protected:
            /// The number of boundary segments on which Dirichlet boundary 
            /// conditions are imposed
            int m_numDirBndCondExpansions;

            /// An object which contains the discretised boundary conditions.
            /**
             * It is an array of size equal to the number of boundary
             * regions and consists of entries of the type
             * MultiRegions#ExpList1D. Every entry corresponds to the
             * one-dimensional spectral/hp expansion on a single
             * boundary region.  The values of the boundary conditions
             * are stored as the coefficients of the one-dimensional
             * expansion.
             */ 
            Array<OneD,MultiRegions::ExpList1DSharedPtr> m_bndCondExpansions;

            /// An array which contains the information about the boundary 
            /// condition on the different boundary regions.
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            /// This function discretises the boundary conditions by setting up
            /// a list of one-dimensional boundary expansions.    
            void GenerateBoundaryConditionExpansion(
                            SpatialDomains::MeshGraph2D &graph2D, 
                            SpatialDomains::BoundaryConditions &bcs, 
                            const std::string variable);

            /// Constructor for use by ContField2D when no boundary conditions
            /// are used.
            DisContField2D( SpatialDomains::MeshGraph2D &graph2D,
                            const GlobalSysSolnType solnType 
                                                = eDirectMultiLevelStaticCond,
                            bool SetUpJustDG = true);
            
        private:
            /// Global boundary matrix.
            GlobalLinSysMapShPtr            m_globalBndMat;
            
            /// Trace space storage for points between elements.
            ExpList1DSharedPtr              m_trace;
            
            /// Local to global DG mapping for trace space.
            LocalToGlobalDGMapSharedPtr     m_traceMap;

            virtual ExpList1DSharedPtr &v_GetTrace();

            virtual LocalToGlobalDGMapSharedPtr &v_GetTraceMap();

            virtual void v_AddTraceIntegral(
                            const Array<OneD, const NekDouble> &Fx, 
                            const Array<OneD, const NekDouble> &Fy, 
                            Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceIntegral(
                            const Array<OneD, const NekDouble> &Fn, 
                            Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceBiIntegral(
                            const Array<OneD, const NekDouble> &Fwd, 
                            const Array<OneD, const NekDouble> &Bwd, 
                            Array<OneD, NekDouble> &outarray);

            virtual void v_GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                            Array<OneD,NekDouble> &Bwd);

            virtual void v_GetFwdBwdTracePhys(
                            const Array<OneD,const NekDouble>  &field, 
                            Array<OneD,NekDouble> &Fwd, 
                            Array<OneD,NekDouble> &Bwd);

            virtual void v_ExtractTracePhys(Array<OneD,NekDouble> &outarray);
            
            virtual void v_ExtractTracePhys(
                            const Array<OneD, const NekDouble> &inarray, 
                            Array<OneD, NekDouble> &outarray);

            virtual const Array<OneD,const MultiRegions::ExpList1DSharedPtr> & 
                                                    v_GetBndCondExpansions();

            virtual const Array<OneD,const SpatialDomains
                                ::BoundaryConditionShPtr>& v_GetBndConditions();

            virtual void v_GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, 
                            Array<OneD,int> &EdgeID);

            virtual void v_EvaluateBoundaryConditions(
                            const NekDouble time = 0.0);

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
        };

        typedef boost::shared_ptr<DisContField2D>   DisContField2DSharedPtr;

        inline ExpList1DSharedPtr& DisContField2D::GetTrace()
        {
            return m_trace;
        }
    
        inline LocalToGlobalDGMapSharedPtr& DisContField2D::GetTraceMap()
        {
            return m_traceMap;
        }

        /**
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
         * @param   time        The time at which the boundary conditions 
         *                      should be evaluated.
         */ 
        inline void DisContField2D::EvaluateBoundaryConditions(
                            const NekDouble time)
        {
            ExpList2D::EvaluateBoundaryConditions(time,m_bndCondExpansions,
                                                       m_bndConditions);
        }

        inline void DisContField2D::ExtractTracePhys()
        {
            ExtractTracePhys(m_trace->UpdatePhys());
        }

        inline void DisContField2D::ExtractTracePhys(
                            Array<OneD,NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            ExtractTracePhys(m_phys, outarray);
        }

        inline const Array<OneD,const MultiRegions::ExpList1DSharedPtr>& 
                                        DisContField2D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }
        
        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& 
                                        DisContField2D::GetBndConditions()
        {
            return m_bndConditions;
        }

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
