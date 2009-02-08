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
	    
	    inline LocalToGlobalDGMapSharedPtr &GetTraceMap(void)
            {
                return m_traceMap;
            }
            
            
            void HelmSolve(const ExpList &Fce, NekDouble lambda, NekDouble tau = 10);
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
        
            /**
             * \brief This method extracts the "forward" and
             * "backward" trace data from the array #m_phys and puts
             * the data into output vectors \a Fwd and \a Bwd.
             *
             * This is a wrapper call.
             *
             * \param field is a NekDouble array which contains the 2D
             * data from which we wish to extract the backward and
             * forward orientated trace/edge arrays.
             *
             * \return Updates  a NekDouble array \a Fwd and \a Bwd
             */ 

            void GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd);

            /**
             * \brief This method extracts the "forward" and
             * "backward" trace data from the array \a field and puts
             * the data into output vectors \a Fwd and \a Bwd.
             *
             * An element unique edge is defined to be #eForwards if
             * the edge is oriented in a counter clockwise sense
             * within the element. Conversley it is defined to be
             * #eBackwards if the elemet edge is orientated in a
             * clockwise sense. Therefore along two intersecting edges
             * one edge is always forwards and the adjacent edge is
             * backwards. We define a unique normal between two
             * adjacent edges as running from the #eFowards edge to the
             * #eBackward edge. 
             *
             * This method collects/interpolates the edge data from
             * the 2D array \a field which contains information over a
             * collection of 2D shapes and puts this edge data into
             * the arrays of trace data \a Bwd or \a Fwd depending on
             * the orientation of the local edge within an element.
             *
             * If an edge is aligned along a boundary we use the
             * method GetBndExpAdjacentOrient() method to determine if
             * an adjacent boundary edge is orientated in a forwards
             * or backwards sense. This method returns an enum
             * #AdjacentTraceOrientation which in 2D has entires of
             * #eAdjacentEdgeIsForwards and #eAdjacentEdgeIsBackwards.
             *
             * \param field is a NekDouble array which contains the 2D
             * data from which we wish to extract the backward and
             * forward orientated trace/edge arrays.
             *
             * \return Updates  a NekDouble array \a Fwd and \a Bwd
             */ 
            void GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field, 
                                    Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd);

            void ExtractTracePhys()
            {
                ExtractTracePhys(m_trace->UpdatePhys());
            }


            void ExtractTracePhys(Array<OneD,NekDouble> &outarray);
            

            /**
             * \brief This method extracts the trace (edges in 2D)
             * from the field \a inarray and puts the values in \a
             * outarray. 

             * It assumes the field is C0 continuous so that
             * it can overwrite the edge data when visited by the two
             * adjacent elements.
             *
             * \param inarray is a NekDouble array which contains the 2D
             * data from which we wish to extract the edge data 
             *
             * \return Updates a NekDouble array \a outarray which
             * contains the edge information
             */ 
            void ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, 
                                  Array<OneD, NekDouble> &outarray);
            

            void AddTraceIntegral(const Array<OneD, const NekDouble> &Fx, 
                                  const Array<OneD, const NekDouble> &Fy, 
                                  Array<OneD, NekDouble> &outarray);

            
            void AddTraceIntegral(const Array<OneD, const NekDouble> &Fn, 
                                  Array<OneD, NekDouble> &outarray);

            inline const Array<OneD,const MultiRegions::ExpList1DSharedPtr>& GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }
            
            inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& GetBndConditions()
            {
                return m_bndConditions;
            }

        protected:

        private:
            Array<OneD,MultiRegions::ExpList1DSharedPtr>       m_bndCondExpansions;
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;
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

            // virtual functions
            inline virtual GenExpList1DSharedPtr &v_GetTrace(void)
            {
                return GetTrace();
            }

	    inline virtual LocalToGlobalDGMapSharedPtr &v_GetTraceMap(void)
            {
                return GetTraceMap();
            }

            virtual void v_AddTraceIntegral(const Array<OneD, const NekDouble> &Fx, 
                                          const Array<OneD, const NekDouble> &Fy, 
                                          Array<OneD, NekDouble> &outarray)
            {
                AddTraceIntegral(Fx,Fy,outarray);
            }

            virtual void v_AddTraceIntegral(const Array<OneD, const NekDouble> &Fn, 
                                          Array<OneD, NekDouble> &outarray)
            {
                AddTraceIntegral(Fn,outarray);
            }


            virtual void v_GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                              Array<OneD,NekDouble> &Bwd)
            {
                GetFwdBwdTracePhys(Fwd,Bwd);
            }

            virtual void v_GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field, 
                                              Array<OneD,NekDouble> &Fwd, 
                                              Array<OneD,NekDouble> &Bwd)
            {
                GetFwdBwdTracePhys(field, Fwd,Bwd);
            }

            virtual void v_ExtractTracePhys(Array<OneD,NekDouble> &outarray)
            {
                ExtractTracePhys(outarray);
            }
            
            virtual void v_ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
            {
                ExtractTracePhys(inarray,outarray);
            }

	    inline virtual const Array<OneD,const MultiRegions::ExpList1DSharedPtr> & v_GetBndCondExpansions(void)
            {
	      return GetBndCondExpansions();
            }

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                return GetBndConditions();
            }

            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0)
            {
                EvaluateBoundaryConditions(time);
            }
	    
	    virtual void v_HelmSolve(const ExpList &In, 
                                     NekDouble lambda,
				     Array<OneD, NekDouble>& dirForcing = NullNekDouble1DArray)
	    {
	      HelmSolve(In,lambda);
            }
        };

        typedef boost::shared_ptr<DisContField2D>   DisContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
