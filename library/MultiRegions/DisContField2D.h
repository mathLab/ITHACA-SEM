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
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/LocalToGlobalDGMap.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/Conditions.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField2D: public ExpList2D
        {
        public:
            MULTI_REGIONS_EXPORT DisContField2D();

            MULTI_REGIONS_EXPORT DisContField2D(
                           LibUtilities::SessionReaderSharedPtr &pSession,
                           SpatialDomains::MeshGraphSharedPtr &graph2D,
                           const std::string variable,
                           bool SetUpJustDG = true,
                           bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT DisContField2D(const DisContField2D &In,
                           SpatialDomains::MeshGraphSharedPtr &graph2D,
                           const std::string variable,
                           bool SetUpJustDG = false,
                           bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT DisContField2D(const DisContField2D &In,
                           bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT ~DisContField2D();

            inline ExpList1DSharedPtr &GetTrace(void)
            {
                return m_trace;
            }

            inline LocalToGlobalDGMapSharedPtr &GetTraceMap(void)
            {
                return m_traceMap;
            }

            /// Determines if another ContField2D shares the same boundary
            /// conditions as this field.
            MULTI_REGIONS_EXPORT bool SameTypeOfBoundaryConditions(const DisContField2D &In);

            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(const GlobalLinSysKey &mkey);

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

            MULTI_REGIONS_EXPORT void GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd,
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
            MULTI_REGIONS_EXPORT void GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field,
                                    Array<OneD,NekDouble> &Fwd,
                                    Array<OneD,NekDouble> &Bwd);

            void ExtractTracePhys()
            {
                ExtractTracePhys(m_trace->UpdatePhys());
            }


            MULTI_REGIONS_EXPORT void ExtractTracePhys(Array<OneD,NekDouble> &outarray);


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
            MULTI_REGIONS_EXPORT void ExtractTracePhys(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray);


            MULTI_REGIONS_EXPORT void AddTraceIntegral(const Array<OneD, const NekDouble> &Fx,
                                  const Array<OneD, const NekDouble> &Fy,
                                  Array<OneD, NekDouble> &outarray);


            MULTI_REGIONS_EXPORT void AddTraceIntegral(const Array<OneD, const NekDouble> &Fn,
                                  Array<OneD, NekDouble> &outarray);

            MULTI_REGIONS_EXPORT void AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd,
                                    const Array<OneD, const NekDouble> &Bwd,
                                    Array<OneD, NekDouble> &outarray);

            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>& GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

	    
            inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& GetBndConditions()
            {
                return m_bndConditions;
            }

            /// \brief Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            MULTI_REGIONS_EXPORT void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                      Array<OneD,int> &EdgeID);

            MULTI_REGIONS_EXPORT NekDouble L2_DGDeriv(const int dir,
                                 const Array<OneD, const NekDouble> &soln);

            /// \brief Set up an stl map containing the information
            /// for a robin aboundary condition in the location of the
            /// element id
            MULTI_REGIONS_EXPORT map<int, RobinBCInfoSharedPtr> GetRobinBCInfo(void);

            MULTI_REGIONS_EXPORT void EvaluateHDGPostProcessing(Array<OneD, NekDouble> &outarray);

            /// Generates a map of periodic edges in the mesh.
            void GetPeriodicEdges(
                        SpatialDomains::MeshGraphSharedPtr &graph2D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        vector<map<int,int> > & periodicVertices,
                        map<int,int>& periodicEdges);

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
            Array<OneD,MultiRegions::ExpListSharedPtr>       m_bndCondExpansions;              
            

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
            void GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraphSharedPtr &graph2D,
                                                    SpatialDomains::BoundaryConditions &bcs,
                                                    const std::string variable,
                                                    bool DeclareCoeffPhysArrays = true);


            virtual void v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                Array<OneD,int> &EdgeID)
            {
                GetBoundaryToElmtMap(ElmtID,EdgeID);
            }

        private:
            GlobalLinSysMapShPtr                               m_globalBndMat;
            ExpList1DSharedPtr                                 m_trace;
            LocalToGlobalDGMapSharedPtr                        m_traceMap;


            // virtual functions
            inline virtual ExpList1DSharedPtr &v_GetTrace(void)
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

            virtual void v_AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd,
                                          const Array<OneD, const NekDouble> &Bwd,
                                          Array<OneD, NekDouble> &outarray)
            {
                AddTraceBiIntegral(Fwd,Bwd,outarray);
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

	    inline virtual const Array<OneD,const MultiRegions::ExpListSharedPtr> & v_GetBndCondExpansions(void)
            {
	      return GetBndCondExpansions();
            }
	    
            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                return GetBndConditions();
            }
	    
            inline MultiRegions::ExpListSharedPtr &v_UpdateBndCondExpansion(int i)
            {
                return m_bndCondExpansions[i];
            }

            inline Array<OneD, SpatialDomains::BoundaryConditionShPtr> &v_UpdateBndConditions()
            {
                return m_bndConditions;
            }


            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0,
													  const NekDouble x2_in = NekConstants::kNekUnsetDouble,
													  const NekDouble x3_in = NekConstants::kNekUnsetDouble);

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo()
            {
                return GetRobinBCInfo();
            }

            virtual void v_GetPeriodicEdges(SpatialDomains::MeshGraphSharedPtr &graph2D,
                                            SpatialDomains::BoundaryConditions &bcs,
                                            const std::string variable,
                                            vector<map<int,int> > & periodicVertices,
                                            map<int,int>& periodicEdges)
            {
                GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);
            }

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

            /// Calculates the result of the multiplication of a global
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            virtual void v_GeneralMatrixOp(
                   const GlobalMatrixKey             &gkey,
                   const Array<OneD,const NekDouble> &inarray,
                   Array<OneD,      NekDouble> &outarray,
                   bool UseContCoeffs = false);

        };

        typedef boost::shared_ptr<DisContField2D>   DisContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
