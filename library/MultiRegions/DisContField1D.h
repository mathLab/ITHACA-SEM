///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.h
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
// Description: Field definition in one-dimension for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList0D.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/LocalToGlobalDGMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// This class is the abstraction of a global discontinuous two-
        /// dimensional spectral/hp element expansion which approximates the
        /// solution of a set of partial differential equations.
        class DisContField1D: public ExpList1D
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT DisContField1D();

            /// Constructs a 1D discontinuous field based on a mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT DisContField1D(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const SpatialDomains::MeshGraphSharedPtr &graph1D,
                    const std::string &variable);

            /// Constructs a 1D discontinuous field based on an existing field.
            MULTI_REGIONS_EXPORT DisContField1D(const DisContField1D &In);

            /// Constructs a 1D discontinuous field based on an existing field.
	    /// (needed in order to use ContField( const ExpList1D &In) constructor
            MULTI_REGIONS_EXPORT DisContField1D(const ExpList1D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~DisContField1D();
			
			////2D
			inline ExpList0DSharedPtr &GetTrace1D(void)
            {
                return m_trace;
            }
			
			////2D
			inline LocalToGlobalDGMapSharedPtr &GetTraceMap(void)
            {
                return m_traceMap;
            }
			
            /// For a given key, returns the associated global linear system.
            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(
                    const GlobalLinSysKey &mkey);

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
			
			
            MULTI_REGIONS_EXPORT void AddTraceIntegral(const Array<OneD, const NekDouble> &Fn,
													   Array<OneD, NekDouble> &outarray);
			
            /// Retrieve the boundary condition expansions.
            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>&GetBndCondExpansions();
			
			inline MultiRegions::ExpListSharedPtr &UpdateBndCondExpansion(int i);
			
			inline Array<OneD, SpatialDomains::BoundaryConditionShPtr> &UpdateBndConditions();

            MULTI_REGIONS_EXPORT void GetBoundaryToElmtMap(Array<OneD,int> &ElmtID, Array<OneD,int> &VertID);

            /// \brief Set up an stl map containing the information
            /// for a robin aboundary condition in the location of the
            /// element id
            MULTI_REGIONS_EXPORT map<int, RobinBCInfoSharedPtr> GetRobinBCInfo(void);
			
        protected:
            /// The number of boundary segments on which Dirichlet boundary
            /// conditions are imposed.
            int m_numDirBndCondExpansions;

            /// Discretised boundary conditions.
            /**
             * It is an array of size equal to the number of boundary points
             * and consists of entries of the type LocalRegions#PointExp. Every
             * entry corresponds to a point on a single boundary region.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr> m_bndCondExpansions;

            /// An array which contains the information about the boundary
            /// condition on the different boundary regions.
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            /// Discretises the boundary conditions.
            void GenerateBoundaryConditionExpansion(
                    const SpatialDomains::MeshGraphSharedPtr &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable);

            /// Generate a associative map of periodic vertices in a mesh.
            void GetPeriodicVertices(
                                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                const SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                      map<int,int>& periodicVertices);

            virtual void v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                Array<OneD,int> &EdgeID)
            {
                GetBoundaryToElmtMap(ElmtID,EdgeID);
            }

        private:

            /// Global boundary matrix.
            GlobalLinSysMapShPtr                               m_globalBndMat;

            /// Trace space storage for points between elements.
			ExpList0DSharedPtr                                 m_trace;
			Array<OneD,NekDouble>                              tmpBndSol;


            /// Local to global DG mapping for trace space.
            LocalToGlobalDGMapSharedPtr                        m_traceMap;
			
			
			////2D
			inline virtual ExpList0DSharedPtr &v_GetTrace1D(void)
            {
                return GetTrace1D();
            }
			
			inline virtual LocalToGlobalDGMapSharedPtr &v_GetTraceMap(void)
			{
				return GetTraceMap();
			}
			 
			virtual void v_AddTraceIntegral(const Array<OneD, const NekDouble> &Fn,
											Array<OneD, NekDouble> &outarray)
			{
				AddTraceIntegral(Fn,outarray);
			}
			 
			/*virtual void v_AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd,
											  const Array<OneD, const NekDouble> &Bwd,
											  Array<OneD, NekDouble> &outarray)
			{
				AddTraceBiIntegral(Fwd,Bwd,outarray);
			}*/
			 
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
			 

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                const SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                Array<OneD, MultiRegions::ExpListSharedPtr>
                                                            &bndCondExpansions,
                                Array<OneD, SpatialDomains
                                    ::BoundaryConditionShPtr> &bndConditions);

            void GenerateFieldBnd1D(
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable);

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo()
            {
                return GetRobinBCInfo();
            }
			
			
			/// Retrieve the boundary condition descriptions.
            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions();
			
			inline Array<OneD, SpatialDomains::BoundaryConditionShPtr> &v_UpdateBndConditions()
            {
                return m_bndConditions;
            }
			
			inline MultiRegions::ExpListSharedPtr &v_UpdateBndCondExpansion(int i)
            {
                return m_bndCondExpansions[i];
            }
			
            /// Evaluate all boundary conditions at a given time..
            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0,
													  const NekDouble x2_in = NekConstants::kNekUnsetDouble,
													  const NekDouble x3_in = NekConstants::kNekUnsetDouble);

            /// Solve the Helmholtz equation.
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const Array<OneD, const NekDouble> &dirForcing);
        };

        typedef boost::shared_ptr<DisContField1D>   DisContField1DSharedPtr;

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr> &DisContField1D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }
		
		inline Array<OneD, SpatialDomains::BoundaryConditionShPtr> &DisContField1D::UpdateBndConditions()
		{
			return m_bndConditions;
		}
		
		inline MultiRegions::ExpListSharedPtr &DisContField1D::UpdateBndCondExpansion(int i)
		{
			return m_bndCondExpansions[i];
		}
		
		

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD1D_H
