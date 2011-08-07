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

            /// Constructs a 1D discontinuous field based on a mesh.
            MULTI_REGIONS_EXPORT DisContField1D(
                    LibUtilities::SessionReaderSharedPtr& pSession,
                    SpatialDomains::MeshGraph1D &graph1D,
                    const bool constructMap = true);

            /// Constructs a 1D discontinuous field based on a mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT DisContField1D(
                    LibUtilities::SessionReaderSharedPtr& pSession,
                    SpatialDomains::MeshGraph1D &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const int bc_loc = 0);

            /// Constructs a 1D discontinuous field based on a mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT DisContField1D(
                    LibUtilities::SessionReaderSharedPtr& pSession,
                    SpatialDomains::MeshGraph1D &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable);

            /// Constructs a 1D discontinuous field based on a mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT DisContField1D(
                    LibUtilities::SessionReaderSharedPtr& pSession,
                    SpatialDomains::MeshGraph1D &graph1D,
                    const std::string variable);

            /// Constructs a 1D discontinuous field based on an existing field.
            MULTI_REGIONS_EXPORT DisContField1D(const DisContField1D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT ~DisContField1D();

            /// For a given key, returns the associated global linear system.
            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(
                    const GlobalLinSysKey &mkey);

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
                    const SpatialDomains::MeshGraph1D &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable);

            /// Generate a associative map of periodic vertices in a mesh.
            void GetPeriodicVertices(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
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
            Array<OneD,NekDouble>                              m_trace;

            /// Local to global DG mapping for trace space.
            LocalToGlobalDGMapSharedPtr                        m_traceMap;

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
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

            /// Solve the Helmholtz equation.
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff);

            /// Solve the Helmholtz equation (DG specific).
            virtual void v_HelmSolveDG(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                          NekDouble tau);

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
