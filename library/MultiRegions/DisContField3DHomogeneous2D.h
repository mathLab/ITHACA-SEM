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
#include <MultiRegions/ExpList3DHomogeneous2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField3DHomogeneous2D: public ExpList3DHomogeneous2D
        {
        public:
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D();

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const LibUtilities::BasisKey &HomoBasis_y,
                      const LibUtilities::BasisKey &HomoBasis_z,
                      const NekDouble lhom_y,
                      const NekDouble lhom_z,
                      const bool useFFT,
                      const bool dealiasing,
                      const Collections::ImplementationType ImpType
                      = Collections::eNoImpType);

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const LibUtilities::BasisKey &HomoBasis_y,
                      const LibUtilities::BasisKey &HomoBasis_z,
                      const NekDouble lhom_y,
                      const NekDouble lhom_z,
                      const bool useFFT,
                      const bool dealiasing,
                      const SpatialDomains::MeshGraphSharedPtr &graph1D,
                      const std::string &variable,
                      const Collections::ImplementationType ImpType
                      = Collections::eNoImpType);
            
            /// Copy constructor.
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(const DisContField3DHomogeneous2D &In,
															 const bool DeclareLinesSetCoeffPhys = true);

            /// Destructor. 
            MULTI_REGIONS_EXPORT virtual ~DisContField3DHomogeneous2D();
            
			
            MULTI_REGIONS_EXPORT void SetupBoundaryConditions(const LibUtilities::BasisKey &HomoBasis_y,
															  const LibUtilities::BasisKey &HomoBasis_z,
															  const NekDouble lhom_y,
															  const NekDouble lhom_z,
															  SpatialDomains::BoundaryConditions &bcs);
            
            MULTI_REGIONS_EXPORT void EvaluateBoundaryConditions(
                const NekDouble time = 0.0,
                const std::string varName = "");
			
			inline const Array<OneD,const MultiRegions::ExpListSharedPtr> &GetBndCondExpansions();
			
			inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &GetBndConditions();
			
			inline std::shared_ptr<ExpList> &UpdateBndCondExpansion(int i);
			
			inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>& UpdateBndConditions();
			
			/// \brief Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            MULTI_REGIONS_EXPORT void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,Array<OneD,int> &EdgeID);
            
            virtual void v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays);
			
			/// Storage space for the boundary to element and boundary to trace map.
            /// This member variable is really allocated just in case a boundary expansion
            /// recasting is required at the solver level. Otherwise is the 2 vectors are not filled up.
            /// If is needed all the funcitons whihc require to use this map do not have to recalculate it anymore.
            Array<OneD, int> m_BCtoElmMap;
            Array<OneD, int> m_BCtoEdgMap;
            
        protected:
           
            Array<OneD, MultiRegions::ExpListSharedPtr>     m_bndCondExpansions;
            
            Array<OneD, NekDouble >                         m_bndCondBndWeight;

            Array<OneD,SpatialDomains::BoundaryConditionShPtr>  m_bndConditions;
			
			virtual void v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                Array<OneD,int> &EdgeID)
            {
                GetBoundaryToElmtMap(ElmtID,EdgeID);
            }

           /// @todo Fix Robin BCs for homogeneous case
           virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo()
           {
               return std::map<int, RobinBCInfoSharedPtr>();
           }

        private:
            // virtual functions
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);
            
            virtual void v_EvaluateBoundaryConditions(
                const NekDouble   time    = 0.0,
                const std::string varName = "",
                const NekDouble   x2_in   = NekConstants::kNekUnsetDouble,
                const NekDouble   x3_in   = NekConstants::kNekUnsetDouble);
			
			virtual const Array<OneD,const std::shared_ptr<ExpList> > &v_GetBndCondExpansions(void);
			
			virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &v_GetBndConditions();
			
            virtual std::shared_ptr<ExpList> &v_UpdateBndCondExpansion(int i);
			
			virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr>& v_UpdateBndConditions();

            inline virtual void v_SetBndCondBwdWeight(
                const int index, 
                const NekDouble value);
        };

        typedef std::shared_ptr<DisContField3DHomogeneous2D>  
            DisContField3DHomogeneous2DSharedPtr;
		
		inline const Array<OneD,const MultiRegions::ExpListSharedPtr> &DisContField3DHomogeneous2D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }
		
		inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &DisContField3DHomogeneous2D::GetBndConditions()
		{
			return m_bndConditions;
		}
		
		inline MultiRegions::ExpListSharedPtr &DisContField3DHomogeneous2D::UpdateBndCondExpansion(int i)
		{
			return m_bndCondExpansions[i];
		}
		
		inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>&DisContField3DHomogeneous2D::UpdateBndConditions()
		{
			return m_bndConditions;
		}

        inline void DisContField3DHomogeneous2D::v_SetBndCondBwdWeight(
            const int index, 
            const NekDouble value)
        {
            m_bndCondBndWeight[index]   =   value;
        }
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3DHOMO2D_H
