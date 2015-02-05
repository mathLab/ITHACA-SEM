///////////////////////////////////////////////////////////////////////////////
//
// File: Mapping.h
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
// Description: Abstract base class for mappings.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_MAPPING
#define NEKTAR_SOLVERUTILS_MAPPING

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace SolverUtils
{
    //  Forward declaration
    class Mapping;

    /// A shared pointer to a Mapping object
    SOLVER_UTILS_EXPORT typedef boost::shared_ptr<Mapping> MappingSharedPtr;

    /// Declaration of the mapping factory
    typedef LibUtilities::NekFactory<std::string, Mapping,
            const LibUtilities::SessionReaderSharedPtr&,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&,
            const TiXmlElement*> MappingFactory;

    /// Declaration of the mapping factory singleton
    SOLVER_UTILS_EXPORT MappingFactory& GetMappingFactory();

    /**
     * @class Mapping
     * @brief Defines a mapping to be applied to the coordinate system
     */
    class Mapping
    {
        public:
            SOLVER_UTILS_EXPORT virtual ~Mapping() {}

            /// Initialise the mapping object
            SOLVER_UTILS_EXPORT void InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>&       pFields,
                const TiXmlElement* pMapping);
            
            SOLVER_UTILS_EXPORT static MappingSharedPtr Load(
                        const LibUtilities::SessionReaderSharedPtr& pSession,
                        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields);
            
            /////////////////////////////////////////////////////////////
            //
            //    Functions for transforming results to and from the
            //          Cartesian coordinate system 
            //          (useful for pre and post processing)
            /////////////////////////////////////////////////////////////
              
            /// Convert a contravariant vector to the Cartesian system
            SOLVER_UTILS_EXPORT void ContravarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ContravarToCartesian( inarray, outarray);
            }
            
            /// Convert a covariant vector to the Cartesian system
            SOLVER_UTILS_EXPORT void CovarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_CovarToCartesian( inarray, outarray);
            }
            
            /// Convert a contravariant vector to the transformed system
            SOLVER_UTILS_EXPORT void ContravarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ContravarFromCartesian( inarray, outarray);
            }
            
            /// Convert a covariant vector to the transformed system
            SOLVER_UTILS_EXPORT void CovarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_CovarFromCartesian( inarray, outarray);
            }
            
            /////////////////////////////////////////////////////////////
            //
            //   Basic tensor calculus functions
            //    
            /////////////////////////////////////////////////////////////   

            /// Get the Jacobian of the transformation
            SOLVER_UTILS_EXPORT void GetJacobian(
                Array<OneD, NekDouble>               &outarray)
            {
                v_GetJacobian( outarray);
            } 
            
            /// Calculate the dot product with the Jacobian gradient
            SOLVER_UTILS_EXPORT void DotGradJacobian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>               &outarray)
            {
                v_DotGradJacobian( inarray, outarray);
            } 
            
            /// Get the metric tensor g_(i,j)
            SOLVER_UTILS_EXPORT void GetMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_GetMetricTensor( outarray);
            }  
            
            /// Get the inverse of metric tensor g^(i,j)
            SOLVER_UTILS_EXPORT void GetInvMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_GetInvMetricTensor( outarray);
            }
            
            /// Lower index using v_(i) = g_(i,j)*v^(j)
            SOLVER_UTILS_EXPORT void LowerIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_LowerIndex( inarray, outarray);
            }     
            
            /// Raise index using v^(i) = g^(i,j)*v_(j)
            SOLVER_UTILS_EXPORT void RaiseIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_RaiseIndex( inarray, outarray);
            }
            
            // Apply the Christoffel symbols to a contravariant vector
            //          outarray = {i,pk}*u^p
            SOLVER_UTILS_EXPORT void ApplyChristoffelContravar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ApplyChristoffelContravar( inarray, outarray);
            }
            
            // Apply the Christoffel symbols to a convariant vector
            //          outarray = {p,ik}*u_p
            SOLVER_UTILS_EXPORT void ApplyChristoffelCovar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ApplyChristoffelCovar( inarray, outarray);
            }

            
            /////////////////////////////////////////////////////////////
            //
            //   Functions used in IncNs
            //    
            /////////////////////////////////////////////////////////////             
            
            // Generalized divergence operator
            //         D = u^i_,i = 1/J*d(J*u^i)/dx^i
            SOLVER_UTILS_EXPORT void Divergence(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>                            &outarray)
            {
                v_Divergence( inarray, outarray);
            }  
            
            // Generalized velocity Laplacian operator
            //         L = g^jk u^i_{,jk}
            SOLVER_UTILS_EXPORT void VelocityLaplacian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_VelocityLaplacian( inarray, outarray);
            } 
            
            // Generalized velocity second order derivatives
            //         ddU = u^i_{,jk}
            SOLVER_UTILS_EXPORT void gradgradU(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_gradgradU( inarray, outarray);
            }
            
            // Correction needed for convective terms = N(u) - ( -(u \nabla) u)
            //     inarray is the velocity field
            SOLVER_UTILS_EXPORT void IncNSAdvectionCorrection(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_IncNSAdvectionCorrection( inarray, outarray);
            } 
            
            // Correction needed for pressure terms   
            //     = -g^(ij)p_j + (\nabla p)/J for variable Jacobian
            //     = -g^(ij)p_j + (\nabla p)   for constant Jacobian
            //         inarray is the pressure field
            SOLVER_UTILS_EXPORT void IncNSPressureCorrection(
                const Array<OneD, NekDouble>                      &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_IncNSPressureCorrection( inarray, outarray);
            }
            
            // Correction needed for viscous terms = g^jk u^i_{,jk}-(\nabla^2 u)
            //     inarray is the velocity field
            SOLVER_UTILS_EXPORT void IncNSViscousCorrection(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_IncNSViscousCorrection( inarray, outarray);
            }
            
            // CurlCurl calculated on the whole field
            //     when the mapping viscous are treated: 
            //          - explicitly it corresponds to the Cartesian curl-curl
            //          - implicitly it corresponds to the Generalized curl-curl
                    
            SOLVER_UTILS_EXPORT void CurlCurlField(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_CurlCurlField( inarray, outarray);
            }
            
            
            /////////////////////////////////////////////////////////////
            //
            //   Functions defining mapping properties
            //    
            /////////////////////////////////////////////////////////////             
            
            // Define if mapping is constant or time-dependent
            SOLVER_UTILS_EXPORT bool IsTimeDependent()
            {
                return v_IsTimeDependent();
            }          
            
            // Define if the Jacobian of the transformation is constant
            SOLVER_UTILS_EXPORT bool HasConstantJacobian()
            {
                return v_HasConstantJacobian();
            }
            
            // Define if the mapping pressure terms should be treated implicitly 
            SOLVER_UTILS_EXPORT bool ImplicitPressure()
            {
                return m_implicitPressure;
            }
            
            // Define if the mapping viscous terms should be treated implicitly
            SOLVER_UTILS_EXPORT bool ImplicitViscous()
            {
                return m_implicitViscous;
            }
            
            // Define relaxation paramter and tolerance for pressure and viscous system
            SOLVER_UTILS_EXPORT NekDouble PressureTolerance()
            {
                return m_pressureTolerance;
            }
            
            SOLVER_UTILS_EXPORT NekDouble ViscousTolerance()
            {
                return m_viscousTolerance;
            }
            
            SOLVER_UTILS_EXPORT NekDouble PressureRelaxation()
            {
                return m_pressureRelaxation;
            }
            
            SOLVER_UTILS_EXPORT NekDouble ViscousRelaxation()
            {
                return m_viscousRelaxation;
            }
            
            //
            //  Function to update time-dependent mappings
            //
            SOLVER_UTILS_EXPORT void UpdateMapping()
            {
                v_UpdateMapping();
            }
            

        protected:
            /// Session reader
            LibUtilities::SessionReaderSharedPtr        m_session;
            /// Fields
            Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
            // Arrays with geometric parameters of the mapping
            Array<OneD, Array<OneD, NekDouble> >        m_GeometricInfo;
            // Number of velocity components
            int                                         m_nConvectiveFields;
            // Flags defining if pressure and viscous mapping terms 
            //should be treated implicitly
            bool                                        m_implicitPressure;
            bool                                        m_implicitViscous;
            // Tolerance and relaxation parameters for pressure and viscous
            //       systems (when solved iteratively)
            NekDouble                                   m_pressureTolerance;
            NekDouble                                   m_viscousTolerance;
            NekDouble                                   m_pressureRelaxation;
            NekDouble                                   m_viscousRelaxation;
            
            

            /// Constructor
            SOLVER_UTILS_EXPORT Mapping(
                const LibUtilities::SessionReaderSharedPtr&,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields);

            // Evaluators
            SOLVER_UTILS_EXPORT void EvaluateFunction(
                    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
                    LibUtilities::SessionReaderSharedPtr        pSession,
                    std::string                                 pFieldName, 
                    Array<OneD, NekDouble>&                     pArray,
                    const std::string& pFunctionName,
                    NekDouble pTime = NekDouble(0));

            SOLVER_UTILS_EXPORT void EvaluateTimeFunction(
                    LibUtilities::SessionReaderSharedPtr        pSession,
                    std::string                                 pFieldName, 
                    Array<OneD, NekDouble>&                     pArray,
                    const std::string&                          pFunctionName,
                    NekDouble pTime = NekDouble(0));            
            
            // Virtual functions
            SOLVER_UTILS_EXPORT virtual void v_InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields,
                const TiXmlElement* pMapping) = 0;

            SOLVER_UTILS_EXPORT virtual void v_ContravarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;

            SOLVER_UTILS_EXPORT virtual void v_CovarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;            

            SOLVER_UTILS_EXPORT virtual void v_ContravarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;

            SOLVER_UTILS_EXPORT virtual void v_CovarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0; 
            
            SOLVER_UTILS_EXPORT virtual void v_GetJacobian(
                Array<OneD, NekDouble>               &outarray)             =0;
            
            SOLVER_UTILS_EXPORT virtual void v_DotGradJacobian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>               &outarray)            =0;
            
            SOLVER_UTILS_EXPORT virtual void v_GetMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            SOLVER_UTILS_EXPORT virtual void v_GetInvMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            SOLVER_UTILS_EXPORT virtual void v_LowerIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            SOLVER_UTILS_EXPORT virtual void v_RaiseIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;

            SOLVER_UTILS_EXPORT virtual void v_ApplyChristoffelContravar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            SOLVER_UTILS_EXPORT virtual void v_ApplyChristoffelCovar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            SOLVER_UTILS_EXPORT virtual void v_Divergence(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>                            &outarray); 
            
            SOLVER_UTILS_EXPORT virtual void v_VelocityLaplacian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);  
            
            SOLVER_UTILS_EXPORT virtual void v_gradgradU(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);
            
            SOLVER_UTILS_EXPORT virtual void v_IncNSAdvectionCorrection(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);
            
            SOLVER_UTILS_EXPORT virtual void v_IncNSPressureCorrection(
                const Array<OneD, NekDouble>                      &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);
            
            SOLVER_UTILS_EXPORT virtual void v_IncNSViscousCorrection(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);

            SOLVER_UTILS_EXPORT virtual void v_CurlCurlField(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);
            
            SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent() =0;  
            
            SOLVER_UTILS_EXPORT virtual bool v_HasConstantJacobian() =0;
            
            SOLVER_UTILS_EXPORT virtual void v_UpdateMapping() =0;
            
    };
}
}

#endif
