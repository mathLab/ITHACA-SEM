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

#ifndef NEKTAR_GLOBALMAPPING_MAPPING
#define NEKTAR_GLOBALMAPPING_MAPPING

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <GlobalMapping/GlobalMappingDeclspec.h>

namespace Nektar
{
namespace GlobalMapping
{
    //  Forward declaration
    class Mapping;

    /// A shared pointer to a Mapping object
    GLOBAL_MAPPING_EXPORT typedef boost::shared_ptr<Mapping> MappingSharedPtr;

    /// Declaration of the mapping factory
    typedef LibUtilities::NekFactory<std::string, Mapping,
            const LibUtilities::SessionReaderSharedPtr&,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&,
            const TiXmlElement*> MappingFactory;

    /// Declaration of the mapping factory singleton
    GLOBAL_MAPPING_EXPORT MappingFactory& GetMappingFactory();

    /**
     * @class Mapping
     * @brief Defines a mapping to be applied to the coordinate system
     */
    class Mapping
    {
        public:
            GLOBAL_MAPPING_EXPORT virtual ~Mapping() {}

            /// Initialise the mapping object
            GLOBAL_MAPPING_EXPORT void InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>&     pFields,
                const TiXmlElement* pMapping)
            {
                v_InitObject( pFields, pMapping);
            }
            
            GLOBAL_MAPPING_EXPORT static MappingSharedPtr Load(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields);
            
            /////////////////////////////////////////////////////////////
            //
            //    Functions for transforming results to and from the
            //          Cartesian coordinate system 
            //          (useful for pre and post processing)
            /////////////////////////////////////////////////////////////
              
            /// Convert a contravariant vector to the Cartesian system
            GLOBAL_MAPPING_EXPORT void ContravarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ContravarToCartesian( inarray, outarray);
            }
            
            /// Convert a covariant vector to the Cartesian system
            GLOBAL_MAPPING_EXPORT void CovarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_CovarToCartesian( inarray, outarray);
            }
            
            /// Convert a contravariant vector to the transformed system
            GLOBAL_MAPPING_EXPORT void ContravarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ContravarFromCartesian( inarray, outarray);
            }
            
            /// Convert a covariant vector to the transformed system
            GLOBAL_MAPPING_EXPORT void CovarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_CovarFromCartesian( inarray, outarray);
            }
            
            /// Get the Cartesian coordinates in the field
            GLOBAL_MAPPING_EXPORT void GetCartesianCoordinates(
                Array<OneD, NekDouble>               &out0,
                Array<OneD, NekDouble>               &out1,
                Array<OneD, NekDouble>               &out2)
            {
                v_GetCartesianCoordinates( out0, out1, out2);
            }
            
            /////////////////////////////////////////////////////////////
            //
            //   Basic tensor calculus functions
            //    
            /////////////////////////////////////////////////////////////   

            /// Get the Jacobian of the transformation
            GLOBAL_MAPPING_EXPORT void GetJacobian(
                Array<OneD, NekDouble>               &outarray)
            {
                v_GetJacobian( outarray);
            } 
            
            /// Calculate the dot product with the Jacobian gradient
            GLOBAL_MAPPING_EXPORT void DotGradJacobian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>               &outarray)
            {
                v_DotGradJacobian( inarray, outarray);
            } 
            
            /// Get the metric tensor g_(i,j)
            GLOBAL_MAPPING_EXPORT void GetMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_GetMetricTensor( outarray);
            }  
            
            /// Get the inverse of metric tensor g^(i,j)
            GLOBAL_MAPPING_EXPORT void GetInvMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_GetInvMetricTensor( outarray);
            }
            
            /// Lower index using v_(i) = g_(i,j)*v^(j)
            GLOBAL_MAPPING_EXPORT void LowerIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_LowerIndex( inarray, outarray);
            }     
            
            /// Raise index using v^(i) = g^(i,j)*v_(j)
            GLOBAL_MAPPING_EXPORT void RaiseIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_RaiseIndex( inarray, outarray);
            }
            
            // Apply the Christoffel symbols to a contravariant vector
            //          outarray = {i,pk}*u^p
            GLOBAL_MAPPING_EXPORT void ApplyChristoffelContravar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ApplyChristoffelContravar( inarray, outarray);
            }
            
            // Apply the Christoffel symbols to a convariant vector
            //          outarray = {p,ik}*u_p
            GLOBAL_MAPPING_EXPORT void ApplyChristoffelCovar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_ApplyChristoffelCovar( inarray, outarray);
            }
            
            // Obtain the velocity of the coordinates for time-dependent mappings
            GLOBAL_MAPPING_EXPORT void GetCoordVelocity(
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_GetCoordVelocity( outarray);
            }
            
            /////////////////////////////////////////////////////////////
            //
            //   Differential operators
            //    
            /////////////////////////////////////////////////////////////             
            
            // Generalized divergence operator
            //         D = u^i_,i = 1/J*d(J*u^i)/dx^i
            GLOBAL_MAPPING_EXPORT void Divergence(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>                            &outarray)
            {
                v_Divergence( inarray, outarray);
            }  
            
            // Generalized velocity Laplacian operator
            //         L = g^jk u^i_{,jk}
            GLOBAL_MAPPING_EXPORT void VelocityLaplacian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_VelocityLaplacian( inarray, outarray);
            } 
            
            // Generalized velocity second order derivatives
            //         ddU = u^i_{,jk}
            GLOBAL_MAPPING_EXPORT void gradgradU(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
            {
                v_gradgradU( inarray, outarray);
            }
            
            // CurlCurl calculated on the whole field
            //     if the flag generalized is: 
            //        - false, this corresponds to the usual Cartesian operator
            //        - true,  it corresponds to the Generalized curl-curl
                    
            GLOBAL_MAPPING_EXPORT void CurlCurlField(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const bool                                        generalized)
            {
                v_CurlCurlField( inarray, outarray, generalized);
            }
            
            
            /////////////////////////////////////////////////////////////
            //
            //   Functions defining mapping properties
            //    
            /////////////////////////////////////////////////////////////             
            
            // Define if mapping is constant or time-dependent
            GLOBAL_MAPPING_EXPORT bool IsTimeDependent()
            {
                return m_timeDependent;
            }          
            
            // Define if the Jacobian of the transformation is constant
            GLOBAL_MAPPING_EXPORT bool HasConstantJacobian()
            {
                return m_constantJacobian;
            }
            
            //
            //  Function to update time-dependent mappings
            //
            GLOBAL_MAPPING_EXPORT void UpdateMapping(const NekDouble time)
            {
                v_UpdateMapping( time);
            }
            
            GLOBAL_MAPPING_EXPORT void UpdateBCs( const NekDouble time)
            {
                v_UpdateBCs(time);
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
            
            // Name of the function containing the coordinates
            string                                      m_funcName;
            
            // Flags to help the solver
            bool                                        m_constantJacobian;
            bool                                        m_timeDependent;
            
            // Static variables to load mapping
            static MappingSharedPtr                     m_mappingPtr;
            static bool                                 m_init;
            
            

            /// Constructor
            GLOBAL_MAPPING_EXPORT Mapping(
                const LibUtilities::SessionReaderSharedPtr&          pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields);

            // Evaluators
            GLOBAL_MAPPING_EXPORT void EvaluateFunction(
                    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
                    LibUtilities::SessionReaderSharedPtr        pSession,
                    std::string                                 pFieldName, 
                    Array<OneD, NekDouble>&                     pArray,
                    const std::string& pFunctionName,
                    NekDouble pTime = NekDouble(0));

            GLOBAL_MAPPING_EXPORT void EvaluateTimeFunction(
                    LibUtilities::SessionReaderSharedPtr        pSession,
                    std::string                                 pFieldName, 
                    Array<OneD, NekDouble>&                     pArray,
                    const std::string&                          pFunctionName,
                    NekDouble pTime = NekDouble(0));            
            
            // Virtual functions
            GLOBAL_MAPPING_EXPORT virtual void v_InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields,
                const TiXmlElement* pMapping);

            GLOBAL_MAPPING_EXPORT virtual void v_ContravarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;

            GLOBAL_MAPPING_EXPORT virtual void v_CovarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;            

            GLOBAL_MAPPING_EXPORT virtual void v_ContravarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;

            GLOBAL_MAPPING_EXPORT virtual void v_CovarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0; 
            
            GLOBAL_MAPPING_EXPORT virtual void v_GetCartesianCoordinates(
                Array<OneD, NekDouble>               &out0,
                Array<OneD, NekDouble>               &out1,
                Array<OneD, NekDouble>               &out2) =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_GetJacobian(
                Array<OneD, NekDouble>               &outarray)             =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_DotGradJacobian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>               &outarray);
            
            GLOBAL_MAPPING_EXPORT virtual void v_GetMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_GetInvMetricTensor(
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_LowerIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);
            
            GLOBAL_MAPPING_EXPORT virtual void v_RaiseIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);

            GLOBAL_MAPPING_EXPORT virtual void v_ApplyChristoffelContravar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_ApplyChristoffelCovar(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray) =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_GetCoordVelocity(
                Array<OneD, Array<OneD, NekDouble> >              &outarray);
            
            GLOBAL_MAPPING_EXPORT virtual void v_Divergence(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble>                            &outarray); 
            
            GLOBAL_MAPPING_EXPORT virtual void v_VelocityLaplacian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);  
            
            GLOBAL_MAPPING_EXPORT virtual void v_gradgradU(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);

            GLOBAL_MAPPING_EXPORT virtual void v_CurlCurlField(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const bool                                      generalized);
            
            GLOBAL_MAPPING_EXPORT virtual void v_UpdateMapping(
                                                const NekDouble time) =0;
            
            GLOBAL_MAPPING_EXPORT virtual void v_UpdateBCs(const NekDouble time);
            
    };
}
}

#endif
