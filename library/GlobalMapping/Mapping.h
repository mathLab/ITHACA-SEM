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
GLOBAL_MAPPING_EXPORT typedef std::shared_ptr<Mapping> MappingSharedPtr;

/// Declaration of the mapping factory
typedef LibUtilities::NekFactory<std::string, Mapping,
        const LibUtilities::SessionReaderSharedPtr&,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&,
        const TiXmlElement*> MappingFactory;

/// Declaration of the mapping factory singleton
GLOBAL_MAPPING_EXPORT MappingFactory& GetMappingFactory();

/**
 * @class Mapping
 * @brief Base class for mapping to be applied to the coordinate system
 */
class Mapping
{
    public:
        /// @brief Destructor
        GLOBAL_MAPPING_EXPORT virtual ~Mapping() {}

        /// @brief Initialise the mapping object
        GLOBAL_MAPPING_EXPORT void InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>&     pFields,
            const TiXmlElement* pMapping)
        {
            v_InitObject( pFields, pMapping);
        }

        /// @brief Replace the Expansion List used by the mapping
        GLOBAL_MAPPING_EXPORT void ReplaceField(
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields);

        /// @brief Return a pointer to the mapping, creating it on first call 
        GLOBAL_MAPPING_EXPORT static MappingSharedPtr Load(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields);

        /// @brief Output function called when a chk or fld file is written
        GLOBAL_MAPPING_EXPORT void Output( 
                LibUtilities::FieldMetaDataMap & fieldMetaDataMap,
                const std::string                    &outname);

        /////////////////////////////////////////////////////////////
        //
        //    Functions for transforming results to and from the
        //          Cartesian coordinate system 
        //          (useful for pre and post processing)
        /////////////////////////////////////////////////////////////

        /// @brief Convert a contravariant vector to the Cartesian system
        GLOBAL_MAPPING_EXPORT void ContravarToCartesian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        /// @brief Convert a covariant vector to the Cartesian system
        GLOBAL_MAPPING_EXPORT void CovarToCartesian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        /// @brief Convert a contravariant vector to the transformed system
        GLOBAL_MAPPING_EXPORT void ContravarFromCartesian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        /// @brief Convert a covariant vector to the transformed system
        GLOBAL_MAPPING_EXPORT void CovarFromCartesian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        /**
         * @brief Get the Cartesian coordinates in the field 
         * 
         * This function is used to obtain the Cartesian coordinates 
         * associated withthe Mapping
         *  
         * @param out0  Coordinates in the x-direction 
         * @param out1  Coordinates in the y-direction 
         * @param out2  Coordinates in the z-direction 
         */
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

        /**
         * @brief Get the Jacobian of the transformation
         * 
         * This function is used to obtain the Jacobian of the Mapping.
         *  
         * @param outarray Array containing the Jacobian
         */
        GLOBAL_MAPPING_EXPORT void GetJacobian(
            Array<OneD, NekDouble>               &outarray)
        {
            v_GetJacobian( outarray);
        } 

        /**
         * @brief Calculate the dot product with the gradient of the Jacobian
         * 
         * This function calculates the dot product of an array against the
         * gradient of the Jacobian of the Mapping.
         *  
         * @param inarray  Input array
         * @param outarray Output array
         */
        GLOBAL_MAPPING_EXPORT void DotGradJacobian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, NekDouble>                            &outarray)
        {
            v_DotGradJacobian( inarray, outarray);
        } 

        /// @brief Get the metric tensor \f$g_{ij}\f$
        GLOBAL_MAPPING_EXPORT void GetMetricTensor(
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            v_GetMetricTensor( outarray);
        }  

        /// @brief Get the inverse of metric tensor \f$g^{ij}\f$
        GLOBAL_MAPPING_EXPORT void GetInvMetricTensor(
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            v_GetInvMetricTensor( outarray);
        }

        /// @brief Lower index of vector: \f$v_{i} = g_{ij}*v^{j}\f$
        GLOBAL_MAPPING_EXPORT void LowerIndex(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);  

        /// @brief Raise index of vector: \f$v^{i} = g^{ij}*v_{j}\f$
        GLOBAL_MAPPING_EXPORT void RaiseIndex(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        /**
         * @brief Apply the Christoffel symbols to a contravariant vector
         * 
         * This function is used apply the Christoffel symbols 
         * \f$ \left(  i,pk\right)\f$ to
         * a contravariant vector \f$u^p\f$, following the relation
         * \f[ (out)^{i}_{k} = \left(  i,pk\right)u^p \f]
         *  
         * @param inarray  Contravariant vector \f$u^p\f$
         * @param outarray Result of applying Christoffel symbols to \f$u^p\f$
         */
        GLOBAL_MAPPING_EXPORT void ApplyChristoffelContravar(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            v_ApplyChristoffelContravar( inarray, outarray);
        }

        /**
         * @brief Apply the Christoffel symbols to a covariant vector
         * 
         * This function is used apply the Christoffel symbols 
         * \f$ \left(  p,ik\right)\f$ to
         * a covariant vector \f$u_p\f$, following the relation
         * \f[ (out)_{ik} = \left(  p,ik\right)u_p \f]
         *  
         * @param inarray  Contravariant vector \f$u_p\f$
         * @param outarray Result of applying Christoffel symbols to \f$u_p\f$
         */
        GLOBAL_MAPPING_EXPORT void ApplyChristoffelCovar(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            v_ApplyChristoffelCovar( inarray, outarray);
        }

        /**
         * @brief Obtain the velocity of the coordinates 
         * 
         * This function is used to obtain the velocity of the coordinates
         * associated with the Mapping
         *  
         * @param outarray  Velocity of the coordinates
         */
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

        /**
         * @brief Calculate the generalised divergence operator
         * 
         * This function is used to calculate the generalised divergence
         * of a contravariant vector, defined as
         * \f[ D = u^i_{,i} = \frac{1}{J}\frac{\partial(Ju^i)}{\partial x^i}\f]
         *  
         * @param inarray   Contravariant vector \f$u^i\f$
         * @param outarray  Divergence of \f$u^i\f$
         */
        GLOBAL_MAPPING_EXPORT void Divergence(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, NekDouble>                            &outarray)
        {
            v_Divergence( inarray, outarray);
        }  

        /**
         * @brief Generalised (correction to the) velocity Laplacian operator
         * 
         * This function is used to calculate a correction defined as
         * the difference between the generalised Laplacian and the 
         * original Laplacian multiplied by a constant \f$\alpha\f$,
         *  resulting in
         * \f[ L^i = g^{jk}u^{i}_{,jk} - \alpha \frac{\partial^2 x^i}
         *           {\partial x^j \partial x^j}\f]
         * By default, \f$\alpha\f$ is zero, resulting in the generalised
         * Laplacian.
         *  
         * @param inarray   Contravariant vector \f$u^i\f$
         * @param outarray  Result of the operation
         * @param alpha     The constant \f$\alpha\f$
         */
        GLOBAL_MAPPING_EXPORT void VelocityLaplacian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            const NekDouble                                   alpha = 0.0)
        {
            v_VelocityLaplacian( inarray, outarray, alpha);
        } 

        /**
         * @brief Second order covariant derivatives of a contravariant vector
         * 
         * This function computes the second order covariant derivatives
         * of a contravariant vector, resulting in \f$u^{i}_{,jk}\f$ 
         * 
         * @param inarray   Contravariant vector \f$u^i\f$
         * @param outarray  Second order derivatives \f$u^{i}_{,jk}\f$ 
         */
        GLOBAL_MAPPING_EXPORT void gradgradU(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            v_gradgradU( inarray, outarray);
        }

        /**
         * @brief CurlCurl calculated on the whole field
         * 
         * This function can be used to compute both the generalised CurlCurl
         * or the typical (Cartesian) CurlCurl, depending on the flag 
         * generalized
         * 
         * @param inarray      Contravariant vector \f$u^i\f$
         * @param outarray     CurlCurl of \f$u\f$ 
         * @param generalized  Flag defining if generalised or typical CurlCurl
         */
        GLOBAL_MAPPING_EXPORT void CurlCurlField(
            Array<OneD, Array<OneD, NekDouble> >              &inarray,
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

        /// @brief Get flag defining if mapping is time-dependent
        GLOBAL_MAPPING_EXPORT bool IsTimeDependent()
        {
            return m_timeDependent;
        } 

        /// @brief Set flag defining if mapping is time-dependent
        GLOBAL_MAPPING_EXPORT void SetTimeDependent( const bool value)
        {
            m_timeDependent = value;
        }

        /// @brief Get flag defining if mapping is defined by a function
        GLOBAL_MAPPING_EXPORT bool IsFromFunction()
        {
            return m_fromFunction;
        } 

        /// @brief Set flag defining if mapping is defined by a function
        GLOBAL_MAPPING_EXPORT void SetFromFunction( const bool value)
        {
            m_fromFunction = value;
        }

        /// @brief Get flag defining if mapping has constant Jacobian
        GLOBAL_MAPPING_EXPORT bool HasConstantJacobian()
        {
            return m_constantJacobian;
        }

        /// @brief Get flag determining if the mapping was defined or is trivial
        GLOBAL_MAPPING_EXPORT bool IsDefined()
        {
            return m_isDefined;
        }

        //
        //  Function to update time-dependent mappings
        //            

        /// @brief Update the Dirichlet Boundary Conditions when using Mappings
        GLOBAL_MAPPING_EXPORT void UpdateBCs( const NekDouble time)
        {
            v_UpdateBCs(time);
        }

        /// @brief Update the Mapping with new coordinates
        GLOBAL_MAPPING_EXPORT void UpdateMapping(const NekDouble time,
            const Array<OneD, Array<OneD, NekDouble> > &coords    
                                                    = NullNekDoubleArrayofArray,
            const Array<OneD, Array<OneD, NekDouble> > &coordsVel 
                                                    = NullNekDoubleArrayofArray)
        {
            v_UpdateMapping( time, coords, coordsVel);
        }

        /// @brief Recompute the metric terms of the Mapping
        GLOBAL_MAPPING_EXPORT void UpdateGeomInfo()
        {
            v_UpdateGeomInfo();
        }


    protected:
        /// Session reader
        LibUtilities::SessionReaderSharedPtr        m_session;
        // FieldIO object used to output mapping
        LibUtilities::FieldIOSharedPtr              m_fld;
        // Fields 
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
        /// Array with the Cartesian coordinates
        Array<OneD, Array<OneD, NekDouble> >        m_coords;
        /// Array with the velocity of the coordinates
        Array<OneD, Array<OneD, NekDouble> >        m_coordsVel;
        /// Array with metric terms of the mapping
        Array<OneD, Array<OneD, NekDouble> >        m_GeometricInfo;
        /// Number of velocity components
        int                                         m_nConvectiveFields;

        /// Name of the function containing the coordinates
        std::string                                 m_funcName;
        /// Name of the function containing the velocity of the coordinates
        std::string                                 m_velFuncName;


        // Flags to help the solver
        /// Flag defining if the Jacobian is constant
        bool                                        m_constantJacobian;
        /// Flag defining if the Mapping is time-dependent
        bool                                        m_timeDependent;
        /// Flag defining if the Mapping is defined by a function
        bool                                        m_fromFunction;

        // Static variables to load mapping
        static MappingSharedPtr                     m_mappingPtr;
        static bool                                 m_init;
        static bool                                 m_isDefined;
        
        // Workspace variables
        Array<OneD, Array<OneD, NekDouble> >        m_wk1;
        Array<OneD, Array<OneD, NekDouble> >        m_wk2;
        Array<OneD, Array<OneD, NekDouble> >        m_tmp;

        /// @brief Constructor
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
            Array<OneD, NekDouble>               &out2);

        GLOBAL_MAPPING_EXPORT virtual void v_GetCoordVelocity(
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

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

        GLOBAL_MAPPING_EXPORT virtual void v_Divergence(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, NekDouble>                            &outarray); 

        GLOBAL_MAPPING_EXPORT virtual void v_VelocityLaplacian(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            const NekDouble                                    alpha);

        GLOBAL_MAPPING_EXPORT virtual void v_gradgradU(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        GLOBAL_MAPPING_EXPORT virtual void v_CurlCurlField(
            Array<OneD, Array<OneD, NekDouble> >              &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            const bool                                         generalized);

        GLOBAL_MAPPING_EXPORT virtual void v_UpdateMapping(
            const NekDouble time,
            const Array<OneD, Array<OneD, NekDouble> > &coords    
                                                = NullNekDoubleArrayofArray,
            const Array<OneD, Array<OneD, NekDouble> > &coordsVel 
                                                = NullNekDoubleArrayofArray);

        GLOBAL_MAPPING_EXPORT virtual void v_UpdateGeomInfo() =0;

        GLOBAL_MAPPING_EXPORT virtual void v_UpdateBCs(const NekDouble time);

};

}
}

#endif
