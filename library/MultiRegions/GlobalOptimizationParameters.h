///////////////////////////////////////////////////////////////////////////////
//
// File GlobalOptimizationParameters.h
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
// Description: Header file of global optimisation parameters class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALOPTIMIZATIONPARAMETERS_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALOPTIMIZATIONPARAMETERS_H
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <StdRegions/StdRegions.hpp>
#include <MultiRegions/MultiRegionsDeclspec.h>

namespace Nektar
{
    namespace NekOptimize
    {
        // enumeration of all operations that are optimisable
        enum OptimizationOperationType
        {
            eBwdTrans,
            eIProductWRTBase,
            eIProductWRTDerivBase,
            eMassMatrixOp,
            eLaplacianMatrixOp,
            eLaplacianMatrixIJOp,
            eWeakDerivMatrixOp,
            eHelmholtzMatrixOp,
            eHybridDGHelmBndLamMatrixOp,
            SIZE_OptimizeOperationType
        };

        const char* const OptimizationOperationTypeMap[] =
        {
            "BwdTrans",
            "IProductWRTBase",
            "IProductWRTDerivBase",
            "MassMatrixOp",
            "LaplacianMatrixOp",
            "LaplacianMatrixIJOp",
            "WeakDerivMatrixOp",
            "HelmholtzMatrixOp",
            "HybridDGHelmBndLamMatrixOp"
        };

        /// Processes global optimisation parameters from a session.
        class GlobalOptParam
        {
        public:
            /// Default constructor requires nel be given 
            MULTI_REGIONS_EXPORT GlobalOptParam(const int nel);
            
            /// Read optimisation parameters from a given file.
            MULTI_REGIONS_EXPORT GlobalOptParam(const LibUtilities::SessionReaderSharedPtr& pSession, const int dim, const Array<OneD, const int> &NumShapeElements);

            /// For a given matrix type, determines if the operation should
            /// be done globally.
            // inline
            MULTI_REGIONS_EXPORT bool DoGlobalMatOp(const StdRegions::MatrixType i) const;
            

            /// For a given matrix type, determines if the operation should be
            /// done with a block matrix
            // inline
            inline const Array<OneD, const bool>  &DoBlockMatOp(const StdRegions::MatrixType i) const;
            
            inline const Array<OneD, const LibUtilities::ShapeType>  &GetShapeList() const;
            inline const Array<OneD, const int>  &GetShapeNumElements() const; 

        private:
            /// Default constructor should not be called
            GlobalOptParam() {};

            /// Flags indicating if different matrices should be evaluated
            /// globally.
            Array<OneD,bool> m_doGlobalMatOp;

            /// Array of Flags of first dimension of the number of
            /// shapes within the space dimension, indicating if
            /// different matrices should be evaluated using a block
            /// matrix
            Array<OneD, Array<OneD,bool> > m_doBlockMatOp; 

            /// A list ExpansionTypes indicating the order in which
            /// shapes are listed to call the appropriate key for the
            /// block matrices.
            Array<OneD, LibUtilities::ShapeType> m_shapeList;

            /// A list of  number of elements contained within each shape type
            Array<OneD, const int> m_shapeNumElements;
        };

        /// Pointer to a GlobalOptParam object.
        typedef  boost::shared_ptr<GlobalOptParam> GlobalOptParamSharedPtr;
        
        /// Pointer to an empty GlobalOptParam object.
        static   GlobalOptParamSharedPtr NullGlobalOptParamSharedPtr;


        /**
         * Determines the elemental optimisation type enum, given the
         * MatrixType and returns the corresponding entry in the table.
         * @param   i           Type of matrix.
         * @returns True if this type of matrix should be evaluated globally.
         */
        inline bool GlobalOptParam::DoGlobalMatOp(const StdRegions::MatrixType i) const
        {
            OptimizationOperationType type;
            switch(i)
            {
            case StdRegions::eBwdTrans:
                {
                    type = eBwdTrans;
                }
                break;
            case StdRegions::eIProductWRTBase:
                {
                    type = eIProductWRTBase;
                }
                break;
            case StdRegions::eMass:
                {
                    type = eMassMatrixOp;
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    type = eHelmholtzMatrixOp;
                }
                break;
            case StdRegions::eLaplacian:
                {
                    type = eLaplacianMatrixOp;
                }
                break;
            case StdRegions::eHybridDGHelmBndLam:
                {
                    type = eHybridDGHelmBndLamMatrixOp;
                }
                break;
            default:
                {
                    ASSERTL0(false,"Optimisation suite not set up for this type"
                                   " of matrix");
                }
            }
            return m_doGlobalMatOp[type];
        }

        /**
         * Determines the elemental optimisation type enum, given the
         * MatrixType and returns the corresponding entry in the table.
         * @param   i           Type of matrix.
         * @returns True if this type of matrix should be evaluated in block
         *          form.
         */
        inline  const Array<OneD, const bool> &GlobalOptParam::DoBlockMatOp(const StdRegions::MatrixType i) const
        {
            OptimizationOperationType type;
            switch(i)
            {
            case StdRegions::eBwdTrans:
                {
                    type = eBwdTrans;
                }
                break;
            case StdRegions::eIProductWRTBase:
                {
                    type = eIProductWRTBase;
                }
                break;
            case StdRegions::eMass:
                {
                    type = eMassMatrixOp;
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    type = eHelmholtzMatrixOp;
                }
                break;
            case StdRegions::eLaplacian:
                {
                    type = eLaplacianMatrixOp;
                }
                break;
            case StdRegions::eHybridDGHelmBndLam:
                {
                    type = eHybridDGHelmBndLamMatrixOp;
                }
                break;
            default:
                {
                    ASSERTL0(false,"Optimisation suite not set up for this type"
                                   " of matrix");
                }
            }

            return m_doBlockMatOp[type];
        }

        inline const Array<OneD, const int>  &GlobalOptParam::GetShapeNumElements() const
        {
            return m_shapeNumElements;
        }

        inline const Array<OneD, const LibUtilities::ShapeType>  &GlobalOptParam::GetShapeList() const
        {
            return m_shapeList;
        }



    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H


