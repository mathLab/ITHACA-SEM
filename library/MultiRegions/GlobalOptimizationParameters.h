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

#include <StdRegions/ElementalOptimizationParameters.hpp>

namespace Nektar
{
    namespace NekOptimize
    {
        /// Processes global optimisation parameters from a session.
        class GlobalOptParam
        {
        public:
            /// Use default values.
            GlobalOptParam();

            /// Read optimisation parameters from a given file.
            GlobalOptParam(const std::string& fileName);

            /// For a given matrix type, determines if the operation should
            /// be done globally.
            // inline
            const bool DoGlobalMatOp(const StdRegions::MatrixType i) const;

            /// For a given matrix type, determines if the operation should be
            /// done
            // inline
            const bool DoBlockMatOp(const StdRegions::MatrixType i) const;

        private:
            /// Flags indicating if different matrices should be evaluated
            /// globally.
            Array<OneD,bool> m_doGlobalMatOp;

            /// Flags indicating if different matrices should be evaluated
            /// in block form.
            Array<OneD,bool> m_doBlockMatOp;
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
        inline const bool GlobalOptParam::DoGlobalMatOp(
                                const StdRegions::MatrixType i) const
        {
            ElementalOptimizationOperationType type;
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
        inline const bool GlobalOptParam::DoBlockMatOp(
                                const StdRegions::MatrixType i) const
        {
            ElementalOptimizationOperationType type;
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
            default:
                {
                    ASSERTL0(false,"Optimisation suite not set up for this type"
                                   " of matrix");
                }
            }
            return m_doBlockMatOp[type];
        }


    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H


