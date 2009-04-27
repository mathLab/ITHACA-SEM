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

#include<StdRegions/ElementalOptimizationParameters.hpp>

namespace Nektar
{
    namespace NekOptimize 
    {        
        class GlobalOptParam 
        {
        public:

            GlobalOptParam(void);

            GlobalOptParam(const std::string& fileName);

            inline const bool DoGlobalMatOp(const StdRegions::MatrixType i) const
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
                        ASSERTL0(false,"Optimisation suite not set up for this type of matrix");
                    }
                }
                return m_doGlobalMatOp[type];
            }
            
            inline const bool DoBlockMatOp(const StdRegions::MatrixType i) const
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
                        ASSERTL0(false,"Optimisation suite not set up for this type of matrix");
                    }
                }
                return m_doBlockMatOp[type];
            }
            
        private:
            Array<OneD,bool> m_doGlobalMatOp;
            Array<OneD,bool> m_doBlockMatOp;
        };    

        typedef  boost::shared_ptr<GlobalOptParam> GlobalOptParamSharedPtr;
        static   GlobalOptParamSharedPtr NullGlobalOptParamSharedPtr;
        
    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H


