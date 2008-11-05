///////////////////////////////////////////////////////////////////////////////
//
// File OptimizationParametersAccess.h
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
// Description: Header file of optimisation parameters class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H
#define NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H

#include <StdRegions/ElementalOptimizationParameters.hpp>


namespace Nektar
{
    namespace NekOptimize 
    {        
        // This class can be used to access the elemental optimisation parameters.
        // No direct access to the ElementalOptimizationParameters class is allowed as
        // they should be accessed through a Loki::Singleton
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype>
            class ElementalOptimization
        {
        public:
            static bool DoMatOp(int nummodes0, int nummodes1) 
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");

                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype> >::Instance()).DoMatOp(nummodes0,nummodes1);
            }
            
            static void DumpParameters(std::ostream &outfile)
            {
                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype> >::Instance()).DumpParameters(outfile);
            }

        private:
            friend class LoadOptimizationParametersInterface;
            static void SetDoMatOp(int nummodes0, int nummodes1, bool a) 
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");

                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype> >::Instance()).SetDoMatOp(nummodes0,nummodes1,a);
            }
        };

        // This is the global function which can be used by the user to load all the
        // optimisation parameters based on an input file
        void LoadOptimizationParameters(const std::string& fileName);

        // This class is merely an interface which collects all the different
        // loading routines
        class LoadOptimizationParametersInterface
        {
        private:
            friend void LoadOptimizationParameters(const std::string& fileName);

            static void LoadElemental2DOptimizationParameters(const std::string& fileName);
        };

        // This is the global function which can be used by the user to dump all the
        // optimisation parameters to a stream
        void DumpOptimizationParameters(std::ostream &outfile);

        // This class is merely an interface which collects all the different
        // dumping routines
        class DumpOptimizationParametersInterface
        {
        private:
            friend void DumpOptimizationParameters(std::ostream &outfile);

            static void DumpElemental2DOptimizationParameters(std::ostream &outfile);
        };

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H


