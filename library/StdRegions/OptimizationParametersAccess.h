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
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace NekOptimize
    {
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype, int dim = 2>
        class ElementalOptimization;

        // This class can be used to access the elemental optimisation parameters.
        // No direct access to the ElementalOptimizationParameters class is allowed as
        // they should be accessed through a Loki::Singleton
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype>
            class ElementalOptimization<etype, mtype, 2>
        {
        public:
            static bool DoMatOp(int nummodes0, int nummodes1)
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");

                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype,2> >::Instance()).DoMatOp(nummodes0,nummodes1);
            }

            static void DumpParameters(std::ostream &outfile)
            {
                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype,2> >::Instance()).DumpParameters(outfile);
            }

        private:
            friend class LoadOptimizationParametersInterface;
            static void SetDoMatOp(int nummodes0, int nummodes1, bool a)
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");

                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype,2> >::Instance()).SetDoMatOp(nummodes0,nummodes1,a);
            }
        };

        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype>
            class ElementalOptimization<etype, mtype, 3>
        {
        public:
            static bool DoMatOp(int nummodes0, int nummodes1, int nummodes2)
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");
                ASSERTL1((nummodes2>1)&&(nummodes2<MaxBoolContainerDim+2),"Invalid number of modes in direction 2");

                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype,3> >::Instance()).DoMatOp(nummodes0,nummodes1,nummodes2);
            }

            static void DumpParameters(std::ostream &outfile)
            {
                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype,3> >::Instance()).DumpParameters(outfile);
            }

        private:
            friend class LoadOptimizationParametersInterface;
            static void SetDoMatOp(int nummodes0, int nummodes1, int nummodes2, bool a)
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");
                ASSERTL1((nummodes2>1)&&(nummodes2<MaxBoolContainerDim+2),"Invalid number of modes in direction 2");

                return (Loki::SingletonHolder<ElementalOptimizationParameters<etype,mtype,3> >::Instance()).SetDoMatOp(nummodes0,nummodes1,nummodes2,a);
            }
        };

        // This is the global function which can be used by the user to load all the
        // optimisation parameters based on an input file
        STD_REGIONS_EXPORT void LoadElementalOptimizationParameters(const std::string& fileName);

        // This class is merely an interface which collects all the different
        // loading routines
        class LoadOptimizationParametersInterface
        {
        private:
                STD_REGIONS_EXPORT friend void LoadElementalOptimizationParameters(const std::string& fileName);

                STD_REGIONS_EXPORT static void LoadElemental2DOptimizationParameters(const std::string& fileName);
        };

        // This is the global function which can be used by the user to dump all the
        // optimisation parameters to a stream
        STD_REGIONS_EXPORT void DumpElementalOptimizationParameters(std::ostream &outfile);

        // This class is merely an interface which collects all the different
        // dumping routines
        class DumpOptimizationParametersInterface
        {
        private:
                STD_REGIONS_EXPORT friend void DumpElementalOptimizationParameters(std::ostream &outfile);

                STD_REGIONS_EXPORT static void DumpElemental2DOptimizationParameters(std::ostream &outfile);
        };

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_OPTIMIZATIONPARAMETERSACCESS_H


