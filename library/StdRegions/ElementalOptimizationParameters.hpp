///////////////////////////////////////////////////////////////////////////////
//
// File ElementalOptimizationParameters.hpp
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

#ifndef NEKTAR_LIB_STDREGIONS_ELEMENTALOPTIMIZATIONPARAMETERS_H
#define NEKTAR_LIB_STDREGIONS_ELEMENTALOPTIMIZATIONPARAMETERS_H

#include <StdRegions/ElementalOptimizationParametersTraits.hpp>
#include <iostream>
#include <fstream>


namespace Nektar
{
    namespace NekOptimize
    {
        // enumeration of all operations that are optimisable
        enum ElementalOptimizationOperationType
        {
            eBwdTrans,
            eIProductWRTBase,
            eIProductWRTDerivBase,
            eMassMatrixOp,
            eLaplacianMatrixOp,
            eLaplacianMatrixIJOp,
            eWeakDerivMatrixOp,
            eHelmholtzMatrixOp,
            SIZE_OptimizeOperationType
        };

        const char* const ElementalOptimizationOperationTypeMap[] =
        {
            "BwdTrans",
            "IProductWRTBase",
            "IProductWRTDerivBase",
            "MassMatrixOp",
            "LaplacianMatrixOp",
            "LaplacianMatrixIJOp",
            "WeakDerivMatrixOp",
            "HelmholtzMatrixOp"
        };

        // This is the class which holds the optimisation parameters
        // for a specific element and method (which are the template parameters)
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype, int d>
        class ElementalOptimizationParametersBase
        {
        public:
            typedef typename ElementalOptimizationParametersTraits<etype>::DataHolderType DataHolderType;

            void DumpParameters(std::ostream &outfile) const
            {
                outfile << "OPTIMIZATION PARAMETERS FOR:" << std::endl;
                outfile << "   Element Type  : " <<  StdRegions::ElementTypeMap[etype] << std::endl;
                outfile << "   Operation Type: " <<  ElementalOptimizationOperationTypeMap[mtype] << std::endl;

                const int dim = ElementalOptimizationParametersTraits<etype>::elementDim;
                ElementalOptimizationParametersDumpPolicy<dim>::DumpElementalOptimizationParameters(outfile,m_data);

                outfile << std::endl;
            }

        protected:
            DataHolderType m_data;

            // This Loki function should be a friend, as only this function
            // can access the private default constructor (used to create the singleton)
            template <typename T> friend struct Loki::CreateUsingNew;

            // The class below should be made friend as well as it is used to load the
            // parameters
            //template<StdRegions::ElementType, ElementalOptimizationOperationType, int> friend class ElementalOptimization;

            // default constructor
            ElementalOptimizationParametersBase():
                m_data()
            {
            }
         };


        // --------------------------------------------------------------------
        // Specialisation for dim = 2
        // --------------------------------------------------------------------
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype, int dim = 2>
        class ElementalOptimizationParameters : public ElementalOptimizationParameters<etype, mtype, dim> {};

        // This is the class which holds the optimisation parameters
        // for a specific element and method (which are the template parameters)
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype>
        class ElementalOptimizationParameters<etype, mtype, 2>
                : public ElementalOptimizationParametersBase<etype, mtype, 2>
        {
        public:
            typedef typename ElementalOptimizationParametersTraits<etype>::DataHolderType DataHolderType;

            bool DoMatOp(int nummodes0, int nummodes1) const
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");

                return ElementalOptimizationParametersBase<etype, mtype, 2>::m_data.GetValue(nummodes0-2,nummodes1-2);
            }

        private:
            // This Loki function should be a friend, as only this function
            // can access the private default constructor (used to create the singleton)
            template <typename T> friend struct Loki::CreateUsingNew;

            // The class below should be made friend as well as it is used to load the
            // parameters
            template<StdRegions::ElementType, ElementalOptimizationOperationType, int> friend class ElementalOptimization;

            // default constructor
            ElementalOptimizationParameters()
            {
            }

            void SetDoMatOp(int nummodes0, int nummodes1, bool a)
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");

                ElementalOptimizationParametersBase<etype, mtype, 2>::m_data.SetValue(nummodes0-2,nummodes1-2,a);
            }
        };

        // This is the class which holds the optimisation parameters
        // for a specific element and method (which are the template parameters)
        template<StdRegions::ElementType etype, ElementalOptimizationOperationType mtype>
        class ElementalOptimizationParameters<etype, mtype, 3>
                : public ElementalOptimizationParametersBase<etype, mtype, 3>
        {
        public:
            typedef typename ElementalOptimizationParametersTraits<etype>::DataHolderType DataHolderType;

            bool DoMatOp(int nummodes0, int nummodes1, int nummodes2) const
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");
                ASSERTL1((nummodes2>1)&&(nummodes2<MaxBoolContainerDim+2),"Invalid number of modes in direction 2");

                return ElementalOptimizationParametersBase<etype, mtype, 3>::m_data.GetValue(nummodes0-2,nummodes1-2,nummodes2-2);
            }

        private:
            // This Loki function should be a friend, as only this function
            // can access the private default constructor (used to create the singleton)
            template <typename T> friend struct Loki::CreateUsingNew;

            // The class below should be made friend as well as it is used to load the
            // parameters
            template<StdRegions::ElementType, ElementalOptimizationOperationType, int> friend class ElementalOptimization;

            // default constructor
            ElementalOptimizationParameters()
            {
            }

            void SetDoMatOp(int nummodes0, int nummodes1, int nummodes2, bool a)
            {
                ASSERTL1((nummodes0>1)&&(nummodes0<MaxBoolContainerDim+2),"Invalid number of modes in direction 0");
                ASSERTL1((nummodes1>1)&&(nummodes1<MaxBoolContainerDim+2),"Invalid number of modes in direction 1");
                ASSERTL1((nummodes2>1)&&(nummodes2<MaxBoolContainerDim+2),"Invalid number of modes in direction 2");

                ElementalOptimizationParametersBase<etype, mtype, 3>::m_data.SetValue(nummodes0-2,nummodes1-2,nummodes2-2,a);
            }
        };

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_ELEMENTALOPTIMIZATIONPARAMETERS_H


