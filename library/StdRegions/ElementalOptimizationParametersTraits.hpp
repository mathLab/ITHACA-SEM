///////////////////////////////////////////////////////////////////////////////
//
// File ElementalOptimizationParametersTraits.hpp
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

#ifndef NEKTAR_LIB_STDREGIONS_ELEMENTALOPTIMIZATIONPARAMETERSTRAITS_H
#define NEKTAR_LIB_STDREGIONS_ELEMENTALOPTIMIZATIONPARAMETERSTRAITS_H

#include <StdRegions/StdRegions.hpp>
#include <iomanip>
#include <iostream>
#include <fstream>

namespace Nektar
{
    namespace NekOptimize
    {

        // ==================================================================
        // Class to hold the boolean values of the optimisation parameters
        // ============
        // if - value = true :  matrix operations
        //    - value = false : dependending on the context:
        //                         - sum-factorisation
        //                         - partitioned evaluation of the operation
        //                         - elemental evaluation (rather than global matrix op)
        // Currently, this container is simply implemented as a multi-dimensional
        // array of bools. When memory might become an issue (too much parameters to store, e.g. 3D),
        // an option might be to replace this by a mechanism which stores boolean value as a bit.
        // (e.g. stl::bitset  -  Although this consumes less memory, its access is slower)
        template<int dim>
        class BoolContainer;

        static const int MaxBoolContainerDim = 18;

        // Implementation for 2D expansions
        template<>
        class BoolContainer<2>
        {
        public:
            BoolContainer()
            {
                std::fill_n(&m_data[0][0],MaxBoolContainerDim*MaxBoolContainerDim,false);
            }

            void SetValue(int i, int j, bool a)
            {
                ASSERTL1((i>-1)&&(i<MaxBoolContainerDim),"Index out of range");
                ASSERTL1((j>-1)&&(j<MaxBoolContainerDim),"Index out of range");

                m_data[i][j] = a;
            }

            bool GetValue(int i, int j) const
            {
                ASSERTL1((i>-1)&&(i<MaxBoolContainerDim),"Index out of range");
                ASSERTL1((j>-1)&&(j<MaxBoolContainerDim),"Index out of range");

                return m_data[i][j];
            }

        private:
            bool m_data[MaxBoolContainerDim][MaxBoolContainerDim];
        };

        // Implementation for 3D expansions
        template<>
        class BoolContainer<3>
        {
        public:
            BoolContainer()
            {
                std::fill_n(&m_data[0][0][0],MaxBoolContainerDim*MaxBoolContainerDim*MaxBoolContainerDim,false);
            }

            void SetValue(int i, int j, int k, bool a)
            {
                ASSERTL1((i>-1)&&(i<MaxBoolContainerDim),"Index out of range");
                ASSERTL1((j>-1)&&(j<MaxBoolContainerDim),"Index out of range");
                ASSERTL1((k>-1)&&(k<MaxBoolContainerDim),"Index out of range");

                m_data[i][j][k] = a;
            }

            bool GetValue(int i, int j, int k) const
            {
                ASSERTL1((i>-1)&&(i<MaxBoolContainerDim),"Index out of range");
                ASSERTL1((j>-1)&&(j<MaxBoolContainerDim),"Index out of range");
                ASSERTL1((k>-1)&&(k<MaxBoolContainerDim),"Index out of range");

                return m_data[i][j][k];
            }

        private:
            bool m_data[MaxBoolContainerDim][MaxBoolContainerDim][MaxBoolContainerDim];
        };


        //default implementation (e.g. for 1D and 3D which are not implemented yet)
        template<>
        class BoolContainer<0>
        {
        public:
            BoolContainer():
                m_data(false)
            {
            }

            void SetValue(bool a)
            {
                m_data = a;
            }

            bool GetValue() const
            {
                return m_data;
            }

        private:
            bool m_data;
        };
        // ===============================================================================



        // ===============================================================================
        // Traits of the ElementalOptimizationParameters class
        // ============

        // Default implementation
        template<StdRegions::ElementType etype>
        class ElementalOptimizationParametersTraits
        {
        public:
            typedef BoolContainer<0> DataHolderType;
            static int const elementDim = 0;
        };

        // template specializations
        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eStdQuadExp>
        {
        public:
            static int const elementDim = 2;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eStdTriExp>
        {
        public:
            static int const elementDim = 2;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eStdNodalTriExp>
        {
        public:
            static int const elementDim = 2;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eQuadExp>
        {
        public:
            static int const elementDim = 2;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eTriExp>
        {
        public:
            static int const elementDim = 2;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eNodalTriExp>
        {
        public:
            static int const elementDim = 2;
            typedef BoolContainer<elementDim> DataHolderType;
        };


        // ----------------------
        // 3D elements - NEEDS FIXING
        // ----------------------
        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eStdHexExp>
        {
        public:
            static int const elementDim = 3;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eHexExp>
        {
        public:
            static int const elementDim = 3;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eStdTetExp>
        {
        public:
            static int const elementDim = 3;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        template<>
        class ElementalOptimizationParametersTraits<StdRegions::eTetExp>
        {
        public:
            static int const elementDim = 3;
            typedef BoolContainer<elementDim> DataHolderType;
        };

        // ===============================================================================

        // ===============================================================================
        // Dump policy of the ElementalOptimizationParameters class
        // ============

        // Default implementation
        template<int dim>
        class ElementalOptimizationParametersDumpPolicy;

        template<>
        class ElementalOptimizationParametersDumpPolicy<2>
        {
        public:
            static void DumpElementalOptimizationParameters(std::ostream &outfile, BoolContainer<2> data)
            {
                int i,j;
                outfile << std::setw(2) << "  " << std::setw(1) << "|";
                for(i = 0; i < MaxBoolContainerDim; i++)
                {
                    outfile << std::right << std::setw(3) << i+2 << std::setw(1) << " ";
                }
                outfile << std::endl;
                outfile << "---";
                for(i = 0; i < MaxBoolContainerDim; i++)
                {
                    outfile << right << std::setw(4) << "----";
                }
                outfile << endl;

                for(i = 0; i < MaxBoolContainerDim; i++)
                {
                    outfile << std::setw(2) << i+2 << std::setw(1) << "|";
                    for( j = 0; j < MaxBoolContainerDim; j++)
                    {
                        outfile << std::right << std::setw(3)  << data.GetValue(i,j) << std::setw(1) <<" ";
                    }
                    outfile << std::endl;
                }
            }
        };

        template<>
        class ElementalOptimizationParametersDumpPolicy<3>
        {
        public:
            static void DumpElementalOptimizationParameters(std::ostream &outfile, BoolContainer<3> data)
            {
                int k = 0;
                for (k = 0; k < MaxBoolContainerDim; ++k)
                {
                    int i,j;
                    outfile << std::setw(2) << "  " << std::setw(1) << "|";
                    for(i = 0; i < MaxBoolContainerDim; i++)
                    {
                        outfile << std::right << std::setw(3) << i+2 << std::setw(1) << " ";
                    }
                    outfile << std::endl;
                    outfile << "---";
                    for(i = 0; i < MaxBoolContainerDim; i++)
                    {
                        outfile << right << std::setw(4) << "----";
                    }
                    outfile << endl;

                    for(i = 0; i < MaxBoolContainerDim; i++)
                    {
                        outfile << std::setw(2) << i+2 << std::setw(1) << "|";
                        for( j = 0; j < MaxBoolContainerDim; j++)
                        {
                            outfile << std::right << std::setw(3)  << data.GetValue(i,j,k) << std::setw(1) <<" ";
                        }
                        outfile << std::endl;
                    }
                }
            }
        };

        // ===============================================================================

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_STDREGIONS_ELEMENTALOPTIMIZATIONPARAMETERSTRAITS_HPP


