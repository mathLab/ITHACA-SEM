///////////////////////////////////////////////////////////////////////////////
//
// File: DiagonalMatrixStoragePolicy.hpp
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
// Description: Interface classes for matrices
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DIAGONAL_MATRIX_STORAGE_POLICY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DIAGONAL_MATRIX_STORAGE_POLICY_HPP

#include <LibUtilities/LinearAlgebra/MatrixStoragePolicy.hpp>
#include <boost/call_traits.hpp>
#include "boost/tuple/tuple.hpp"

namespace Nektar
{
    template<>
    class MatrixStoragePolicy<DiagonalMatrixTag>
    {
        public:
            
            static boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(const unsigned int totalRows, const unsigned int totalColumns,
                    const unsigned int curRow, const unsigned int curColumn)
            {
                ASSERTL0(curRow == curColumn, "Iteration of a diagonal matrix is only valid along the diagonal.");

                unsigned int nextRow = curRow;
                unsigned int nextColumn = curColumn;

                if( nextRow < totalRows )
                {
                    ++nextRow;
                    ++nextColumn;
                }

                if( nextRow >= totalRows )
                {
                    nextRow = std::numeric_limits<unsigned int>::max();
                    nextColumn = std::numeric_limits<unsigned int>::max();
                }

                return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
            }

            template<typename DataType>
            static void Invert(unsigned int rows, unsigned int columns,
                               Array<OneD, DataType>& data,
                               const char transpose)
            {
                ASSERTL0(rows==columns, "Only square matrices can be inverted.");
                for(unsigned int i = 0; i < rows; ++i)
                {
                    data[i] = 1.0/data[i];
                }
            }
            
            static unsigned int GetRequiredStorageSize(unsigned int rows, unsigned int columns)
            {
                ASSERTL0(rows==columns, "Diagonal matrices must be square.");
                return rows;
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DIAGONAL_MATRIX_STORAGE_POLICY_HPP
