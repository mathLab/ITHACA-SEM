///////////////////////////////////////////////////////////////////////////////
//
// File: BandedMatrixStoragePolicy.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BANDED_MATRIX_STORAGE_POLICY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BANDED_MATRIX_STORAGE_POLICY_HPP

#include <LibUtilities/LinearAlgebra/MatrixStoragePolicy.hpp>
#include <boost/call_traits.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/optional/optional.hpp>
#include <limits>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace Nektar
{
    template<>
    class MatrixStoragePolicy<BandedMatrixTag>
    {
        public:                        
            /// \brief Calculates and returns the storage size required.
            ///
            /// This method assumes that the matrix will be used with LU factorizationa and 
            /// allocates additional storage as appropriate.
            static unsigned int GetRequiredStorageSize(unsigned int totalRows, unsigned int totalColumns,
                                                       unsigned int subDiags, unsigned int superDiags)
            {
                return CalculateNumberOfRows(totalRows, subDiags, superDiags)*totalColumns;
            }

            static unsigned int CalculateNumberOfDiags(unsigned int totalRows, unsigned int diags)
            {
                if( diags != std::numeric_limits<unsigned int>::max() )
                {
                    return diags;
                }
                else if( totalRows > 0 )
                {
                    return totalRows-1;
                }
                else
                {
                    return 0;
                }
            }
                
            static unsigned int CalculateNumberOfRows(unsigned int totalRows, unsigned int subDiags, unsigned int superDiags)
            {
                return CalculateNumberOfDiags(totalRows, subDiags) + CalculateNumberOfDiags(totalRows, superDiags) + 1;
            }

            static boost::optional<unsigned int> CalculateIndex(unsigned int totalRows, 
                                                                unsigned int totalColumns,
                                                                unsigned int row, unsigned int column,
                                                                unsigned int sub, unsigned int super)
            {
                if( (column <= row && (row - column) <= CalculateNumberOfDiags(totalRows, sub)) ||
                    (column > row && (column - row) <= CalculateNumberOfDiags(totalRows, super)) )
                {
                    unsigned int elementRow = CalculateNumberOfDiags(totalRows, super)+row-column;
                    unsigned int elementColumn = column;

                    return elementRow + elementColumn*CalculateNumberOfRows(totalRows, sub, super);
                }
                else
                {
                    return boost::optional<unsigned int>();
                }
            }


            static boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(const unsigned int totalRows, const unsigned int totalColumns,
                    const unsigned int curRow, const unsigned int curColumn)
            {
                unsigned int nextRow = curRow;
                unsigned int nextColumn = curColumn;

                //if( transpose == 'N' )
                //{
                //    if( (column <= row && (row - column) <= data.GetNumberOfSubDiagonals(totalRows)) ||
                //    (column > row && (column - row) <= data.GetNumberOfSuperDiagonals(totalRows)) )
                //    {
                //        unsigned int arrayColumns = totalColumns;

                //        unsigned int elementRow = data.GetNumberOfSuperDiagonals(totalRows)+row-column;
                //        unsigned int elementColumn = column;

                //        return elementRow + elementColumn*CalculateNumberOfRows(totalRows, data);
                //    }
                //    else
                //    {
                //        return boost::optional<unsigned int>();
                //    }
                //}
                //else
                //{
                //}
                return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BANDED_MATRIX_STORAGE_POLICY_HPP
