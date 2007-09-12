///////////////////////////////////////////////////////////////////////////////
//
// File: SymmetricMatrixStoragePolicy.hpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SYMMETRIC_MATRIX_STORAGE_POLICY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SYMMETRIC_MATRIX_STORAGE_POLICY_HPP

#include <LibUtilities/LinearAlgebra/MatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/TriangularMatrixStoragePolicy.hpp>
#include <boost/call_traits.hpp>
#include <boost/tuple/tuple.hpp>

namespace Nektar
{
    template<typename DataType>
    class MatrixStoragePolicy<DataType, SymmetricMatrixTag> : private TriangularMatrixStoragePolicy<DataType>
    {
        public:
            typedef TriangularMatrixStoragePolicy<DataType> BaseType;
            typedef typename BaseType::GetValueReturnType GetValueReturnType;
            typedef typename TriangularMatrixStoragePolicy<DataType>::PolicySpecificDataHolderType PolicySpecificDataHolderType;
            
            using BaseType::Initialize;
            static typename boost::call_traits<DataType>::const_reference GetValue(unsigned int totalRows, unsigned int totalColumns,
                                                                                   unsigned int curRow, unsigned int curColumn,
                                                                                   const ConstArray<OneD, DataType>& data,
                                                                                   const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Symmetric matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to get a value from row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to get a value from column " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");

                if( curRow <= curColumn )
                {
                    return data[CalculateIndex(curRow, curColumn)];
                }
                else
                {
                    return data[CalculateIndex(curColumn, curRow)];
                }
            }            
            
            static void SetValue(unsigned int totalRows, unsigned int totalColumns,
                                 unsigned int curRow, unsigned int curColumn,
                                 Array<OneD, DataType>& data, typename boost::call_traits<DataType>::const_reference d,
                                 const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Symmetric matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to set a value from row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to set a value from column " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");

                if( curRow <= curColumn )
                {
                    data[CalculateIndex(curRow, curColumn)] = d;
                }
                else
                {
                    data[CalculateIndex(curColumn, curRow)] = d;
                }   
            }
                    
            static unsigned int CalculateIndex(unsigned int curRow, unsigned int curColumn)
            {
                return curRow + curColumn*(curColumn+1)/2;
            }

            static boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(const unsigned int totalRows, const unsigned int totalColumns,
                    const unsigned int curRow, const unsigned int curColumn,
                    const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Symmetric matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to iterate through an element on row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to iterate through an element on row " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");

                unsigned int nextRow = curRow;
                unsigned int nextColumn = curColumn;

                if( nextRow < totalRows )
                {
                    ++nextRow;
                }

                if( nextRow >= totalRows )
                {
                    nextRow = 0;
                    ++nextColumn;
                }

                if( nextColumn >= totalColumns )
                {
                    nextRow = std::numeric_limits<unsigned int>::max();
                    nextColumn = std::numeric_limits<unsigned int>::max();
                }

                return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SYMMETRIC_MATRIX_STORAGE_POLICY_HPP
