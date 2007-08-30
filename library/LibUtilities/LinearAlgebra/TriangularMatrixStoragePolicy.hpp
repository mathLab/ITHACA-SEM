///////////////////////////////////////////////////////////////////////////////
//
// File: TriangularMatrixStoragePolicy.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_TRIANGULAR_MATRIX_STORAGE_POLICY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_TRIANGULAR_MATRIX_STORAGE_POLICY_HPP

#include <LibUtilities/LinearAlgebra/MatrixStoragePolicy.hpp>
#include <boost/call_traits.hpp>
#include "boost/tuple/tuple.hpp"

namespace Nektar
{
    template<typename DataType>
    class TriangularMatrixStoragePolicy
    {
        public:
            typedef typename boost::call_traits<DataType>::value_type GetValueReturnType;
            typedef DefaultPolicySpecificDataHolder PolicySpecificDataHolderType;
            static DataType ZeroElement;
            
            static Array<OneD, DataType> Initialize()
            {
                return Array<OneD, DataType>();
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, const PolicySpecificDataHolderType&)
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                return Array<OneD, DataType>(rows*(rows+1)/2);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    typename boost::call_traits<DataType>::const_reference d,
                                                    const PolicySpecificDataHolderType& )
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                return Array<OneD, DataType>(rows*(rows+1)/2, d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const DataType* d, const PolicySpecificDataHolderType&)
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                return Array<OneD, DataType>(rows*(rows+1)/2, d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const ConstArray<OneD, DataType>& d, const PolicySpecificDataHolderType&)
            {
                unsigned int size = rows*(rows+1)/2;

                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                ASSERTL0(size <= d.num_elements(), 
                    std::string("An attempt has been made to create a triangular matrix of size ") +
                    boost::lexical_cast<std::string>(rows*(rows+1)/2) + 
                    std::string(" but the array being used to populate it only has ") + 
                    boost::lexical_cast<std::string>(d.num_elements()) + 
                    std::string(" elements."));
                
                Array<OneD, DataType> result(size, d);
                return result;
            }
    };
    
    template<typename DataType>
    DataType TriangularMatrixStoragePolicy<DataType>::ZeroElement = DataType(0);

    template<typename DataType>
    class MatrixStoragePolicy<DataType, UpperTriangularMatrixTag> : public TriangularMatrixStoragePolicy<DataType>
    {
        public:
            typedef TriangularMatrixStoragePolicy<DataType> BaseType;
            typedef typename BaseType::GetValueReturnType GetValueReturnType;
            typedef typename BaseType::PolicySpecificDataHolderType PolicySpecificDataHolderType;

            static typename boost::call_traits<DataType>::const_reference GetValue(unsigned int totalRows, unsigned int totalColumns,
                                                                             unsigned int curRow, unsigned int curColumn,
                                                                             const Array<OneD, DataType>& data,
                                                                             const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to retrieve a value from row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to retrieve a value from column " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");

                if( curRow <= curColumn )
                {
                    return data[CalculateIndex(totalRows, curRow, curColumn)];
                }
                else
                {
                    return BaseType::ZeroElement;
                }
            }            
            
            static void SetValue(unsigned int totalRows, unsigned int totalColumns,
                                 unsigned int curRow, unsigned int curColumn,
                                 Array<OneD, DataType>& data, typename boost::call_traits<DataType>::const_reference d,
                                 const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to set a value from row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to set a value from column " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curRow <= curColumn, "Attemping to set element (" +
                    boost::lexical_cast<std::string>(curRow) + ", " +
                    boost::lexical_cast<std::string>(curColumn) + ") of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");

                data[CalculateIndex(totalRows, curRow, curColumn)] = d;
            }
                    
            static unsigned int CalculateIndex(unsigned int totalRows, unsigned int curRow, unsigned int curColumn)
            {
                return curColumn + (2*totalRows - curRow - 1)*(curRow)/2;
            }

            static boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(const unsigned int totalRows, const unsigned int totalColumns,
                    const unsigned int curRow, const unsigned int curColumn, 
                    const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to iterate through an element on row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to iterate through an element on row " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curRow <= curColumn, "Attemping to iterate through element (" +
                    boost::lexical_cast<std::string>(curRow) + ", " +
                    boost::lexical_cast<std::string>(curColumn) + ") of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");

                unsigned int nextRow = curRow;
                unsigned int nextColumn = curColumn;

                if( nextColumn < totalColumns )
                {
                    ++nextColumn;
                }

                if( nextColumn >= totalColumns )
                {
                    ++nextRow;
                    nextColumn = nextRow;
                }
                
                if( nextRow >= totalRows )
                {
                    nextRow = std::numeric_limits<unsigned int>::max();
                    nextColumn = std::numeric_limits<unsigned int>::max();
                }

                return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
            }
    };
    
    
    
    
    template<typename DataType>
    class MatrixStoragePolicy<DataType, LowerTriangularMatrixTag> : public TriangularMatrixStoragePolicy<DataType>
    {
        public:
            typedef TriangularMatrixStoragePolicy<DataType> BaseType;
            typedef typename BaseType::GetValueReturnType GetValueReturnType;
            typedef typename BaseType::PolicySpecificDataHolderType PolicySpecificDataHolderType;

            static typename boost::call_traits<DataType>::const_reference GetValue(unsigned int totalRows, unsigned int totalColumns,
                                                                             unsigned int curRow, unsigned int curColumn,
                                                                             const Array<OneD, DataType>& data,
                                                                             const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to retrieve a value from row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to retrieve a value from column " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");

                if( curRow >= curColumn )
                {
                    return data[CalculateIndex(totalRows, curRow, curColumn)];
                }
                else
                {
                    return BaseType::ZeroElement;
                }
            }            
            
            static void SetValue(unsigned int totalRows, unsigned int totalColumns,
                                 unsigned int curRow, unsigned int curColumn,
                                 Array<OneD, DataType>& data, typename boost::call_traits<DataType>::const_reference d,
                                 const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to set a value from row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to set a value from column " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
                ASSERTL1(curRow >= curColumn, "Attemping to set element (" +
                    boost::lexical_cast<std::string>(curRow) + ", " +
                    boost::lexical_cast<std::string>(curColumn) + ") of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");

                data[CalculateIndex(totalRows, curRow, curColumn)] = d;
            }
                    
            static unsigned int CalculateIndex(unsigned int totalRows, unsigned int curRow, unsigned int curColumn)
            {
                return curColumn + curRow*(curRow+1)/2;
            }

            static boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(const unsigned int totalRows, const unsigned int totalColumns,
                    const unsigned int curRow, const unsigned int curColumn,
                    const PolicySpecificDataHolderType&)
            {
                ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
                ASSERTL1(curRow < totalRows, "Attemping to iterate through an element on row " +
                    boost::lexical_cast<std::string>(curRow) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");
                ASSERTL1(curColumn < totalColumns, "Attemping to iterate through an element on row " +
                    boost::lexical_cast<std::string>(curColumn) + " of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");
                ASSERTL1(curRow >= curColumn, "Attemping to iterate through element (" +
                    boost::lexical_cast<std::string>(curRow) + ", " +
                    boost::lexical_cast<std::string>(curColumn) + ") of a (" +
                    boost::lexical_cast<std::string>(totalRows) + ", " +
                    boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");

                unsigned int nextRow = curRow;
                unsigned int nextColumn = curColumn;

                if( nextColumn < curRow )
                {
                    ++nextColumn;
                }

                if( nextColumn >= curRow )
                {
                    ++nextRow;
                    nextColumn = 0;
                }

                if( nextRow >= totalRows )
                {
                    nextRow = std::numeric_limits<unsigned int>::max();
                    nextColumn = std::numeric_limits<unsigned int>::max();
                }

                return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_TRIANGULAR_MATRIX_STORAGE_POLICY_HPP
