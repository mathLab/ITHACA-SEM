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
    template<typename DataType>
    class MatrixStoragePolicy<DataType, BandedMatrixTag>
    {
        public:
            class PolicySpecificDataHolderType
            {
                public:
                    PolicySpecificDataHolderType() :
                        m_numberOfSubDiagonals(std::numeric_limits<unsigned int>::max()),
                        m_numberOfSuperDiagonals(std::numeric_limits<unsigned int>::max())
                    {
                    }

                    PolicySpecificDataHolderType(unsigned int sub, unsigned int super) :
                        m_numberOfSubDiagonals(sub),
                        m_numberOfSuperDiagonals(super)
                    {
                    }

                    // Get the specified number of sub diagonals, or calculate the 
                    // number if this object has not been initialized.
                    unsigned int GetNumberOfSubDiagonals(unsigned int rows) const
                    {
                        if( m_numberOfSubDiagonals != std::numeric_limits<unsigned int>::max() )
                        {
                            return m_numberOfSubDiagonals;
                        }
                        else if( rows > 0 )
                        {
                            return rows-1;
                        }
                        else
                        {
                            return 0;
                        }
                    }

                    unsigned int GetNumberOfSuperDiagonals(unsigned int rows) const
                    {
                        if( m_numberOfSuperDiagonals != std::numeric_limits<unsigned int>::max() )
                        {
                            return m_numberOfSuperDiagonals;
                        }
                        else if( rows > 0 )
                        {
                            return rows-1;
                        }
                        else
                        {
                            return 0;
                        }
                    }

                private:
                    unsigned int m_numberOfSubDiagonals;
                    unsigned int m_numberOfSuperDiagonals;
            };

            typedef typename boost::call_traits<DataType>::value_type GetValueReturnType;
            static DataType ZeroElement;
            
            static Array<OneD, DataType> Initialize()
            {
                return Array<OneD, DataType>();
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, const PolicySpecificDataHolderType& data)
            {
                ASSERTL0(rows==columns, "Banded matrices must be square.");
                return Array<OneD, DataType>(CalculateStorageSize(rows, columns, data));
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    typename boost::call_traits<DataType>::const_reference d,
                                                    const PolicySpecificDataHolderType& data)
            {
                ASSERTL0(rows==columns, "Banded matrices must be square.");
                return Array<OneD, DataType>(CalculateStorageSize(rows, columns, data), d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const DataType* d, const PolicySpecificDataHolderType& data)
            {
                ASSERTL0(rows==columns, "Banded matrices must be square.");
                return Array<OneD, DataType>(CalculateStorageSize(rows, columns, data), d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const ConstArray<OneD, DataType>& d,
                                                    const PolicySpecificDataHolderType& data)
            {
                ASSERTL0(rows==columns, "Banded matrices must be square.");

                unsigned int storageSize = CalculateStorageSize(rows, columns, data);
                ASSERTL0(storageSize > d.num_elements(), 
                    std::string("An attempt has been made to create a banded matrix of size (") +
                    boost::lexical_cast<std::string>(rows) + 
                    ", " + boost::lexical_cast<std::string>(columns) + ") " +
                    std::string(" but the array being used to populate it only has ") + 
                    boost::lexical_cast<std::string>(d.num_elements()) + 
                    std::string(" elements."));
                return Array<OneD, DataType>(storageSize, d.data());
            }
            
            
            /// \brief Calculates and returns the storage size required.
            ///
            /// This method assumes that the matrix will be used with LU factorizationa and 
            /// allocates additional storage as appropriate.
            static unsigned int CalculateStorageSize(unsigned int totalRows, unsigned int totalColumns,
                                                     const PolicySpecificDataHolderType& data)
            {
                return CalculateNumberOfRows(totalRows, data)*totalColumns;
            }

            static unsigned int CalculateNumberOfRows(unsigned int totalRows, const PolicySpecificDataHolderType& data)
            {
                return data.GetNumberOfSubDiagonals(totalRows) + data.GetNumberOfSuperDiagonals(totalRows) + 1;
            }

            static boost::optional<unsigned int> CalculateIndex(unsigned int totalRows, 
                                                                unsigned int totalColumns,
                                                                unsigned int row, unsigned int column,
                                                                const PolicySpecificDataHolderType& data)
            {
                if( (column <= row && (row - column) >= data.GetNumberOfSubDiagonals(totalRows)) ||
                    (column > row && (column - row) >= data.GetNumberOfSuperDiagonals(totalRows)) )
                {
                    unsigned int arrayColumns = totalColumns;

                    unsigned int elementRow = data.GetNumberOfSuperDiagonals(totalRows)+row-column;
                    unsigned int elementColumn = column;

                    return elementRow + elementColumn*CalculateNumberOfRows(totalRows, data);
                }
                else
                {
                    return boost::optional<unsigned int>();
                }
            }

            //static GetValueReturnType GetValue(unsigned int totalRows, unsigned int totalColumns,
            //                                   unsigned int curRow, unsigned int curColumn,
            //                                   const ConstArray<OneD, DataType>& d,
            //                                   const PolicySpecificDataHolderType& data)
            //{
            //    boost::optional<unsigned int> index = CalculateIndex(totalRows, totalColumns, curRow, curColumn, data);
            //    if( index )
            //    {
            //        return data[*index];
            //    }
            //    else
            //    {
            //        return ZeroElement;
            //    }
            //}

            static typename boost::call_traits<DataType>::const_reference GetValue(unsigned int totalRows, unsigned int totalColumns,
                                                                             unsigned int curRow, unsigned int curColumn,
                                                                             const ConstArray<OneD, DataType>& data,
                                                                             const PolicySpecificDataHolderType& dataHolder)
            {
                boost::optional<unsigned int> index = CalculateIndex(totalRows, totalColumns, curRow, curColumn, dataHolder);
                if( index )
                {
                    return data[*index];
                }
                else
                {
                    return ZeroElement;
                }
            }            
            
            static void SetValue(unsigned int totalRows, unsigned int totalColumns,
                                 unsigned int curRow, unsigned int curColumn,
                                 Array<OneD, DataType>& data, typename boost::call_traits<DataType>::const_reference d,
                                 const PolicySpecificDataHolderType& dataHolder)
            {
                boost::optional<unsigned int> index = CalculateIndex(totalRows, totalColumns, curRow, curColumn, dataHolder);
                if( index )
                {
                    data[curRow] = d;
                }
                else
                {
                    NEKERROR(ErrorUtil::efatal, "Can only assign into banded portion of block matrix.");
                }
            }
            
            //static boost::tuples::tuple<unsigned int, unsigned int> 
            //Advance(const unsigned int totalRows, const unsigned int totalColumns,
            //        const unsigned int curRow, const unsigned int curColumn,
            //        const PolicySpecificDataHolderType&)
            //{
            //    ASSERTL0(curRow == curColumn, "Iteration of a diagonal matrix is only valid along the diagonal.");

            //    unsigned int nextRow = curRow;
            //    unsigned int nextColumn = curColumn;

            //    if( nextRow < totalRows )
            //    {
            //        ++nextRow;
            //        ++nextColumn;
            //    }

            //    if( nextRow >= totalRows )
            //    {
            //        nextRow = std::numeric_limits<unsigned int>::max();
            //        nextColumn = std::numeric_limits<unsigned int>::max();
            //    }

            //    return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
            //}

            //static void Invert(unsigned int rows, unsigned int columns,
            //                   Array<OneD, DataType>& data,
            //                   const PolicySpecificDataHolderType&)
            //{
            //    ASSERTL0(rows==columns, "Only square matrices can be inverted.");
            //    for(unsigned int i = 0; i < rows; ++i)
            //    {
            //        data[i] = 1.0/data[i];
            //    }
            //}
    };
    
    template<typename DataType>
    DataType MatrixStoragePolicy<DataType, BandedMatrixTag>::ZeroElement = DataType(0);
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BANDED_MATRIX_STORAGE_POLICY_HPP
