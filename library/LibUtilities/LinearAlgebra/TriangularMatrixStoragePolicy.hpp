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

namespace Nektar
{
    template<typename DataType>
    class MatrixStoragePolicy<DataType, UpperTriangularMatrixTag>
    {
        public:
            typedef typename boost::call_traits<DataType>::value_type GetValueReturnType;
            static DataType ZeroElement;
            
            static Array<OneD, DataType> Initialize()
            {
                return Array<OneD, DataType>();
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns)
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                return Array<OneD, DataType>(rows*(rows+1)/2);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    typename boost::call_traits<DataType>::const_reference d)
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                return Array<OneD, DataType>(rows*(rows+1)/2, d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const DataType* d)
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                return Array<OneD, DataType>(rows*(rows+1)/2, d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const ConstArray<OneD, DataType>& d)
            {
                ASSERTL0(rows==columns, "Triangular matrices must be square.");
                ASSERTL0(rows*(rows+1)/2 <= d.num_elements(), 
                    std::string("An attempt has been made to create a triangular matrix of size ") +
                    boost::lexical_cast<std::string>(rows*(rows+1)/2) + 
                    std::string(" but the array being used to populate it only has ") + 
                    boost::lexical_cast<std::string>(d.num_elements()) + 
                    std::string(" elements."));
                Array<OneD, DataType> result;
                CopyArray(d, result);
                return result;
            }
            
            static GetValueReturnType GetValue(unsigned int totalRows, unsigned int totalColumns,
                                               unsigned int curRow, unsigned int curColumn,
                                               Array<OneD, DataType>& data)
            {
                if( curRow <= curColumn )
                {
                    return data[CalculateIndex(totalRows, curRow, curColumn)];
                }
                else
                {
                    return ZeroElement;
                }
            }
            
            static typename boost::call_traits<DataType>::const_reference GetValue(unsigned int totalRows, unsigned int totalColumns,
                                                                             unsigned int curRow, unsigned int curColumn,
                                                                             const Array<OneD, DataType>& data)
            {
                if( curRow <= curColumn )
                {
                    return data[CalculateIndex(totalRows, curRow, curColumn)];
                }
                else
                {
                    return ZeroElement;
                }
            }            
            
            static void SetValue(unsigned int totalRows, unsigned int totalColumns,
                                 unsigned int curRow, unsigned int curColumn,
                                 Array<OneD, DataType>& data, typename boost::call_traits<DataType>::const_reference d)
            {
                ASSERTL0(curRow <= curColumn, "Can only assign into the upper triangular portion of an upper triangular matrix.");
                data[CalculateIndex(totalRows, curRow, curColumn)] = d;
            }
                    
            static unsigned int CalculateIndex(unsigned int totalRows, unsigned int curRow, unsigned int curColumn)
            {
                unsigned int base = curRow*(-static_cast<int>(curRow) + 1 + 2*totalRows)/2;
                return static_cast<unsigned int>(base + curColumn - curRow);
            }
    };
    
    template<typename DataType>
    DataType MatrixStoragePolicy<DataType, UpperTriangularMatrixTag>::ZeroElement = DataType(0);
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_TRIANGULAR_MATRIX_STORAGE_POLICY_HPP
