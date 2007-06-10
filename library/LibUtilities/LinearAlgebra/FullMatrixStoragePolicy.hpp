///////////////////////////////////////////////////////////////////////////////
//
// File: FullMatrixStoragePolicy.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_FULL_MATRIX_STORAGE_POLICY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_FULL_MATRIX_STORAGE_POLICY_HPP

#include <LibUtilities/LinearAlgebra/MatrixStoragePolicy.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    template<typename DataType>
    class MatrixStoragePolicy<DataType, FullMatrixTag>
    {
        public:
            typedef typename boost::call_traits<DataType>::reference GetValueReturnType;
        
            static Array<OneD, DataType> Initialize()
            {
                return Array<OneD, DataType>();
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns)
            {
                return Array<OneD, DataType>(rows*columns);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    typename boost::call_traits<DataType>::const_reference d)
            {
                return Array<OneD, DataType>(rows*columns, d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const DataType* d)
            {
                return Array<OneD, DataType>(rows*columns, d);
            }
            
            static Array<OneD, DataType> Initialize(unsigned int rows, unsigned int columns, 
                                                    const ConstArray<OneD, DataType>& d)
            {
                ASSERTL0(rows*columns <= d.num_elements(), 
                    std::string("An attempt has been made to create a full matrix of size ") +
                    boost::lexical_cast<std::string>(rows*columns) + 
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
                return data[curRow*totalColumns + curColumn];
            }
            
            static typename boost::call_traits<DataType>::const_reference GetValue(unsigned int totalRows, unsigned int totalColumns,
                                                                             unsigned int curRow, unsigned int curColumn,
                                                                             const Array<OneD, DataType>& data)
            {
                return data[curRow*totalColumns + curColumn];
            }
            
            static void SetValue(unsigned int totalRows, unsigned int totalColumns,
                                 unsigned int curRow, unsigned int curColumn,
                                 Array<OneD, DataType>& data, typename boost::call_traits<DataType>::const_reference d)
            {
                data[curRow*totalColumns + curColumn] = d;
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_FULL_MATRIX_STORAGE_POLICY_HPP
