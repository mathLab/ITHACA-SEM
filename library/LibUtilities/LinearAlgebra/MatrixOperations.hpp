///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixOperations.hpp
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
// Description: Defines the global functions needed for matrix operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

// Since this file defines all of the operations for all combination of matrix types, 
// we have to include all matrix specializations first.
#include <LibUtilities/LinearAlgebra/NormalMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>

namespace Nektar
{
    template<typename DataType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& lhs,
                const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
//         ASSERTL2(result.GetRows() == lhs.GetRows() && result.GetColumns() == rhs.GetColumns(),
//             std::string("The result matrix passed to NekAdd has the wrong number of rows and columns.  It has (") +
//             boost::lexical_cast<std::string>(result.GetRows()) + std::string(",") + 
//             boost::lexical_cast<std::string>(result.GetColumns()) + std::string("), but ") + 
//             boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(",") + 
//             boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string("), was expected.")); 
            
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs);
        typename NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>::iterator result_iter = result.begin();
        typename NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>::const_iterator rhs_iter = rhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) += (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
        }
    }
    
    template<typename DataType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& lhs,
                        const NekMatrix<DataType, DiagonalMatrixTag, StandardMatrixTag>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
        
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) += rhs(i,i);
        }
    }
    
    template<typename DataType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<DataType, DiagonalMatrixTag, StandardMatrixTag>& lhs,
                        const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
        
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(rhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) += lhs(i,i);
        }
    }
    
    
    template<typename DataType>
    void NekMultiply(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                            const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& lhs,
                            const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        ASSERTL0(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));
           
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs.GetRows(), rhs.GetColumns());

        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < result.GetColumns(); ++j)
            {
                DataType t = DataType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < lhs.GetColumns(); ++k)
                {
                    t += lhs(i,k)*rhs(k,j);
                }
                result(i,j) = t;
            }
        }
    }
    
    #ifdef NEKTAR_USING_BLAS    
    /// \brief Floating point specialization when blas is in use.
    ///
    /// The signature of this function should match the unspecialized NekMultiply.
    void NekMultiply(NekMatrix<NekDouble, FullMatrixTag, StandardMatrixTag>& result,
                            const NekMatrix<NekDouble, FullMatrixTag, StandardMatrixTag>& lhs,
                            const NekMatrix<NekDouble, FullMatrixTag, StandardMatrixTag>& rhs);
    #endif //NEKTAR_USING_BLAS
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

