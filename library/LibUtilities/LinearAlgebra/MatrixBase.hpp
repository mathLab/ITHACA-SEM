///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixBase.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP

#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/LinearAlgebra/MatrixType.h>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#ifdef max
#undef max
#endif

namespace Nektar
{    
    template<typename DataType>
    class ConstMatrix
    {
        public:
            LIB_UTILITIES_EXPORT virtual ~ConstMatrix();
            
        public:
            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::value_type operator()(unsigned int row, unsigned int column) const;
            
            LIB_UTILITIES_EXPORT unsigned int GetStorageSize() const;
            
            LIB_UTILITIES_EXPORT inline MatrixStorage GetType() const
            {
                return m_storageType;
            }

            LIB_UTILITIES_EXPORT inline MatrixStorage GetStorageType() const
            {
                return m_storageType;
            }

            LIB_UTILITIES_EXPORT unsigned int GetRows() const;

            LIB_UTILITIES_EXPORT unsigned int GetTransposedRows(char transpose) const ;
            
            LIB_UTILITIES_EXPORT unsigned int GetColumns() const;

            LIB_UTILITIES_EXPORT unsigned int GetTransposedColumns(char transpose) const ;

            LIB_UTILITIES_EXPORT const unsigned int* GetSize() const;
            
            LIB_UTILITIES_EXPORT void Transpose() ;

            LIB_UTILITIES_EXPORT char GetTransposeFlag() const;
            
            LIB_UTILITIES_EXPORT static unsigned int CalculateIndex(MatrixStorage type, 
                unsigned int row, unsigned int col, 
                unsigned int numRows, unsigned int numColumns, const char transpose =  'N',
                unsigned int numSubDiags = 0, unsigned int numSuperDiags = 0) ;
            
            LIB_UTILITIES_EXPORT static unsigned int GetRequiredStorageSize(MatrixStorage type, unsigned int rows, 
                unsigned int columns, unsigned int subDiags = 0, unsigned int superDiags = 0);

        protected:
            
            // All constructors are private to enforce the abstract nature of ConstMatrix without
            // resorting to pure virtual functions.
            LIB_UTILITIES_EXPORT ConstMatrix(unsigned int rows, unsigned int columns,
                                             MatrixStorage policy = eFULL);
            
            LIB_UTILITIES_EXPORT ConstMatrix(const ConstMatrix<DataType>& rhs);

            LIB_UTILITIES_EXPORT ConstMatrix<DataType>& operator=(const ConstMatrix<DataType>& rhs);
            
            /// \brief Resets the rows and columns in the array.
            /// This method does not update the data storage to match the new row and column counts.
            LIB_UTILITIES_EXPORT void Resize(unsigned int rows, unsigned int columns);

            LIB_UTILITIES_EXPORT void SetTransposeFlag(char newValue);
            LIB_UTILITIES_EXPORT inline char GetRawTransposeFlag() const
            {
                return m_transpose;
            }

        private:
            
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const = 0;            
            virtual unsigned int v_GetStorageSize() const = 0;            
            LIB_UTILITIES_EXPORT virtual void v_Transpose();
            LIB_UTILITIES_EXPORT virtual char v_GetTransposeFlag() const;
            unsigned int m_size[2];
            char m_transpose;
            MatrixStorage m_storageType;
    };
    
    template<typename DataType>
    class Matrix : public ConstMatrix<DataType>
    {  
        public:
            LIB_UTILITIES_EXPORT virtual ~Matrix();
            
            LIB_UTILITIES_EXPORT void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d);
            
        protected:            
            // All constructors are private to enforce the abstract nature of ConstMatrix without
            // resorting to pure virtual functions.
            LIB_UTILITIES_EXPORT Matrix(unsigned int rows, unsigned int columns,
                                        MatrixStorage policy = eFULL);
            
            LIB_UTILITIES_EXPORT Matrix(const Matrix<DataType>& rhs);

            LIB_UTILITIES_EXPORT Matrix<DataType>& operator=(const Matrix<DataType>& rhs);
            
            LIB_UTILITIES_EXPORT Matrix<DataType>& operator=(const ConstMatrix<DataType>& rhs);
            
        private:
            LIB_UTILITIES_EXPORT virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d) = 0;
    };
    
    
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP
