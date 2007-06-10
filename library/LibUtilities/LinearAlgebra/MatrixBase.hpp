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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP

#include <LibUtilities/LinearAlgebra/MatrixType.h>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <boost/lexical_cast.hpp>

namespace Nektar
{

    template<typename DataType>
    class ConstMatrix
    {
        public:
            typedef ConstMatrix<DataType> ThisType;
            
        public:
            typename boost::call_traits<DataType>::value_type operator()(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetColumns()) +
                    std::string(" columns"));
                return v_GetValue(row, column);
            }
            
            unsigned int GetStorageSize() const
            {
                return v_GetStorageSize();
            }
            
            MatrixStorage GetStorageType() const 
            {
                return v_GetStorageType();
            }
            
            unsigned int GetRows() const { return m_size[0]; }
            unsigned int GetColumns() const { return m_size[1]; }
            const unsigned int* GetSize() const { return m_size; }
            
        protected:
            
            
            
        protected:
            
            // All constructors are private to enforce the abstract nature of ConstMatrix without
            // resorting to pure virtual functions.
            ConstMatrix(unsigned int rows, unsigned int columns) :
                m_size()
            {
                m_size[0] = rows;
                m_size[1] = columns;
            }
            
            ConstMatrix(const ThisType& rhs) :
                m_size()
            {
                m_size[0] = rhs.GetRows();
                m_size[1] = rhs.GetColumns();
            }

            ThisType& operator=(const ThisType& rhs)
            {
                m_size[0] = rhs.m_size[0];
                m_size[1] = rhs.m_size[1];
                return *this;
            }
            
        private:
            
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const = 0;            
            virtual unsigned int v_GetStorageSize() const = 0;            
            virtual MatrixStorage v_GetStorageType() const = 0;
            
            unsigned int m_size[2];
    };
    
    template<typename DataType>
    class Matrix : public ConstMatrix<DataType>
    {
        public:
            typedef Matrix<DataType> ThisType;
            typedef ConstMatrix<DataType> BaseType;
            
        public:
            void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                v_SetValue(row, column, d);
            }
            
        protected:            
            // All constructors are private to enforce the abstract nature of ConstMatrix without
            // resorting to pure virtual functions.
            Matrix(unsigned int rows, unsigned int columns) :
                BaseType(rows, columns)
            {
            }
            
            Matrix(const ThisType& rhs) :
                BaseType(rhs)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
            
        private:
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d) = 0;
    };
    
    
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP
