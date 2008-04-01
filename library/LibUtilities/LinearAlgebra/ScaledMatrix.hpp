///////////////////////////////////////////////////////////////////////////////
//
// File: ScaledMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SCALED_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SCALED_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/Memory/DeleteNothing.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

#include <boost/shared_ptr.hpp>

namespace Nektar
{
    template<typename DataType, typename InnerStorageType, typename InnerMatrixType, typename StorageType>
    class NekMatrix<NekMatrix<DataType, InnerStorageType, InnerMatrixType>, StorageType, ScaledMatrixTag> : public ConstMatrix<DataType>
    {
        public:
            typedef ConstMatrix<DataType> BaseType;
            typedef NekMatrix<DataType, InnerStorageType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, StorageType, ScaledMatrixTag> ThisType;
            typedef typename InnerType::NumberType NumberType;
            
            typedef NumberType GetValueType;
            typedef NumberType ConstGetValueType;
            
        public:
            /// \brief Iterator through elements of the matrix.
            /// Not quite a real iterator in the C++ sense.
            class const_iterator
            {
                public:
                    const_iterator(typename InnerType::const_iterator iter,
                                   const NumberType& scale) :
                        m_iter(iter),
                        m_scale(scale)
                    {
                    }
                    
                    const_iterator operator++(int)
                    {
                        const_iterator out = *this;
                        m_iter++;
                        return out;
                    }
                    
                    const_iterator& operator++()
                    {
                        ++m_iter;
                        return *this;
                    }
                    
                    NumberType operator*()
                    {
                        return m_scale*(*m_iter);
                    }
                    
                    bool operator==(const const_iterator& rhs)
                    {
                        return m_iter == rhs.m_iter;
                    }
                    
                    bool operator!=(const const_iterator& rhs)
                    {
                        return !(*this == rhs);
                    }
                    
                private:
                    typename InnerType::const_iterator m_iter;
                    NumberType m_scale;
            };
            
            
            
        public:
            NekMatrix() :
                BaseType(0,0),
                m_matrix(new InnerType()),
                m_scale(0)
            {
            }
            
            NekMatrix(typename boost::call_traits<NumberType>::const_reference scale,
                      boost::shared_ptr<const InnerType> m) :
                BaseType(m->GetRows(), m->GetColumns()),
                m_matrix(m),
                m_scale(scale)
            {
            }
            
            NekMatrix(const NekMatrix<InnerType, StorageType, ScaledMatrixTag>& rhs) :
                BaseType(rhs),
                m_matrix(rhs.m_matrix),
                m_scale(rhs.m_scale)
            {
            }
            
            typename boost::call_traits<NumberType>::value_type operator()(unsigned int row, unsigned int col) const
            {
                return m_scale*m_matrix->GetValue(row, col, this->GetTransposeFlag());
            }
            
            unsigned int GetStorageSize() const 
            {
                return m_matrix->GetStorageSize();
            }
            
            MatrixStorage GetStorageType() const
            {
                return m_matrix->GetStorageType();
            }
                        
            NumberType Scale() const
            {
                return m_scale*m_matrix->Scale();
            }
            
            boost::shared_ptr<const InnerType> GetOwnedMatrix() const
            {
                return m_matrix; 
            }
            
            const_iterator begin() const { return const_iterator(m_matrix->begin(this->GetTransposeFlag()), m_scale); }
            const_iterator end() const { return const_iterator(m_matrix->end(this->GetTransposeFlag()), m_scale); }
            
        public:
        
        private:
            virtual typename boost::call_traits<NumberType>::value_type 
            v_GetValue(unsigned int row, unsigned int column) const 
            {
                return ThisType::operator()(row, column);
            }
            
            virtual unsigned int v_GetStorageSize() const 
            {
                return ThisType::GetStorageSize();
            }
            
            virtual MatrixStorage v_GetStorageType() const
            {
                return ThisType::GetStorageType();
            }
            
            virtual char v_GetTransposeFlag() const
            {
                if( this->GetRawTransposeFlag() == 'N' )
                {
                    return m_matrix->GetTransposeFlag();
                }
                else
                {
                    if( m_matrix->GetTransposeFlag() == 'N' )
                    {
                        return 'T';
                    }
                    else
                    {
                        return 'N';
                    }
                }
            }

            boost::shared_ptr<const InnerType> m_matrix;
            NumberType m_scale;
    };
    
    template<typename DataType, typename StorageType>
    NekMatrix<DataType, StorageType, ScaledMatrixTag>
    Transpose(NekMatrix<DataType, StorageType, ScaledMatrixTag>& rhs)
    {
        NekMatrix<DataType, StorageType, ScaledMatrixTag> result(rhs);
        result.Transpose();
        return result;
    }
    

}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SCALED_MATRIX_HPP
