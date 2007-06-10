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
#include <LibUtilities/BasicUtils/BinaryExpressionTraits.hpp>

#include <boost/shared_ptr.hpp>

namespace Nektar
{
    template<typename DataType, typename StorageType, typename OwnedMatrixType>
    class NekMatrix<NekMatrix<DataType, StorageType, OwnedMatrixType>, StorageType, ScaledMatrixTag> : public ConstMatrix<DataType>
    {
        public:
            typedef ConstMatrix<DataType> BaseType;
            typedef NekMatrix<NekMatrix<DataType, StorageType, OwnedMatrixType>, StorageType, ScaledMatrixTag> ThisType;
            
            NekMatrix(typename boost::call_traits<DataType>::const_reference scale,
                      boost::shared_ptr<NekMatrix<DataType, StorageType, OwnedMatrixType> > m) :
                BaseType(m->GetRows(), m->GetColumns()),
                m_matrix(m),
                m_scale(scale)
            {
            }
            
            typename boost::call_traits<DataType>::value_type operator()(unsigned int row, unsigned int col) const
            {
                return m_scale*(*m_matrix)(row, col);
            }
            
            unsigned int GetStorageSize() const 
            {
                return m_matrix->GetStorageSize();
            }
            
            MatrixStorage GetStorageType() const
            {
                return m_matrix->GetStorageType();
            }
            
            typename boost::call_traits<DataType>::reference Scale()
            {
                return m_scale;
            }
            
            typename boost::call_traits<DataType>::const_reference Scale() const
            {
                return m_scale;
            }
            
            boost::shared_ptr<const NekMatrix<DataType, StorageType, OwnedMatrixType> > GetOwnedMatrix() const
            {
                return m_matrix; 
            }
            
        public:
        
        private:
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const 
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
            

            boost::shared_ptr<const NekMatrix<DataType, StorageType, OwnedMatrixType> > m_matrix;
            DataType m_scale;
    };
    
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SCALED_MATRIX_HPP
