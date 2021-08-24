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

#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

namespace Nektar
{
    template<typename DataType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, InnerMatrixType>, ScaledMatrixTag> : public ConstMatrix<DataType>
    {
        public:
            typedef ConstMatrix<DataType> BaseType;
            typedef NekMatrix<DataType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, ScaledMatrixTag> ThisType;
            typedef typename std::remove_const<typename InnerType::NumberType>::type NumberType;
            
            typedef NumberType GetValueType;
            typedef NumberType ConstGetValueType;
            
        public:
            /// \internal
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
                        ++m_iter;
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
            LIB_UTILITIES_EXPORT NekMatrix() ;
            
            LIB_UTILITIES_EXPORT NekMatrix(typename boost::call_traits<NumberType>::const_reference scale,
                      std::shared_ptr<const InnerType> m);
            
            LIB_UTILITIES_EXPORT NekMatrix(const ThisType& rhs);


            LIB_UTILITIES_EXPORT NekMatrix(typename boost::call_traits<NumberType>::const_reference scale, const ThisType& rhs);
           

            LIB_UTILITIES_EXPORT ThisType& operator=(const ThisType&) = default;

            LIB_UTILITIES_EXPORT ConstGetValueType operator()(unsigned int row, unsigned int col) const;

            LIB_UTILITIES_EXPORT unsigned int GetStorageSize() const ;
                                    
            LIB_UTILITIES_EXPORT NumberType Scale() const;
            LIB_UTILITIES_EXPORT void SetScale(const NumberType&);
            
            LIB_UTILITIES_EXPORT const NumberType* GetRawPtr() const;
            
            LIB_UTILITIES_EXPORT std::shared_ptr<const InnerType> GetOwnedMatrix() const;
            
            LIB_UTILITIES_EXPORT unsigned int GetNumberOfSubDiagonals() const;
            LIB_UTILITIES_EXPORT unsigned int GetNumberOfSuperDiagonals() const;
            
            LIB_UTILITIES_EXPORT const_iterator begin() const;
            LIB_UTILITIES_EXPORT const_iterator end() const;
            
            LIB_UTILITIES_EXPORT static ThisType CreateWrapper(const ThisType& rhs);
            
            LIB_UTILITIES_EXPORT static std::shared_ptr<ThisType> CreateWrapper(const std::shared_ptr<ThisType>& rhs);
            
        public:
        
        private:
            LIB_UTILITIES_EXPORT virtual typename boost::call_traits<NumberType>::value_type
            v_GetValue(unsigned int row, unsigned int column) const;
            
            LIB_UTILITIES_EXPORT virtual unsigned int v_GetStorageSize() const ;
            
            LIB_UTILITIES_EXPORT virtual char v_GetTransposeFlag() const;

            std::shared_ptr<const InnerType> m_matrix;
            NumberType m_scale;
    };

    template<typename DataType>
    LIB_UTILITIES_EXPORT void NegateInPlace(NekMatrix<DataType, ScaledMatrixTag>& v);
    
    template<typename DataType>
    LIB_UTILITIES_EXPORT NekMatrix<DataType, ScaledMatrixTag>
    Transpose(NekMatrix<DataType, ScaledMatrixTag>& rhs);

}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SCALED_MATRIX_HPP
