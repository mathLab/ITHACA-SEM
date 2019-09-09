///////////////////////////////////////////////////////////////////////////////
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

#include <LibUtilities/LinearAlgebra/MatrixType.h>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>

namespace Nektar
{
    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NekMatrix() :
        NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::BaseType(0,0,eFULL),
        m_matrix(new typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::InnerType()),
        m_scale(0)
    {
    }

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NekMatrix(typename boost::call_traits<NumberType>::const_reference scale,
              std::shared_ptr<const InnerType> m) :
        NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::BaseType(m->GetRows(), m->GetColumns(),m->GetStorageType()),
        m_matrix(m),
        m_scale(scale)
    {
    }

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NekMatrix(const NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>& rhs) :
        NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::BaseType(rhs),
        m_matrix(rhs.m_matrix),
        m_scale(rhs.m_scale)
    {
    }


    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NekMatrix(
            typename boost::call_traits<typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NumberType>::const_reference scale, const NekMatrix<InnerType, ScaledMatrixTag>& rhs) :
        NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::BaseType(rhs),
        m_matrix(rhs.m_matrix),
        m_scale(scale)
    {
    }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::ConstGetValueType
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::operator()(unsigned int row, unsigned int col) const
    {
        return m_scale*m_matrix->GetValue(row, col, this->GetTransposeFlag());
    }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::GetStorageSize() const
    {
        return m_matrix->GetStorageSize();
    }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NumberType
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::Scale() const
    {
        return m_scale*m_matrix->Scale();
    }

    template<typename DataType, typename InnerMatrixType>
    void NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::SetScale(const NumberType& value)
    {
        m_scale = value;
    }
    
    template<typename DataType, typename InnerMatrixType>
    const typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NumberType*
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::GetRawPtr() const { return m_matrix->GetRawPtr(); }

    template<typename DataType, typename InnerMatrixType>
    std::shared_ptr<const typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::InnerType>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::GetOwnedMatrix() const
    {
        return m_matrix;
    }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::GetNumberOfSubDiagonals() const { return m_matrix->GetNumberOfSubDiagonals(); }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::GetNumberOfSuperDiagonals() const { return m_matrix->GetNumberOfSuperDiagonals(); }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::begin() const { return const_iterator(m_matrix->begin(this->GetTransposeFlag()), m_scale); }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::end() const { return const_iterator(m_matrix->end(this->GetTransposeFlag()), m_scale); }

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::CreateWrapper(const NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>& rhs)
    {
        return NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>(rhs);
    }

    template<typename DataType, typename InnerMatrixType>
    std::shared_ptr<NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag> >
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::CreateWrapper(const std::shared_ptr<NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag> >& rhs)
    {
        return std::shared_ptr<ThisType>(new ThisType(*rhs));
    }


    template<typename DataType, typename InnerMatrixType>
    typename boost::call_traits<typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NumberType>::value_type
    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::v_GetValue(unsigned int row, unsigned int column) const
    {
        return ThisType::operator()(row, column);
    }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::v_GetStorageSize() const
    {
        return ThisType::GetStorageSize();
    }

    template<typename DataType, typename InnerMatrixType>
    char NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::v_GetTransposeFlag() const
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




    template<typename DataType>
    NekMatrix<DataType, ScaledMatrixTag>
    Transpose(NekMatrix<DataType, ScaledMatrixTag>& rhs)
    {
        NekMatrix<DataType, ScaledMatrixTag> result(rhs);
        result.Transpose();
        return result;
    }


//    template<typename DataType, typename InnerMatrixType>
//    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator::const_iterator(
//            typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::InnerType::const_iterator iter,
//            const typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NumberType& scale) :
//        m_iter(iter),
//        m_scale(scale)
//    {
//    }

//    template<typename DataType, typename InnerMatrixType>
//    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator::operator++(int)
//    {
//        const_iterator out = *this;
//        ++m_iter;
//        return out;
//    }

//    template<typename DataType, typename InnerMatrixType>
//    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator& NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator::operator++()
//    {
//        ++m_iter;
//        return *this;
//    }

//    template<typename DataType, typename InnerMatrixType>
//    typename NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::NumberType
//    NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator::operator*()
//    {
//        return m_scale*(*m_iter);
//    }

//    template<typename DataType, typename InnerMatrixType>
//    bool NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator::operator==(const NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator& rhs)
//    {
//        return m_iter == rhs.m_iter;
//    }

//    template<typename DataType, typename InnerMatrixType>
//    bool NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator::operator!=(const NekMatrix<NekMatrix< DataType, InnerMatrixType>, ScaledMatrixTag>::const_iterator& rhs)
//    {
//        return !(*this == rhs);
//    }



    template LIB_UTILITIES_EXPORT class NekMatrix<NekMatrix< NekDouble, StandardMatrixTag>, ScaledMatrixTag>;

    template
    LIB_UTILITIES_EXPORT NekMatrix< NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>
    Transpose(NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void NegateInPlace(NekMatrix<DataType, ScaledMatrixTag>& v)
    {
        v.SetScale(-1.0*v.Scale());
    }

    template LIB_UTILITIES_EXPORT void NegateInPlace(NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>& v);
}

