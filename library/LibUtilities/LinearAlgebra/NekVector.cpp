///////////////////////////////////////////////////////////////////////////////
//
// File: NekVector.cpp
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

#include <LibUtilities/LinearAlgebra/NekVector.hpp>

namespace Nektar
{

    template<typename DataType>
    NekVector<DataType>::NekVector() :
        m_size(0),
        m_data(),
        m_wrapperType(eCopy)
    {
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(unsigned int size) :
        m_size(size),
        m_data(size),
        m_wrapperType(eCopy)
    {
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(unsigned int size, typename boost::call_traits<DataType>::const_reference a) :
        m_size(size),
        m_data(size),
        m_wrapperType(eCopy)
    {
        std::fill_n(m_data.get(), m_size, a);
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(const std::string& vectorValues) :
        m_size(0),
        m_data(),
        m_wrapperType(eCopy)
    {
        try
        {
            std::vector<DataType> values = FromString<DataType>(vectorValues);
            m_size = values.size();
            m_data = Array<OneD, DataType>(m_size);
            std::copy(values.begin(), values.end(), m_data.begin());

            ASSERTL0(m_size > 0, "Error converting string values to vector");
        }
        catch(std::runtime_error& e)
        {
            NEKERROR(ErrorUtil::efatal, e.what());
        }
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(typename boost::call_traits<DataType>::const_reference x,
              typename boost::call_traits<DataType>::const_reference y,
              typename boost::call_traits<DataType>::const_reference z) :
        m_size(3),
        m_data(m_size),
        m_wrapperType(eCopy)
    {
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(const NekVector<DataType>& rhs) :
        m_size(rhs.GetDimension()),
        m_data(rhs.m_data),
        m_wrapperType(rhs.m_wrapperType)
    {
        if( m_wrapperType == eCopy )
        {
            m_data = Array<OneD, DataType>(m_size);
            std::copy(rhs.begin(), rhs.end(), m_data.get());
        }
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(unsigned int size, const DataType* const ptr) :
        m_size(size),
        m_data(size, ptr),
        m_wrapperType(eCopy)
    {
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(const Array<OneD, DataType>& ptr, PointerWrapper h) :
        m_size(ptr.size()),
        m_data(ptr),
        m_wrapperType(h)
    {
        if( h == eCopy )
        {
            m_data = Array<OneD, DataType>(m_size);
            CopyArray(ptr, m_data);
        }
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(unsigned int size, Array<OneD, DataType>& ptr, PointerWrapper h) :
        m_size(size),
        m_data(ptr),
        m_wrapperType(h)
    {
        if( h == eCopy )
        {
            ASSERTL0(size <= ptr.size(), "Attempting to populate a vector of size " +
                std::to_string(size) + " but the incoming array only has " +
                std::to_string(ptr.size()) + " elements.");

            m_data = Array<OneD, DataType>(size);
            std::copy(ptr.begin(), ptr.begin()+size, m_data.begin());
        }
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(const Array<OneD, const DataType>& ptr, PointerWrapper h) :
        m_size(ptr.size()),
        m_data(ptr, eVECTOR_WRAPPER),
        m_wrapperType(h)
    {
        if( h == eCopy )
        {
            m_data = Array<OneD, DataType>(m_size);
            CopyArray(ptr, m_data);
        }
    }

    template<typename DataType>
    NekVector<DataType>::NekVector(unsigned int size, const Array<OneD, const DataType>& ptr, PointerWrapper h) :
        m_size(size),
        m_data(ptr, eVECTOR_WRAPPER),
        m_wrapperType(h)
    {
        if( h == eCopy )
        {
            ASSERTL0(size <= ptr.size(), "Attempting to populate a vector of size " +
                std::to_string(size) + " but the incoming array only has " +
                std::to_string(ptr.size()) + " elements.");

            m_data = Array<OneD, DataType>(size);
            std::copy(ptr.begin(), ptr.begin()+size, m_data.begin());
        }
    }

    template<typename DataType>
    NekVector<DataType>::~NekVector() {}

    template<typename DataType>
    NekVector<DataType>& NekVector<DataType>::operator=(const NekVector<DataType>& rhs)
    {
        if( m_wrapperType == eCopy  )
        {
            // If the current vector is a copy, then regardless of the rhs type
            // we just copy over the values, resizing if needed.
            if( GetDimension() != rhs.GetDimension() )
            {
                m_size = rhs.GetDimension();
                m_data = Array<OneD, DataType>(m_size);
            }
        }
        else if( m_wrapperType == eWrapper )
        {
            // If the current vector is wrapped, then just copy over the top,
            // but the sizes of the two vectors must be the same.
            ASSERTL0(GetDimension() == rhs.GetDimension(), "Wrapped NekVectors must have the same dimension in operator=");
        }

        std::copy(rhs.begin(), rhs.end(), m_data.get());
        return *this;
    }


    template<typename DataType>
    unsigned int NekVector<DataType>::GetDimension() const
    {
        return m_size;
    }

    template<typename DataType>
    unsigned int NekVector<DataType>::GetRows() const
    {
        return m_size;
    }

    template<typename DataType>
    DataType* NekVector<DataType>::GetRawPtr()
    {
        return this->GetData().get();
    }

    template<typename DataType>
    Array<OneD, DataType>& NekVector<DataType>::GetPtr() { return this->GetData(); }

    template<typename DataType>
    const DataType* NekVector<DataType>::GetRawPtr() const
    {
        return m_data.get();
    }

    template<typename DataType>
    const Array<OneD, const DataType>& NekVector<DataType>::GetPtr() const { return m_data; }

    template<typename DataType>
    typename NekVector<DataType>::iterator NekVector<DataType>::begin() { return GetRawPtr(); }

    template<typename DataType>
    typename NekVector<DataType>::iterator NekVector<DataType>::end() { return GetRawPtr() + this->GetDimension(); }

    template<typename DataType>
    typename NekVector<DataType>::const_iterator NekVector<DataType>::begin() const { return GetRawPtr(); }

    template<typename DataType>
    typename NekVector<DataType>::const_iterator NekVector<DataType>::end() const { return GetRawPtr() + GetDimension(); }

    template<typename DataType>
    typename boost::call_traits<DataType>::reference NekVector<DataType>::operator()(unsigned int i)
    {
        ASSERTL1(i < this->GetDimension(),
                 "Invalid access to m_data via parenthesis operator");
        return this->GetData()[i];
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::reference NekVector<DataType>::operator[](unsigned int i)
    {
        return this->GetData()[i];
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::reference NekVector<DataType>::x()
    {
        ASSERTL1(this->GetDimension() >= 1, "Invalid use of NekVector::x");
        return (*this)(0);
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::reference NekVector<DataType>::y()
    {
        ASSERTL1(this->GetDimension() >= 2, "Invalid use of NekVector::y");
        return (*this)(1);
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::reference NekVector<DataType>::z()
    {
        ASSERTL1(this->GetDimension() >= 3, "Invalid use of NekVector::z");
        return (*this)(2);
    }

    template<typename DataType>
    void NekVector<DataType>::SetX(typename boost::call_traits<DataType>::const_reference val)
    {
        ASSERTL1(this->GetDimension() >= 1, "Invalid use of NekVector::SetX");
        this->GetData()[0] = val;
    }

    template<typename DataType>
    void NekVector<DataType>::SetY(typename boost::call_traits<DataType>::const_reference val)
    {
        ASSERTL1(this->GetDimension() >= 2, "Invalid use of NekVector::SetX");
        this->GetData()[1] = val;
    }

    template<typename DataType>
    void NekVector<DataType>::SetZ(typename boost::call_traits<DataType>::const_reference val)
    {
        ASSERTL1(this->GetDimension() >= 3, "Invalid use of NekVector::SetX");
        this->GetData()[2] = val;
    }

    template<typename DataType>
    NekVector<DataType>& NekVector<DataType>::operator+=(const NekVector<DataType>& rhs)
    {
        AddEqual(*this, rhs);
        return *this;
    }

    template<typename DataType>
    NekVector<DataType>& NekVector<DataType>::operator-=(const NekVector<DataType>& rhs)
    {
        SubtractEqual(*this, rhs);
        return *this;
    }

    template<typename DataType>
    NekVector<DataType>& NekVector<DataType>::operator*=(typename boost::call_traits<DataType>::const_reference rhs)
    {
        MultiplyEqual(*this, rhs);
        return *this;
    }

    template<typename DataType>
    NekVector<DataType>& NekVector<DataType>::operator/=(typename boost::call_traits<DataType>::const_reference rhs)
    {
        DivideEqual(*this, rhs);
        return *this;
    }

    template<typename DataType>
    void NekVector<DataType>::Normalize() { return Nektar::Normalize(*this); }

    template<typename DataType>
    typename boost::call_traits<DataType>::const_reference NekVector<DataType>::operator()(unsigned int i) const
    {
        ASSERTL1(i < GetDimension(),
                 "Invalid access to m_data via parenthesis operator");
        return m_data[i];
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::const_reference NekVector<DataType>::operator[](unsigned int i) const
    {
        return m_data[i];
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::const_reference NekVector<DataType>::x() const
    {
        ASSERTL1( GetDimension() >= 1, "Invalid use of NekVector::x");
        return (*this)(0);
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::const_reference NekVector<DataType>::y() const
    {
        ASSERTL1( GetDimension() >= 2, "Invalid use of NekVector::y");
        return (*this)(1);
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::const_reference NekVector<DataType>::z() const
    {
        ASSERTL1( GetDimension() >= 3, "Invalid use of NekVector::z");
        return (*this)(2);
    }

    template<typename DataType>
    NekVector<DataType> NekVector<DataType>::operator-() const { return Negate(*this); }

    template<typename DataType>
    DataType NekVector<DataType>::Magnitude() const { return Nektar::Magnitude(*this); }

    template<typename DataType>
    DataType NekVector<DataType>::Dot(const NekVector<DataType>& rhs) const { return Nektar::Dot(*this, rhs); }

    template<typename DataType>
    NekVector<DataType> NekVector<DataType>::Cross(const NekVector<DataType>& rhs) const
    {
        return Nektar::Cross(*this, rhs);
    }

    template<typename DataType>
    std::string NekVector<DataType>::AsString() const { return Nektar::AsString(*this); }

    // Norms
    template<typename DataType>
    DataType NekVector<DataType>::L1Norm() const { return Nektar::L1Norm(*this); }

    template<typename DataType>
    DataType NekVector<DataType>::L2Norm() const { return Nektar::L2Norm(*this); }

    template<typename DataType>
    DataType NekVector<DataType>::InfinityNorm() const { return Nektar::InfinityNorm(*this); }

    template<typename DataType>
    PointerWrapper NekVector<DataType>::GetWrapperType() const { return m_wrapperType; }

    template<typename DataType>
    Array<OneD, DataType>& NekVector<DataType>::GetData() { return m_data; }

    template<typename DataType>
    void NekVector<DataType>::SetSize(unsigned int s) { m_size = s; }

    template<typename DataType>
    void NekVector<DataType>::SetWrapperType(PointerWrapper p) { m_wrapperType = p; }

    template<typename DataType>
    void NekVector<DataType>::SetData(const Array<OneD, DataType>& newData) { m_data = newData; }

    template<typename DataType>
    void NekVector<DataType>::Resize(unsigned int newSize)
    {
        if(m_data.size() < newSize )
        {
            m_data = Array<OneD, DataType>(newSize);
        }
        m_size = newSize;
    }

    template LIB_UTILITIES_EXPORT class NekVector<NekDouble>;

    template<typename DataType>
    void Add(NekVector<DataType>& result,
           const NekVector<DataType>& lhs,
           const NekVector<DataType>& rhs)
    {
        DataType* r_buf = result.GetRawPtr();
        const DataType* lhs_buf = lhs.GetRawPtr();
        const DataType* rhs_buf = rhs.GetRawPtr();
        const unsigned int ldim = lhs.GetDimension();
        for(int i = 0; i < ldim; ++i)
        {
            r_buf[i] = lhs_buf[i] + rhs_buf[i];
        }
    }

    template<typename DataType>
    void AddNegatedLhs(NekVector<DataType>& result,
           const NekVector<DataType>& lhs,
           const NekVector<DataType>& rhs)
    {
        DataType* r_buf = result.GetRawPtr();
        const DataType* lhs_buf = lhs.GetRawPtr();
        const DataType* rhs_buf = rhs.GetRawPtr();
        const unsigned int ldim = lhs.GetDimension();
        for(int i = 0; i < ldim; ++i)
        {
            r_buf[i] = -lhs_buf[i] + rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT void Add(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& lhs,
           const NekVector<NekDouble>& rhs);
    template LIB_UTILITIES_EXPORT void AddNegatedLhs(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& lhs,
           const NekVector<NekDouble>& rhs);

    template<typename DataType>
    void AddEqual(NekVector<DataType>& result,
           const NekVector<DataType>& rhs)
    {
        DataType* r_buf = result.GetRawPtr();
        const DataType* rhs_buf = rhs.GetRawPtr();
        const unsigned int rdim = rhs.GetDimension();
        for(int i = 0; i < rdim; ++i)
        {
            r_buf[i] += rhs_buf[i];
        }
    }

    template<typename DataType>
    void AddEqualNegatedLhs(NekVector<DataType>& result,
           const NekVector<DataType>& rhs)
    {
         DataType* r_buf = result.GetRawPtr();
        const DataType* rhs_buf = rhs.GetRawPtr();
        const unsigned int rdim = rhs.GetDimension();
        for(int i = 0; i < rdim; ++i)
        {
            r_buf[i] = -r_buf[i] + rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void AddEqual(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& rhs);
    template LIB_UTILITIES_EXPORT
    void AddEqualNegatedLhs(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& rhs);

    template<typename LhsDataType,
             typename RhsDataType>
    NekVector<LhsDataType> Add(const NekVector<LhsDataType>& lhs,
                               const NekVector<RhsDataType>& rhs)
    {
        NekVector<LhsDataType> result(lhs.GetDimension());
        Add(result, lhs, rhs);
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble> Add(const NekVector<NekDouble>& lhs,
                             const NekVector<NekDouble>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void Subtract(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekVector<InputDataType>& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename std::add_const<InputDataType>::type* lhs_buf = lhs.GetRawPtr();
        typename std::add_const<InputDataType>::type* rhs_buf = rhs.GetRawPtr();
        const unsigned int ldim = lhs.GetDimension();
        for(int i = 0; i < ldim; ++i)
        {
            r_buf[i] = lhs_buf[i] - rhs_buf[i];
        }
    }

    template<typename ResultDataType, typename InputDataType>
    void SubtractNegatedLhs(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekVector<InputDataType>& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename std::add_const<InputDataType>::type* lhs_buf = lhs.GetRawPtr();
        typename std::add_const<InputDataType>::type* rhs_buf = rhs.GetRawPtr();
        const unsigned int ldim = lhs.GetDimension();
        for(int i = 0; i < ldim; ++i)
        {
            r_buf[i] = -lhs_buf[i] - rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void Subtract(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& lhs,
           const NekVector<NekDouble>& rhs);

    template LIB_UTILITIES_EXPORT
    void SubtractNegatedLhs(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& lhs,
           const NekVector<NekDouble>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void SubtractEqual(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename std::add_const<InputDataType>::type* rhs_buf = rhs.GetRawPtr();
        const unsigned int rdim = rhs.GetDimension();
        for(int i = 0; i < rdim; ++i)
        {
            r_buf[i] -= rhs_buf[i];
        }
    }

    template<typename ResultDataType, typename InputDataType>
    void SubtractEqualNegatedLhs(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& rhs)
    {
                ResultDataType* r_buf = result.GetRawPtr();
        typename std::add_const<InputDataType>::type* rhs_buf = rhs.GetRawPtr();
        const unsigned int rdim = rhs.GetDimension();
        for(int i = 0; i < rdim; ++i)
        {
            r_buf[i] = -r_buf[i] - rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void SubtractEqual(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& rhs);

    template LIB_UTILITIES_EXPORT
    void SubtractEqualNegatedLhs(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& rhs);

    template<typename DataType>
    NekVector<DataType>
    Subtract(const NekVector<DataType>& lhs,
                const NekVector<DataType>& rhs)
    {
        NekVector<DataType> result(lhs.GetDimension());
        Subtract(result, lhs, rhs);
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble>
    Subtract(const NekVector<NekDouble>& lhs,
                const NekVector<NekDouble>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void Divide(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename std::add_const<InputDataType>::type* lhs_buf = lhs.GetRawPtr();

        const unsigned int ldim = lhs.GetDimension();
        for(int i = 0; i < ldim; ++i)
        {
            r_buf[i] = lhs_buf[i] / rhs;
        }
    }

    template LIB_UTILITIES_EXPORT
    void Divide(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& lhs,
           const NekDouble& rhs);

    template<typename ResultDataType>
    void DivideEqual(NekVector<ResultDataType>& result,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();

        const unsigned int resdim = result.GetDimension();
        for(int i = 0; i < resdim; ++i)
        {
            r_buf[i] /= rhs;
        }
    }

    template LIB_UTILITIES_EXPORT
    void DivideEqual(NekVector<NekDouble>& result,
                     const NekDouble& rhs);

    template<typename DataType>
    NekVector<DataType>
    Divide(const NekVector<DataType>& lhs,
                const NekDouble& rhs)
    {
        NekVector<DataType> result(lhs.GetDimension());
        Divide(result, lhs, rhs);
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble>
    Divide(const NekVector<NekDouble>& lhs,
                const NekDouble& rhs);


    template<typename ResultDataType, typename InputDataType>
    void Multiply(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekVector<InputDataType>& rhs)
    {
        ResultDataType* result_buf = result.GetRawPtr();
        const InputDataType* rhs_buf = rhs.GetRawPtr();
        const InputDataType* lhs_buf = lhs.GetRawPtr();
        const unsigned int resdim = result.GetDimension();
        for(int i = 0; i < resdim; ++i)
        {
            result_buf[i] = lhs_buf[i] * rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void Multiply(NekVector<NekDouble>& result, const NekVector<NekDouble>& lhs, const NekVector<NekDouble>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void MultiplyEqual(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& rhs)
    {
        ResultDataType* result_buf = result.GetRawPtr();
        const InputDataType* rhs_buf = rhs.GetRawPtr();
        const unsigned int resdim = result.GetDimension();
        for(int i = 0; i < resdim; ++i)
        {
            result_buf[i] *= rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void MultiplyEqual(NekVector<NekDouble>& result, const NekVector<NekDouble>& rhs);

    template<typename DataType, typename InputDataType>
    NekVector<DataType>
    Multiply(const NekVector<DataType>& lhs,
                const NekVector<InputDataType>& rhs)
    {
        NekVector<DataType> result(lhs.GetDimension());
        Multiply(result, lhs, rhs);
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble>
    Multiply(const NekVector<NekDouble>& lhs, const NekVector<NekDouble>& rhs);


    template<typename ResultDataType, typename InputDataType>
    void Multiply(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        const InputDataType* lhs_buf = lhs.GetRawPtr();

        const unsigned int ldim = lhs.GetDimension();
        for(int i = 0; i < ldim; ++i)
        {
            r_buf[i] = lhs_buf[i] * rhs;
        }
    }

    template LIB_UTILITIES_EXPORT
    void Multiply(NekVector<NekDouble>& result,
           const NekVector<NekDouble>& lhs,
           const NekDouble& rhs);

    template<typename ResultDataType>
    void MultiplyEqual(NekVector<ResultDataType>& result,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        const unsigned int rdim = result.GetDimension();
        for(unsigned int i = 0; i < rdim; ++i)
        {
            r_buf[i] *= rhs;
        }
    }
    template LIB_UTILITIES_EXPORT
    void MultiplyEqual(NekVector<NekDouble>& result,
           const NekDouble& rhs);

    template<typename DataType>
    NekVector<DataType>
    Multiply(const NekVector<DataType>& lhs,
                const NekDouble& rhs)
    {
        NekVector<DataType> result(lhs.GetDimension());
        Multiply(result, lhs, rhs);
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble>
    Multiply(const NekVector<NekDouble>& lhs,
                const NekDouble& rhs);

    template<typename ResultDataType, typename InputDataType>
    void Multiply(NekVector<ResultDataType>& result,
                        const NekDouble& lhs,
                        const NekVector<InputDataType>& rhs)
    {
                Multiply(result, rhs, lhs);
    }

    template<typename ResultDataType, typename InputDataType>
    void MultiplyInvertedLhs(NekVector<ResultDataType>& result,
                  const NekDouble& lhs,
                  const NekVector<InputDataType>& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        const InputDataType* rhs_buf = rhs.GetRawPtr();
        NekDouble inverse = 1.0/lhs;

        const unsigned int rdim = rhs.GetDimension();
        for(int i = 0; i < rdim; ++i)
        {
            r_buf[i] = inverse * rhs_buf[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void MultiplyInvertedLhs(NekVector<NekDouble>& result,
                        const NekDouble& lhs,
                        const NekVector<NekDouble>& rhs);

    template LIB_UTILITIES_EXPORT
    void Multiply(NekVector<NekDouble>& result,
                        const NekDouble& lhs,
                        const NekVector<NekDouble>& rhs);

    template<typename DataType>
    NekVector<DataType> Multiply(const DataType& lhs,
                                 const NekVector<DataType>& rhs)
    {
                return Multiply(rhs, lhs);
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble> Multiply(const NekDouble& lhs,
                                  const NekVector<NekDouble>& rhs);


    template<typename DataType>
    std::ostream& operator<<(std::ostream& os, const NekVector<DataType>& rhs)
    {
        os << rhs.AsString();
        return os;
    }

    template LIB_UTILITIES_EXPORT
    std::ostream& operator<<(std::ostream& os, const NekVector<NekDouble>& rhs);

    template<typename DataType>
    NekVector<DataType> createVectorFromPoints(const NekPoint<DataType>& source,
                                               const NekPoint<DataType>& dest)
    {
        NekVector<DataType> result(3, 0.0);
        for(unsigned int i = 0; i < 3; ++i)
        {
            result[i] = dest[i]-source[i];
        }
        return result;
    }


    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble> createVectorFromPoints(const NekPoint<NekDouble>& source,
                                               const NekPoint<NekDouble>& dest);

    template<typename DataType>
    NekPoint<DataType> findPointAlongVector(const NekVector<DataType>& lhs,
                                            const DataType& t)
    {
        NekPoint<DataType> result;
        for(unsigned int i = 0; i < 3; ++i)
        {
            result[i] = lhs[i]*t;
        }

        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekPoint<NekDouble> findPointAlongVector(const NekVector<NekDouble>& lhs,
                                            const NekDouble& t);

    template<typename DataType>
    bool operator==(const NekVector<DataType>& lhs,
                    const NekVector<DataType>& rhs)
    {
        if( lhs.GetDimension() != rhs.GetDimension() )
        {
            return false;
        }

        return std::equal(lhs.begin(), lhs.end(), rhs.begin());
    }

    template LIB_UTILITIES_EXPORT
    bool operator==(const NekVector<NekDouble>& lhs,
                    const NekVector<NekDouble>& rhs);

    template<typename DataType>
    bool operator!=(const NekVector<DataType>& lhs,
                    const NekVector<DataType>& rhs)
    {
        return !(lhs == rhs);
    }

    template LIB_UTILITIES_EXPORT
    bool operator!=(const NekVector<NekDouble>& lhs,
                    const NekVector<NekDouble>& rhs);

    template<typename DataType>
    std::vector<DataType> FromString(const std::string& str)
    {
        std::vector<DataType> result;

        try
        {
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            boost::char_separator<char> sep("(<,>) ");
            tokenizer tokens(str, sep);
            for( tokenizer::iterator strIter = tokens.begin(); strIter != tokens.end(); ++strIter)
            {
                result.push_back(boost::lexical_cast<DataType>(*strIter));
            }
        }
        catch(boost::bad_lexical_cast&)
        {
        }

        return result;
    }

    template LIB_UTILITIES_EXPORT
    std::vector<NekDouble> FromString(const std::string& str);

    template<typename DataType>
    DataType L1Norm(const NekVector<DataType>& v)
    {
        typedef NekVector<DataType> VectorType;

        DataType result(0);
        for(typename VectorType::const_iterator iter = v.begin(); iter != v.end(); ++iter)
        {
            result += fabs(*iter);
        }

        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekDouble L1Norm(const NekVector<NekDouble>& v);

    template<typename DataType>
    DataType L2Norm(const NekVector<DataType>& v)
    {
        typedef NekVector<DataType> VectorType;

        DataType result(0);
        for(typename VectorType::const_iterator iter = v.begin(); iter != v.end(); ++iter)
        {
            DataType v = fabs(*iter);
            result += v*v;
        }
        return sqrt(result);
    }

    template LIB_UTILITIES_EXPORT
    NekDouble L2Norm(const NekVector<NekDouble>& v);

    template<typename DataType>
    DataType InfinityNorm(const NekVector<DataType>& v)
    {
        DataType result = fabs(v[0]);
        const unsigned int vdim = v.GetDimension();
        for(unsigned int i = 0; i < vdim; ++i)
        {
            result = std::max(fabs(v[i]), result);
        }
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekDouble InfinityNorm(const NekVector<NekDouble>& v);

    template<typename DataType>
    NekVector<DataType> Negate(const NekVector<DataType>& v)
    {
        NekVector<DataType> temp(v);
        const unsigned int tdim = temp.GetDimension();
        for(unsigned int i = 0; i < tdim; ++i)
        {
            temp(i) = -temp(i);
        }
        return temp;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble> Negate(const NekVector<NekDouble>& v);

    template<typename DataType>
    void NegateInPlace(NekVector<DataType>& v)
    {
        DataType* data = v.GetRawPtr();
        const unsigned int vdim = v.GetDimension();
        for(unsigned int i = 0; i < vdim; ++i)
        {
            data[i] = -data[i];
        }
    }

    template LIB_UTILITIES_EXPORT
    void NegateInPlace(NekVector<NekDouble>& v);

    template<typename DataType>
    DataType Magnitude(const NekVector<DataType>& v)
    {
        DataType result = DataType(0);

        const unsigned int vdim = v.GetDimension();
        for(unsigned int i = 0; i < vdim; ++i)
        {
            result += v[i]*v[i];
        }
        return sqrt(result);
    }

    template LIB_UTILITIES_EXPORT
    NekDouble Magnitude(const NekVector<NekDouble>& v) ;

    template<typename DataType>
    DataType Dot(const NekVector<DataType>& lhs,
                 const NekVector<DataType>& rhs)
    {
        ASSERTL1( lhs.GetDimension() == rhs.GetDimension(), "Dot, dimension of the two operands must be identical.");

        DataType result = DataType(0);
        const unsigned int ldim = lhs.GetDimension();
        for(unsigned int i = 0; i < ldim; ++i)
        {
            result += lhs[i]*rhs[i];
        }

        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekDouble Dot(const NekVector<NekDouble>& lhs,
                 const NekVector<NekDouble>& rhs) ;

    template<typename DataType>
    void Normalize(NekVector<DataType>& v)
    {
        DataType m = v.Magnitude();
        if( m > DataType(0) )
        {
            v /= m;
        }
    }

    template LIB_UTILITIES_EXPORT
    void Normalize(NekVector<NekDouble>& v);

    void NegateInPlace(NekDouble& v) { v = -v; }
    void InvertInPlace(NekDouble& v) { v = 1.0/v; }

    template<typename DataType>
    NekVector<DataType> Cross(const NekVector<DataType>& lhs,
                                          const NekVector<DataType>& rhs)
    {
        ASSERTL1(lhs.GetDimension() == 3 && rhs.GetDimension() == 3, "Cross is only valid for 3D vectors.");

        DataType first = lhs.y()*rhs.z() - lhs.z()*rhs.y();
        DataType second = lhs.z()*rhs.x() - lhs.x()*rhs.z();
        DataType third = lhs.x()*rhs.y() - lhs.y()*rhs.x();

        NekVector<DataType> result(first, second, third);
        return result;
    }

    template LIB_UTILITIES_EXPORT
    NekVector<NekDouble> Cross(const NekVector<NekDouble>& lhs, const NekVector<NekDouble>& rhs);

    template<typename DataType>
    std::string AsString(const NekVector<DataType>& v)
    {
        unsigned int d = v.GetRows();
        std::string result = "(";
        for(unsigned int i = 0; i < d; ++i)
        {
            result += boost::lexical_cast<std::string>(v[i]);
            if( i < v.GetDimension()-1 )
            {
                result += ", ";
            }
        }
        result += ")";
        return result;
    }

    template LIB_UTILITIES_EXPORT
    std::string AsString(const NekVector<NekDouble>& v);
}

