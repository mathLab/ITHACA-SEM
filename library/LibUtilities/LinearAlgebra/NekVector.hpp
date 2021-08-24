///////////////////////////////////////////////////////////////////////////////
//
// File: NekVector.hpp
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
// Description: Generic N-Dimensional Vector.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP
#define NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/LinearAlgebra/NekPoint.hpp>

#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>

#include <functional>
#include <algorithm>
#include <cmath>

#include <type_traits>
#include <boost/call_traits.hpp>

namespace Nektar
{
    template<typename DataType>
    class NekVector
    {
        public:
            /// \brief Creates an empty vector.
            LIB_UTILITIES_EXPORT NekVector();

            /// \brief Creates a vector of given size.  The elements are not initialized.
            LIB_UTILITIES_EXPORT explicit NekVector(unsigned int size);

            /// \brief Creates a vector with given size and initial value.
            LIB_UTILITIES_EXPORT NekVector(unsigned int size, typename boost::call_traits<DataType>::const_reference a);


            LIB_UTILITIES_EXPORT explicit NekVector(const std::string& vectorValues);

            LIB_UTILITIES_EXPORT NekVector(typename boost::call_traits<DataType>::const_reference x,
                      typename boost::call_traits<DataType>::const_reference y,
                      typename boost::call_traits<DataType>::const_reference z);

            LIB_UTILITIES_EXPORT NekVector(const NekVector<DataType>& rhs);

            LIB_UTILITIES_EXPORT NekVector(unsigned int size, const DataType* const ptr);
            LIB_UTILITIES_EXPORT explicit NekVector(const Array<OneD, DataType>& ptr, PointerWrapper h = eCopy);
            LIB_UTILITIES_EXPORT NekVector(unsigned int size, Array<OneD, DataType>& ptr, PointerWrapper h = eCopy);

            LIB_UTILITIES_EXPORT NekVector(unsigned int size, const Array<OneD, const DataType>& ptr, PointerWrapper h = eCopy);

            LIB_UTILITIES_EXPORT ~NekVector();

            LIB_UTILITIES_EXPORT NekVector<DataType>& operator=(const NekVector<DataType>& rhs);


            /// \brief Returns the number of dimensions for the point.
            LIB_UTILITIES_EXPORT unsigned int GetDimension() const;

            LIB_UTILITIES_EXPORT unsigned int GetRows() const;

            LIB_UTILITIES_EXPORT DataType* GetRawPtr();

            LIB_UTILITIES_EXPORT Array<OneD, DataType>& GetPtr();

            LIB_UTILITIES_EXPORT const DataType* GetRawPtr() const;

            LIB_UTILITIES_EXPORT const Array<OneD, const DataType>& GetPtr() const;

            typedef DataType* iterator;
            LIB_UTILITIES_EXPORT iterator begin();
            LIB_UTILITIES_EXPORT iterator end();

            typedef const DataType* const_iterator;
            LIB_UTILITIES_EXPORT const_iterator begin() const;
            LIB_UTILITIES_EXPORT const_iterator end() const;

            /// \brief Returns i^{th} element.
            /// \param i The element to return.
            /// \pre i < dim
            /// \return A reference to the i^{th} element.
            ///
            /// Retrieves the i^{th} element.  Since it returns a reference you may
            /// assign a new value (i.e., p(2) = 3.2;)
            ///
            /// This operator performs range checking.
            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::reference operator()(unsigned int i);

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::reference operator[](unsigned int i);

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::reference x();

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::reference y();

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::reference z();

            LIB_UTILITIES_EXPORT void SetX(typename boost::call_traits<DataType>::const_reference val);

            LIB_UTILITIES_EXPORT void SetY(typename boost::call_traits<DataType>::const_reference val);

            LIB_UTILITIES_EXPORT void SetZ(typename boost::call_traits<DataType>::const_reference val);

            LIB_UTILITIES_EXPORT NekVector<DataType>& operator+=(const NekVector<DataType>& rhs);

            LIB_UTILITIES_EXPORT NekVector<DataType>& operator-=(const NekVector<DataType>& rhs);

            LIB_UTILITIES_EXPORT NekVector<DataType>& operator*=(typename boost::call_traits<DataType>::const_reference rhs);

            LIB_UTILITIES_EXPORT NekVector<DataType>& operator/=(typename boost::call_traits<DataType>::const_reference rhs);

            LIB_UTILITIES_EXPORT void Normalize();

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const;

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference operator[](unsigned int i) const;

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference x() const;

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference y() const;

            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference z() const;

            LIB_UTILITIES_EXPORT NekVector<DataType> operator-() const;

            LIB_UTILITIES_EXPORT DataType Magnitude() const;
            LIB_UTILITIES_EXPORT DataType Dot(const NekVector<DataType>& rhs) const;

            LIB_UTILITIES_EXPORT NekVector<DataType> Cross(const NekVector<DataType>& rhs) const;

            LIB_UTILITIES_EXPORT std::string AsString() const;

            // Norms
            LIB_UTILITIES_EXPORT DataType L1Norm() const;
            LIB_UTILITIES_EXPORT DataType L2Norm() const;
            LIB_UTILITIES_EXPORT DataType InfinityNorm() const;


            LIB_UTILITIES_EXPORT PointerWrapper GetWrapperType() const;

        protected:

            LIB_UTILITIES_EXPORT Array<OneD, DataType>& GetData();
            LIB_UTILITIES_EXPORT void SetSize(unsigned int s);
            LIB_UTILITIES_EXPORT void SetWrapperType(PointerWrapper p);
            LIB_UTILITIES_EXPORT void SetData(const Array<OneD, DataType>& newData);
            LIB_UTILITIES_EXPORT void Resize(unsigned int newSize);

        private:
            // Prevents accidental use of wrapped mode around ConstArrays.
            NekVector(const Array<OneD, const DataType>& ptr, PointerWrapper h);

            unsigned int m_size;
            Array<OneD, DataType> m_data;
            PointerWrapper m_wrapperType;
    };

    template<typename DataType>
    LIB_UTILITIES_EXPORT void Add(NekVector<DataType>& result,
           const NekVector<DataType>& lhs,
           const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void AddNegatedLhs(NekVector<DataType>& result,
           const NekVector<DataType>& lhs,
           const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void AddEqual(NekVector<DataType>& result,
           const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void AddEqualNegatedLhs(NekVector<DataType>& result,
           const NekVector<DataType>& rhs);

    template<typename LhsDataType,
             typename RhsDataType>
    LIB_UTILITIES_EXPORT NekVector<LhsDataType> Add(const NekVector<LhsDataType>& lhs,
                               const NekVector<RhsDataType>& rhs);



    template<typename ResultDataType, typename InputDataType>
    LIB_UTILITIES_EXPORT void Subtract(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekVector<InputDataType>& rhs);

    template<typename ResultDataType, typename InputDataType>
    LIB_UTILITIES_EXPORT void SubtractNegatedLhs(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekVector<InputDataType>& rhs);

    template<typename ResultDataType, typename InputDataType>
    LIB_UTILITIES_EXPORT void SubtractEqual(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& rhs);

    template<typename ResultDataType, typename InputDataType>
    LIB_UTILITIES_EXPORT void SubtractEqualNegatedLhs(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& rhs);

    template<typename DataType>
    NekVector<DataType>
    LIB_UTILITIES_EXPORT Subtract(const NekVector<DataType>& lhs,
                const NekVector<DataType>& rhs);





    template<typename ResultDataType, typename InputDataType>
    void LIB_UTILITIES_EXPORT Divide(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekDouble& rhs);

    template<typename ResultDataType>
    void LIB_UTILITIES_EXPORT DivideEqual(NekVector<ResultDataType>& result,
           const NekDouble& rhs);

    template<typename DataType>
    NekVector<DataType>
    LIB_UTILITIES_EXPORT Divide(const NekVector<DataType>& lhs,
                const NekDouble& rhs);


    template<typename ResultDataType, typename InputDataType>
    void LIB_UTILITIES_EXPORT Multiply(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekVector<InputDataType>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void LIB_UTILITIES_EXPORT MultiplyEqual(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& rhs);

    template<typename DataType, typename InputDataType>
    NekVector<DataType>
    LIB_UTILITIES_EXPORT Multiply(const NekVector<DataType>& lhs,
                const NekVector<InputDataType>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void LIB_UTILITIES_EXPORT Multiply(NekVector<ResultDataType>& result,
           const NekVector<InputDataType>& lhs,
           const NekDouble& rhs);

    template<typename ResultDataType>
    void LIB_UTILITIES_EXPORT MultiplyEqual(NekVector<ResultDataType>& result,
           const NekDouble& rhs);

    template<typename DataType>
    NekVector<DataType>
    LIB_UTILITIES_EXPORT Multiply(const NekVector<DataType>& lhs,
                const NekDouble& rhs);

    template<typename ResultDataType, typename InputDataType>
    void LIB_UTILITIES_EXPORT Multiply(NekVector<ResultDataType>& result,
                  const NekDouble& lhs,
                  const NekVector<InputDataType>& rhs);

    template<typename ResultDataType, typename InputDataType>
    void LIB_UTILITIES_EXPORT MultiplyInvertedLhs(NekVector<ResultDataType>& result,
                  const NekDouble& lhs,
                  const NekVector<InputDataType>& rhs);

    template<typename DataType>
    NekVector<DataType>
    LIB_UTILITIES_EXPORT Multiply(const DataType& lhs,
                  const NekVector<DataType>& rhs);

    template<typename DataType>
    NekVector<DataType> operator*(
        const NekVector<DataType> &lhs, const NekDouble &rhs)
    {
        return Multiply(lhs, rhs);
    }

    template<typename DataType>
    NekVector<DataType> operator*(
        const NekDouble &lhs, const NekVector<DataType>& rhs)
    {
        return Multiply(lhs, rhs);
    }

    template<typename DataType>
    NekVector<DataType> operator*(
        const NekVector<DataType> &lhs, const NekVector<DataType> &rhs)
    {
        return Multiply(lhs, rhs);
    }

    template<typename DataType>
    NekVector<DataType> operator/(
        const NekVector<DataType> &lhs, const NekDouble &rhs)
    {
        return Divide(lhs, rhs);
    }

    template<typename DataType>
    NekVector<DataType> operator+(
        const NekVector<DataType> &lhs, const NekVector<DataType> &rhs)
    {
        return Add(lhs, rhs);
    }

    template<typename DataType>
    NekVector<DataType> operator-(
        const NekVector<DataType> &lhs, const NekVector<DataType> &rhs)
    {
        return Subtract(lhs, rhs);
    }

    template<typename DataType>
    LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT NekVector<DataType> createVectorFromPoints(const NekPoint<DataType>& source, const NekPoint<DataType>& dest);

    template<typename DataType>
    LIB_UTILITIES_EXPORT NekPoint<DataType> findPointAlongVector(const NekVector<DataType>& lhs, const DataType& t);

    template<typename DataType>
    LIB_UTILITIES_EXPORT bool operator==(const NekVector<DataType>& lhs, const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT bool operator!=(const NekVector<DataType>& lhs, const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT DataType Magnitude(const NekVector<DataType>& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT DataType Dot(const NekVector<DataType>& lhs,
                 const NekVector<DataType>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT std::vector<DataType> FromString(const std::string& str);

    /// \todo Do the Norms with Blas where applicable.
    template<typename DataType>
    LIB_UTILITIES_EXPORT DataType L1Norm(const NekVector<DataType>& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT DataType L2Norm(const NekVector<DataType>& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT DataType InfinityNorm(const NekVector<DataType>& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT NekVector<DataType> Negate(const NekVector<DataType>& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void NegateInPlace(NekVector<DataType>& v);

    LIB_UTILITIES_EXPORT void NegateInPlace(NekDouble& v);
    LIB_UTILITIES_EXPORT void InvertInPlace(NekDouble& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void Normalize(NekVector<DataType>& v);

    template<typename DataType>
    LIB_UTILITIES_EXPORT NekVector<DataType> Cross(const NekVector<DataType>& lhs,
                                          const NekVector<DataType>& rhs);
    template<typename DataType>
    LIB_UTILITIES_EXPORT std::string AsString(const NekVector<DataType>& v);

}

#endif // NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

