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
// Description: Generic N-Dimensional Vector.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP
#define NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/LinearAlgebra/NekPoint.hpp>

#include <LibUtilities/LinearAlgebra/NekVectorMetadata.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorConstantSized.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorVariableSized.hpp>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>

#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>

#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>

#include <functional>
#include <algorithm>
#include <math.h>

#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_array.hpp>
#include <boost/typeof/typeof.hpp>

#include  BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

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

            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                template<typename L, typename Op, typename R>
                NekVector(const expt::Node<L, Op, R>& rhs) :
                    m_size(MatrixSize<expt::Node<L, Op, R>, typename expt::Node<L, Op, R>::Indices, 0>::template GetRequiredSize(rhs.GetData()).template get<0>()),
                    m_data(m_size),
                    m_wrapperType(eCopy)
                {
                    expt::ExpressionEvaluator::EvaluateWithoutAliasingCheck(rhs, *this);
                }

                template<typename L, typename Op, typename R>
                NekVector<DataType>& operator=(const expt::Node<L, Op, R>& rhs)
                {
                    boost::tuple<unsigned int, unsigned int, unsigned int> sizes =
                            MatrixSize<expt::Node<L, Op, R>, typename expt::Node<L, Op, R>::Indices, 0>::GetRequiredSize(rhs.GetData());
                    unsigned int newRows = sizes.get<0>();

                    this->SetSize(newRows);
                    if( this->GetWrapperType() == eCopy )
                    {
                        if( this->GetData().num_elements() < newRows )
                        {
                            this->SetData(Array<OneD, DataType>(newRows));
                        }
                    }
                    else
                    {
                        ASSERTL0(this->GetData().num_elements() >= newRows,
                                 "Attempting to store too many elements in a wrapped vector.");
                    }

                    expt::ExpressionEvaluator::Evaluate(rhs, *this);
                    return *this;
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

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

            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            expt::Node<expt::Node<NekVector<DataType> >, expt::NegateOp, void > operator-() const
            {
                expt::Node<NekVector<DataType> > leafNode(*this);
                return expt::Node<expt::Node<NekVector<DataType> >, expt::NegateOp, void >(leafNode);
            }
            #else
            LIB_UTILITIES_EXPORT NekVector<DataType> operator-() const;
            #endif

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

    BOOST_TYPEOF_REGISTER_TEMPLATE(NekVector, 1);

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

    GENERATE_MULTIPLICATION_OPERATOR(NekVector, 1, NekDouble, 0);
    GENERATE_MULTIPLICATION_OPERATOR(NekDouble, 0, NekVector, 1);
    GENERATE_MULTIPLICATION_OPERATOR(NekVector, 1, NekVector, 1);
    
    GENERATE_DIVISION_OPERATOR(NekVector, 1, NekDouble, 0);
    GENERATE_ADDITION_OPERATOR(NekVector, 1, NekVector, 1);
    GENERATE_SUBTRACTION_OPERATOR(NekVector, 1, NekVector, 1);


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

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
namespace expt
{
    template<typename DataType, typename L, typename Op, typename R>
    DataType Dot(const Nektar::NekVector<DataType>& lhs,
                 const expt::Node<L, Op, R>& expr)
    {
        // Evalaute the expression first, expression templates don't chain past
        // this point since the return value is a scalar.
        typename expt::Node<L, Op, R>::ResultType rhs = expt::ExpressionEvaluator::Evaluate(expr);
        return Dot(lhs, rhs);
    }

    template<typename DataType>
    struct IsAlias<Nektar::NekVector<DataType>, Nektar::NekVector<DataType> >
    {
        static bool Apply(const Nektar::NekVector<DataType>& lhs, const Nektar::NekVector<DataType>& rhs)
        {
            return lhs.GetPtr().Overlaps(rhs.GetPtr());
        }
    };

    template<typename DataType, typename NodeType, typename Indices, unsigned int StartIndex>
    struct CreateFromTree<Nektar::NekVector<DataType>, NodeType, Indices, StartIndex>
    {
        template<typename ArgVectorType>
        static Nektar::NekVector<DataType> Apply(const ArgVectorType& tree)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> sizes =
                Nektar::MatrixSize<NodeType, Indices, StartIndex>::GetRequiredSize(tree);

            unsigned int rows = sizes.get<0>();
            return Nektar::NekVector<DataType>(rows);
        }
    };

    // Override default expression handling for addition/subtraction of vectors.
    // Optimal execution is obtained by loop unrolling.

    template<typename NodeType, typename enabled = void>
    struct NodeCanUnroll : public boost::false_type {};

    template<typename Type>
    struct NodeCanUnroll<expt::Node<Type, void, void>,
        typename boost::enable_if
        <
            Nektar::IsVector<typename expt::Node<Type, void, void>::ResultType>
        >::type > : public boost::true_type
    {
    };
        
    template<typename LhsType, typename OpType, typename RhsType>
    struct NodeCanUnroll<expt::Node<LhsType, OpType, RhsType>,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                Nektar::IsVector<typename LhsType::ResultType>,
                Nektar::IsVector<typename RhsType::ResultType>,
                NodeCanUnroll<LhsType>,
                NodeCanUnroll<RhsType>,
                //boost::mpl::or_
                //<
                    boost::is_same<OpType, expt::AddOp>
                //    boost::is_same<OpType, expt::SubtractOp>
                //>
        > >::type >: public boost::true_type
    {
    };

    template<typename NodeType, typename IndicesType, unsigned int index>
    struct Accumulate;

    template<typename LhsType, typename IndicesType, unsigned int index>
    struct Accumulate<expt::Node<LhsType, void, void>, IndicesType, index>
    {
        static const unsigned int MappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;

        template<typename ResultType, typename ArgumentVectorType>
        static void Execute(ResultType& accumulator, const ArgumentVectorType& args, unsigned int i)
        {
            accumulator = boost::fusion::at_c<MappedIndex>(args)[i];
        }
    };

    template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index>
    struct Accumulate<expt::Node<LhsType, Op, RhsType>, IndicesType, index>
    {
        static const int rhsNodeIndex = index + LhsType::TotalCount;

        template<typename ResultType, typename ArgumentVectorType>
        static void Execute(ResultType& accumulator, const ArgumentVectorType& args, unsigned int i)
        {
            Accumulate<LhsType, IndicesType, index>::Execute(accumulator, args, i);
            ResultType rhs;
            Accumulate<RhsType, IndicesType, rhsNodeIndex>::Execute(rhs, args, i);
            Op::OpEqual(accumulator, rhs);
        }
    };



    template<typename IndicesType, unsigned int startIndex, unsigned int endIndex, typename enabled=void>
    struct Unroll;

    #ifndef NEKTAR_NEKVECTOR_MAX_UNROLL_ARGS
    #define NEKTAR_NEKVECTOR_MAX_UNROLL_ARGS 10
    #endif

    #define NEKTAR_NEKVECTOR_UNROLL_GENERATE_INDEX(z, n, IndexName) \
        static const unsigned int BOOST_PP_CAT(IndexName, n) = boost::mpl::at_c<IndicesType, startIndex+n>::type::value;

    #define NEKTAR_NEKVECTOR_UNROLL_GENERATE_VARIABLE(z, n, VariableName) \
        BOOST_AUTO(BOOST_PP_CAT(VariableName, n), boost::fusion::at_c<BOOST_PP_CAT(index, n)>(args).GetRawPtr());

    #define NEKTAR_NEKVECTOR_UNROLL_GENERATE_VARIABLE_NAME_IN_ADDITION_SEQUENCE(z, n, VariableName) \
        + BOOST_PP_CAT(VariableName, n)[i]

    #define NEKTAR_NEKVECTOR_UNROLL_IMPL(z, n, ClassName) \
    template<typename IndicesType, unsigned int startIndex, unsigned int endIndex> \
    struct ClassName<IndicesType, startIndex, endIndex, \
        typename boost::enable_if_c \
        < \
            endIndex-startIndex == BOOST_PP_CAT(n, u) \
        >::type> \
    { \
        BOOST_PP_REPEAT_FROM_TO(0, n, NEKTAR_NEKVECTOR_UNROLL_GENERATE_INDEX, index)\
        \
        template<typename AccumulatorType, typename ArgumentVectorType> \
        static inline void Execute(AccumulatorType& accumulator, const ArgumentVectorType& args) \
        { \
            BOOST_AUTO(a, accumulator.GetRawPtr()); \
            BOOST_PP_REPEAT_FROM_TO(0, n, NEKTAR_NEKVECTOR_UNROLL_GENERATE_VARIABLE, t) \
            \
            const unsigned int r = accumulator.GetRows(); \
            for(unsigned int i = 0; i < r; ++i) \
            { \
                a[i] = t0[i] \
                BOOST_PP_REPEAT_FROM_TO(1, n, NEKTAR_NEKVECTOR_UNROLL_GENERATE_VARIABLE_NAME_IN_ADDITION_SEQUENCE, t); \
            } \
        } \
    };

    BOOST_PP_REPEAT_FROM_TO(2, NEKTAR_NEKVECTOR_MAX_UNROLL_ARGS, NEKTAR_NEKVECTOR_UNROLL_IMPL, Unroll);

    //template<typename IndicesType, unsigned int startIndex, unsigned int endIndex>
    //struct Unroll<IndicesType, startIndex, endIndex,
    //    typename boost::enable_if_c
    //    <
    //        endIndex-startIndex == 2u
    //    >::type>
    //{
    //    static const unsigned int index0 = boost::mpl::at_c<IndicesType, startIndex>::type::value;
    //    static const unsigned int index1 = boost::mpl::at_c<IndicesType, startIndex+1>::type::value;
    //    ;
    //    template<typename AccumulatorType, typename ArgumentVectorType>
    //    static inline void Execute(AccumulatorType& accumulator, const ArgumentVectorType& args)
    //    {
    //        BOOST_AUTO(a, accumulator.GetRawPtr());
    //        BOOST_AUTO(t0, boost::fusion::at_c<index0>(args).GetRawPtr());
    //        BOOST_AUTO(t1, boost::fusion::at_c<index1>(args).GetRawPtr());

    //        const unsigned int r = accumulator.GetRows();
    //        for(unsigned int i = 0; i < r; ++i)
    //        {
    //            accumulator[i] = t0[i] + t1[i];
    //        }
    //    }
    //};

    //template<typename IndicesType, unsigned int startIndex, unsigned int endIndex>
    //struct Unroll<IndicesType, startIndex, endIndex,
    //    typename boost::enable_if_c
    //    <
    //        endIndex-startIndex == 3u
    //    >::type>
    //{
    //    static const unsigned int index0 = boost::mpl::at_c<IndicesType, startIndex>::type::value;
    //    static const unsigned int index1 = boost::mpl::at_c<IndicesType, startIndex+1>::type::value;
    //    static const unsigned int index2 = boost::mpl::at_c<IndicesType, startIndex+2>::type::value;

    //    template<typename AccumulatorType, typename ArgumentVectorType>
    //    static inline void Execute(AccumulatorType& accumulator, const ArgumentVectorType& args)
    //    {
    //        BOOST_AUTO(a, accumulator.GetRawPtr());
    //        BOOST_AUTO(t0, boost::fusion::at_c<index0>(args).GetRawPtr());
    //        BOOST_AUTO(t1, boost::fusion::at_c<index1>(args).GetRawPtr());
    //        BOOST_AUTO(t2, boost::fusion::at_c<index2>(args).GetRawPtr());

    //        const unsigned int r = accumulator.GetRows();
    //        for(unsigned int i = 0; i < r; ++i)
    //        {
    //            accumulator[i] = t0[i] + t1[i] + t2[i];
    //        }
    //    }
    //};

    //template<typename IndicesType, unsigned int startIndex, unsigned int endIndex>
    //struct Unroll<IndicesType, startIndex, endIndex,
    //    typename boost::enable_if_c
    //    <
    //        endIndex-startIndex == 4u
    //    >::type>
    //{
    //    static const unsigned int index0 = boost::mpl::at_c<IndicesType, startIndex>::type::value;
    //    static const unsigned int index1 = boost::mpl::at_c<IndicesType, startIndex+1>::type::value;
    //    static const unsigned int index2 = boost::mpl::at_c<IndicesType, startIndex+2>::type::value;
    //    static const unsigned int index3 = boost::mpl::at_c<IndicesType, startIndex+3>::type::value;

    //    template<typename AccumulatorType, typename ArgumentVectorType>
    //    static inline void Execute(AccumulatorType& accumulator, const ArgumentVectorType& args)
    //    {
    //        BOOST_AUTO(a, accumulator.GetRawPtr());
    //        BOOST_AUTO(t0, boost::fusion::at_c<index0>(args).GetRawPtr());
    //        BOOST_AUTO(t1, boost::fusion::at_c<index1>(args).GetRawPtr());
    //        BOOST_AUTO(t2, boost::fusion::at_c<index2>(args).GetRawPtr());
    //        BOOST_AUTO(t3, boost::fusion::at_c<index3>(args).GetRawPtr());

    //        const unsigned int r = accumulator.GetRows();
    //        for(unsigned int i = 0; i < r; ++i)
    //        {
    //            accumulator[i] = t0[i] + t1[i] + t2[i] + t3[i];
    //        }
    //    }
    //};

        // Conditions
        // Lhs and Rhs must result in a vector.
        // Op must be Plus or Minus
        // Must apply recursively.
    template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<LhsType, Op, RhsType, IndicesType, index,
        typename boost::enable_if
        <
            NodeCanUnroll<expt::Node<LhsType, Op, RhsType> >
        >::type> : public boost::true_type 
    {
        static const int endIndex = index + LhsType::TotalCount + RhsType::TotalCount;

        template<typename ResultType, typename ArgumentVectorType>
        static inline void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            Unroll<IndicesType, index, endIndex>::Execute(accumulator, args);
        }
    };
}
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES


#endif // NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

