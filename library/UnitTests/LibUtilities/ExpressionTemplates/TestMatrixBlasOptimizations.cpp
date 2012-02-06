///////////////////////////////////////////////////////////////////////////////
//
// File: TestMatrixBlasOptimizations.cpp
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

#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <UnitTests/CountedObject.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/timer.hpp>

using namespace expt;

namespace Nektar
{
    BOOST_AUTO_TEST_CASE(TestTestBinaryNode)
    {
        typedef Node<NekMatrix<double> > A;
        typedef Node<NekVector<double> > B;
        typedef Node<double> C;

        typedef Node<A, MultiplyOp, B> Test1;
        typedef Node<A, MultiplyOp, C> Test2;
        typedef Node<Node<A, MultiplyOp, A>, MultiplyOp, B> Test3;

        BOOST_MPL_ASSERT(( TestBinaryNode<Test1, CanGetRawPtr, MultiplyOp, IsVector >));
        BOOST_MPL_ASSERT(( boost::mpl::not_<TestBinaryNode<Test1, CanGetRawPtr, AddOp, IsVector > >));
        BOOST_MPL_ASSERT(( boost::mpl::not_<TestBinaryNode<Test1, CanGetRawPtr, MultiplyOp, impl::IsDouble > >));
        BOOST_MPL_ASSERT(( boost::mpl::not_<TestBinaryNode<Test1, impl::IsDouble, MultiplyOp, CanGetRawPtr > >));
        BOOST_MPL_ASSERT(( TestBinaryNode<Test3, CanGetRawPtr, MultiplyOp, IsVector >));
        BOOST_MPL_ASSERT(( TestBinaryNode<Test2, CanGetRawPtr, MultiplyOp, impl::IsDouble >));
        
    }

    BOOST_AUTO_TEST_CASE(TestTest3ArgumentAssociativeNode)
    {
        typedef Node<NekMatrix<double> > A;
        typedef Node<NekVector<double> > B;
        typedef Node<double> C;

        typedef Node<Node<A, MultiplyOp, A>, MultiplyOp, C> Test1;

        BOOST_MPL_ASSERT(( Test3ArgumentAssociativeNode<Test1, CanGetRawPtr, MultiplyOp, CanGetRawPtr, MultiplyOp, impl::IsDouble> ));
        BOOST_MPL_ASSERT(( TestBinaryNode<Test1, CanGetRawPtr, MultiplyOp, impl::IsDouble> ));

        typedef Node<A, MultiplyOp, B> Test2;
        BOOST_MPL_ASSERT(( boost::mpl::and_<TestBinaryNode<Test2, CanGetRawPtr, MultiplyOp, IsVector >, boost::mpl::not_<Test3ArgumentAssociativeNode<Test2, CanGetRawPtr, MultiplyOp, impl::IsDouble, MultiplyOp, IsVector > > > ));
    }

    BOOST_AUTO_TEST_CASE(TestDgemmAlphaABDetection)
    {
        typedef Node<NekMatrix<double> > M;
        typedef Node<double> D;
        typedef Node<M, MultiplyOp, M> Test1;

        BOOST_MPL_ASSERT(( impl::AlphaABParameterAccess<Test1, Test1::Indices, 0> ));

        typedef Node< Node<M, MultiplyOp, M>, MultiplyOp, D> Test2;
        typedef Node< Node<D, MultiplyOp, M>, MultiplyOp, M> Test3;
        typedef Node< Node<M, MultiplyOp, D>, MultiplyOp, M> Test4;

        BOOST_MPL_ASSERT(( impl::AlphaABParameterAccess<Test2, Test2::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::AlphaABParameterAccess<Test3, Test3::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::AlphaABParameterAccess<Test4, Test4::Indices, 0> ));

        typedef Node< Node<M, AddOp, M>, MultiplyOp, D> Test5;
        typedef Node< Node<D, MultiplyOp, M>, MultiplyOp, D> Test7;
        typedef Node< Node<M, MultiplyOp, D>, MultiplyOp, D> Test8;

        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaABParameterAccess<Test5, Test5::Indices, 0> > ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaABParameterAccess<Test7, Test7::Indices, 0> > ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaABParameterAccess<Test8, Test8::Indices, 0> > ));

        typedef Node<Test1, NegateOp> Test9;
        typedef Node<Test2, NegateOp> Test10;
        typedef Node<Test5, NegateOp> Test11;
        typedef Node<Test7, NegateOp> Test12;

        BOOST_MPL_ASSERT(( impl::AlphaABParameterAccess<Test9, Test9::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::AlphaABParameterAccess<Test10, Test10::Indices, 0> ));

        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaABParameterAccess<Test11, Test11::Indices, 0> > ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaABParameterAccess<Test12, Test12::Indices, 0> > ));

    }

    BOOST_AUTO_TEST_CASE(TestDgemvAlphaAXDetection)
    {
        typedef Node<NekMatrix<double> > M;
        typedef Node<double> D;
        typedef Node<NekVector<double> > V;
        typedef Node<M, MultiplyOp, V> Test1;

        BOOST_MPL_ASSERT(( impl::AlphaAXParameterAccess<Test1, Test1::Indices, 0> ));

        typedef Node< Node<M, MultiplyOp, V>, MultiplyOp, D> Test2;
        typedef Node< Node<D, MultiplyOp, M>, MultiplyOp, V> Test3;
        typedef Node< Node<M, MultiplyOp, D>, MultiplyOp, V> Test4;

        BOOST_MPL_ASSERT(( impl::AlphaAXParameterAccess<Test2, Test2::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::AlphaAXParameterAccess<Test3, Test3::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::AlphaAXParameterAccess<Test4, Test4::Indices, 0> ));

        typedef Node< Node<M, AddOp, M>, MultiplyOp, D> Test5;
        typedef Node< Node<D, MultiplyOp, M>, MultiplyOp, D> Test7;
        typedef Node< Node<M, MultiplyOp, D>, MultiplyOp, D> Test8;

        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaAXParameterAccess<Test5, Test5::Indices, 0> > ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaAXParameterAccess<Test7, Test7::Indices, 0> > ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaAXParameterAccess<Test8, Test8::Indices, 0> > ));

        typedef Node<Test1, NegateOp> Test9;
        typedef Node<Test2, NegateOp> Test10;
        typedef Node<Test5, NegateOp> Test11;
        typedef Node<Test7, NegateOp> Test12;

        BOOST_MPL_ASSERT(( impl::AlphaAXParameterAccess<Test9, Test9::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::AlphaAXParameterAccess<Test10, Test10::Indices, 0> ));

        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaAXParameterAccess<Test11, Test11::Indices, 0> > ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::AlphaAXParameterAccess<Test12, Test12::Indices, 0> > ));
    }

    BOOST_AUTO_TEST_CASE(TestDgemmBetaCDetection)
    {
        typedef Node<NekMatrix<double> > M;
        typedef Node<NekVector<double> > V;
        typedef Node<double> D;

        typedef M Test1;
        typedef Node<M, MultiplyOp, D> Test2;
        typedef Node<D, MultiplyOp, M> Test3;
        typedef Node<V, MultiplyOp, D> Test4;

        BOOST_MPL_ASSERT(( impl::BetaCParameterAccess<Test1, Test1::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::BetaCParameterAccess<Test2, Test2::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::BetaCParameterAccess<Test3, Test3::Indices, 0> ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::BetaCParameterAccess<Test4, Test4::Indices, 0> > ));
    }

    BOOST_AUTO_TEST_CASE(TestDgemvBetaYDetection)
    {
        typedef Node<NekMatrix<double> > M;
        typedef Node<NekVector<double> > V;
        typedef Node<double> D;

        typedef M Test1;
        typedef Node<M, MultiplyOp, D> Test2;
        typedef Node<D, MultiplyOp, M> Test3;
        typedef Node<V, MultiplyOp, D> Test4;

        BOOST_MPL_ASSERT(( impl::BetaCParameterAccess<Test1, Test1::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::BetaCParameterAccess<Test2, Test2::Indices, 0> ));
        BOOST_MPL_ASSERT(( impl::BetaCParameterAccess<Test3, Test3::Indices, 0> ));
        BOOST_MPL_ASSERT(( boost::mpl::not_<impl::BetaCParameterAccess<Test4, Test4::Indices, 0> > ));
    }

    BOOST_AUTO_TEST_CASE(TestFullDgemmDetection)
    {
        typedef Node<NekMatrix<double> > M;
        typedef Node<NekVector<double> > V;
        typedef Node<double> D;

        typedef Node<M, MultiplyOp, M> Lhs1;
        typedef Node< Node<M, MultiplyOp, M>, MultiplyOp, D> Lhs2;
        typedef Node< Node<D, MultiplyOp, M>, MultiplyOp, M> Lhs3;
        typedef Node< Node<M, MultiplyOp, D>, MultiplyOp, M> Lhs4;
        typedef Node<Lhs1, NegateOp> Lhs5;
        typedef Node<Lhs2, NegateOp> Lhs6;
        typedef Node<Lhs3, NegateOp> Lhs7;
        typedef Node<Lhs4, NegateOp> Lhs8;

        typedef M Rhs1;
        typedef Node<M, MultiplyOp, D> Rhs2;
        typedef Node<D, MultiplyOp, M> Rhs3;
        typedef Node<Rhs1, NegateOp> Rhs4;
        typedef Node<Rhs2, NegateOp> Rhs5;
        typedef Node<Rhs3, NegateOp> Rhs6;

        typedef Node<Lhs1, AddOp, Rhs1> Test1;

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, AddOp, Rhs1, Node<Lhs1, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, AddOp, Rhs2, Node<Lhs1, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, AddOp, Rhs3, Node<Lhs1, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, AddOp, Rhs4, Node<Lhs1, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, AddOp, Rhs5, Node<Lhs1, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, AddOp, Rhs6, Node<Lhs1, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, AddOp, Rhs1, Node<Lhs2, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, AddOp, Rhs2, Node<Lhs2, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, AddOp, Rhs3, Node<Lhs2, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, AddOp, Rhs4, Node<Lhs2, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, AddOp, Rhs5, Node<Lhs2, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, AddOp, Rhs6, Node<Lhs2, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, AddOp, Rhs1, Node<Lhs3, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, AddOp, Rhs2, Node<Lhs3, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, AddOp, Rhs3, Node<Lhs3, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, AddOp, Rhs4, Node<Lhs3, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, AddOp, Rhs5, Node<Lhs3, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, AddOp, Rhs6, Node<Lhs3, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, AddOp, Rhs1, Node<Lhs4, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, AddOp, Rhs2, Node<Lhs4, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, AddOp, Rhs3, Node<Lhs4, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, AddOp, Rhs4, Node<Lhs4, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, AddOp, Rhs5, Node<Lhs4, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, AddOp, Rhs6, Node<Lhs4, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, AddOp, Rhs1, Node<Lhs5, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, AddOp, Rhs2, Node<Lhs5, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, AddOp, Rhs3, Node<Lhs5, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, AddOp, Rhs4, Node<Lhs5, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, AddOp, Rhs5, Node<Lhs5, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, AddOp, Rhs6, Node<Lhs5, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, AddOp, Rhs1, Node<Lhs6, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, AddOp, Rhs2, Node<Lhs6, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, AddOp, Rhs3, Node<Lhs6, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, AddOp, Rhs4, Node<Lhs6, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, AddOp, Rhs5, Node<Lhs6, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, AddOp, Rhs6, Node<Lhs6, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs1, Node<Lhs7, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs2, Node<Lhs7, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs3, Node<Lhs7, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs4, Node<Lhs7, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs5, Node<Lhs7, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs6, Node<Lhs7, AddOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs1, Node<Lhs8, AddOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs2, Node<Lhs8, AddOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs3, Node<Lhs8, AddOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs4, Node<Lhs8, AddOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs5, Node<Lhs8, AddOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, AddOp, Rhs6, Node<Lhs8, AddOp, Rhs6>::Indices, 0> ));




        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, SubtractOp, Rhs1, Node<Lhs1, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, SubtractOp, Rhs2, Node<Lhs1, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, SubtractOp, Rhs3, Node<Lhs1, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, SubtractOp, Rhs4, Node<Lhs1, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, SubtractOp, Rhs5, Node<Lhs1, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs1, SubtractOp, Rhs6, Node<Lhs1, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, SubtractOp, Rhs1, Node<Lhs2, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, SubtractOp, Rhs2, Node<Lhs2, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, SubtractOp, Rhs3, Node<Lhs2, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, SubtractOp, Rhs4, Node<Lhs2, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, SubtractOp, Rhs5, Node<Lhs2, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs2, SubtractOp, Rhs6, Node<Lhs2, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, SubtractOp, Rhs1, Node<Lhs3, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, SubtractOp, Rhs2, Node<Lhs3, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, SubtractOp, Rhs3, Node<Lhs3, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, SubtractOp, Rhs4, Node<Lhs3, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, SubtractOp, Rhs5, Node<Lhs3, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs3, SubtractOp, Rhs6, Node<Lhs3, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, SubtractOp, Rhs1, Node<Lhs4, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, SubtractOp, Rhs2, Node<Lhs4, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, SubtractOp, Rhs3, Node<Lhs4, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, SubtractOp, Rhs4, Node<Lhs4, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, SubtractOp, Rhs5, Node<Lhs4, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs4, SubtractOp, Rhs6, Node<Lhs4, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, SubtractOp, Rhs1, Node<Lhs5, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, SubtractOp, Rhs2, Node<Lhs5, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, SubtractOp, Rhs3, Node<Lhs5, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, SubtractOp, Rhs4, Node<Lhs5, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, SubtractOp, Rhs5, Node<Lhs5, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs5, SubtractOp, Rhs6, Node<Lhs5, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, SubtractOp, Rhs1, Node<Lhs6, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, SubtractOp, Rhs2, Node<Lhs6, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, SubtractOp, Rhs3, Node<Lhs6, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, SubtractOp, Rhs4, Node<Lhs6, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, SubtractOp, Rhs5, Node<Lhs6, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs6, SubtractOp, Rhs6, Node<Lhs6, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs1, Node<Lhs7, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs2, Node<Lhs7, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs3, Node<Lhs7, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs4, Node<Lhs7, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs5, Node<Lhs7, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs6, Node<Lhs7, SubtractOp, Rhs6>::Indices, 0> ));

        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs1, Node<Lhs8, SubtractOp, Rhs1>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs2, Node<Lhs8, SubtractOp, Rhs2>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs3, Node<Lhs8, SubtractOp, Rhs3>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs4, Node<Lhs8, SubtractOp, Rhs4>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs5, Node<Lhs8, SubtractOp, Rhs5>::Indices, 0> ));
        BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Lhs7, SubtractOp, Rhs6, Node<Lhs8, SubtractOp, Rhs6>::Indices, 0> ));

    }


    BOOST_AUTO_TEST_CASE(TestAlphaABOptimization)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);

        NekMatrix<double> test1 = 2.0*a*b;
        NekMatrix<double> test2 = a*2.0*b;
        NekMatrix<double> test3 = a*b*2.0;

        double expected_result_buf[] = {38, 56, 54, 80};
        NekMatrix<double> expected_result(2,2,expected_result_buf);

        BOOST_CHECK_EQUAL(test1, expected_result);
        BOOST_CHECK_EQUAL(test2, expected_result);
        BOOST_CHECK_EQUAL(test3, expected_result);

        NekMatrix<double> test4 = 3.0*(-a)*b;
        NekMatrix<double> test5 = (-a)*3.0*b;
        NekMatrix<double> test6 = (-a)*b*3.0;

        double expected_result_buf1[] = {-57, -84, -81, -120};
        NekMatrix<double> expected_result1(2,2,expected_result_buf1);

        BOOST_CHECK_EQUAL(test4, expected_result1);
        BOOST_CHECK_EQUAL(test5, expected_result1);
        BOOST_CHECK_EQUAL(test6, expected_result1);

        NekMatrix<double> test7 = 4.0*(-a)*(-b);
        NekMatrix<double> test8 = (-a)*4.0*(-b);
        NekMatrix<double> test9 = (-a)*(-b)*4.0;

        double expected_result_buf2[] = {76, 112, 108, 160};
        NekMatrix<double> expected_result2(2,2,expected_result_buf2);

        BOOST_CHECK_EQUAL(test7, expected_result2);
        BOOST_CHECK_EQUAL(test8, expected_result2);
        BOOST_CHECK_EQUAL(test9, expected_result2);

        NekMatrix<double> test10 = (-(a*b))*5.0;
        NekMatrix<double> test11 = 5.0*(-(a*b));

        double expected_result_buf3[] = {-95, -140, -135, -200};
        NekMatrix<double> expected_result3(2,2,expected_result_buf3);

        BOOST_CHECK_EQUAL(test10, expected_result3);
        BOOST_CHECK_EQUAL(test11, expected_result3);
    }

    BOOST_AUTO_TEST_CASE(TestAlphaABBetaCOptimization)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);

        NekMatrix<double> test1 = 2.0*a*b + 3.0*c;
        NekMatrix<double> test2 = 2.0*a*b + c*3.0;
        NekMatrix<double> test3 = a*2.0*b + c*3.0;
        NekMatrix<double> test4 = a*b*2.0 + c*3.0;

        NekMatrix<double> test1i = 3.0*c + 2.0*a*b;
        NekMatrix<double> test2i = 3.0*c + 2.0*a*b;
        NekMatrix<double> test3i = 3.0*c + a*2.0*b;
        NekMatrix<double> test4i = 3.0*c + a*b*2.0;

        double expected_result_buf1[] = {62, 83, 84, 113};
        NekMatrix<double> expected_result1(2,2, expected_result_buf1);

        BOOST_CHECK_EQUAL(test1, expected_result1);
        BOOST_CHECK_EQUAL(test2, expected_result1);
        BOOST_CHECK_EQUAL(test3, expected_result1);
        BOOST_CHECK_EQUAL(test4, expected_result1);
        BOOST_CHECK_EQUAL(test1i, expected_result1);
        BOOST_CHECK_EQUAL(test2i, expected_result1);
        BOOST_CHECK_EQUAL(test3i, expected_result1);
        BOOST_CHECK_EQUAL(test4i, expected_result1);

        NekMatrix<double> test5 = 4.0*a*b - 7.0*c;
        NekMatrix<double> test6 = 4.0*a*b - c*7.0;
        NekMatrix<double> test7 = a*4.0*b - c*7.0;
        NekMatrix<double> test8 = a*b*4.0 - c*7.0;
        NekMatrix<double> test9 = a*b*4.0 + (-(c*7.0));

        double expected_result_buf2[] = {20, 49, 38, 83};
        NekMatrix<double> expected_result2(2,2, expected_result_buf2);

        BOOST_CHECK_EQUAL(test5, expected_result2);
        BOOST_CHECK_EQUAL(test6, expected_result2);
        BOOST_CHECK_EQUAL(test7, expected_result2);
        BOOST_CHECK_EQUAL(test8, expected_result2);
        BOOST_CHECK_EQUAL(test9, expected_result2);

        NekMatrix<double> test10 = -(a*b*4.0) - c*7.0;
        NekMatrix<double> test11 = -a*b*4.0 - c*7.0;
        NekMatrix<double> test12 = a*(-b)*4.0 - c*7.0;
        NekMatrix<double> test13 = a*b*(-4.0) - c*7.0;

        double expected_result_buf3[] = {-132, -175, -178, -237};
        NekMatrix<double> expected_result3(2,2, expected_result_buf3);

        BOOST_CHECK_EQUAL(test10, expected_result3);
        BOOST_CHECK_EQUAL(test11, expected_result3);
        BOOST_CHECK_EQUAL(test12, expected_result3);
        BOOST_CHECK_EQUAL(test13, expected_result3);

        NekMatrix<double> test14 = -a*b*4.0;
        NekMatrix<double> test15 = a*(-b)*4.0;
        NekMatrix<double> test16 = a*b*(-4.0);

        double expected_result_buf4[] = {-76, -112, -108, -160};
        NekMatrix<double> expected_result4(2,2, expected_result_buf4);

        BOOST_CHECK_EQUAL(test14, expected_result4);
        BOOST_CHECK_EQUAL(test15, expected_result4);
        BOOST_CHECK_EQUAL(test16, expected_result4);

    }

    BOOST_AUTO_TEST_CASE(TestIsAlphaABNodeNegateNodeInTree)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);

        typedef Node<NekMatrix<double> > L;
        typedef Node<double> D;
        typedef Node<L, MultiplyOp, L> MultiplyNode;

        typedef Node< Node<MultiplyNode, NegateOp>, MultiplyOp, D> TreeType;
        
        TreeType expression = -(a*b)*3.0;
        
        typedef RemoveUnecessaryTemporaries<TreeType>::TransformedNodeType OptimizedTreeType;
        BOOST_MPL_ASSERT(( boost::is_same<TreeType, OptimizedTreeType> ));

        NekMatrix<double> result = expression;

        NekMatrix<double> result1 = (-a)*b;
        double expected_result_buf[] = {-57, -84, -81, -120};
        NekMatrix<double> expected_result(2,2,expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }

    BOOST_AUTO_TEST_CASE(DgemmTest)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);
        
        NekMatrix<double> result = a*b + c;
        
        double expected_result_buf[] = {27, 37, 37, 51};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }
    

    BOOST_AUTO_TEST_CASE(TestNegatedRhs)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {5, 6, 7, 8};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);
        
        NekMatrix<double> result = a*b - c;

        double expected_result_buf[] = {15, 25, 21, 35};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }

    BOOST_AUTO_TEST_CASE(TestAB)
    {
        typedef Node<NekMatrix<double> > Leaf;
        typedef Node<Leaf, MultiplyOp, Leaf> MultiplyNode;

        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        
        boost::shared_ptr<NekMatrix<double> > ia(new NekMatrix<double>(2, 2, a_buf));
        boost::shared_ptr<NekMatrix<double> > ib(new NekMatrix<double>(2, 2, b_buf));
        
        MultiplyNode nodeTest = (*ia)*(*ib);
        NekMatrix<double> result1 = 2.0*(*ia)*(*ib);
        NekMatrix<double> result2 = (*ia)*3.0*(*ib);
        NekMatrix<double> result3 = (*ia)*(*ib)*4.0;

        double expectedResult1Buf[] = {38, 56, 54, 80};
        double expectedResult2Buf[] = {57, 84, 81, 120};
        double expectedResult3Buf[] = {76, 112, 108, 160};
        NekMatrix<double> expectedResult1(2,2,expectedResult1Buf);
        NekMatrix<double> expectedResult2(2,2,expectedResult2Buf);
        NekMatrix<double> expectedResult3(2,2,expectedResult3Buf);

        BOOST_CHECK_EQUAL(expectedResult1, result1);
        BOOST_CHECK_EQUAL(expectedResult2, result2);
        BOOST_CHECK_EQUAL(expectedResult3, result3);
       
    }


    BOOST_AUTO_TEST_CASE(DgemmScaledMatrixTest)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        boost::shared_ptr<NekMatrix<double> > ia(new NekMatrix<double>(2, 2, a_buf));
        boost::shared_ptr<NekMatrix<double> > ib(new NekMatrix<double>(2, 2, b_buf));
        boost::shared_ptr<NekMatrix<double> > ic(new NekMatrix<double>(2, 2, c_buf));
        
        NekMatrix<NekMatrix<double>, ScaledMatrixTag> a(2.0, ia);
        NekMatrix<NekMatrix<double>, ScaledMatrixTag> b(3.0, ib);
        NekMatrix<NekMatrix<double>, ScaledMatrixTag> c(-6.0, ic);
        
        NekMatrix<double> result = a*b + c;
        
        double expected_result_buf[] = {66, 114, 102, 174};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }

    //
    //BOOST_AUTO_TEST_CASE(DgemmDiagonalTest)
    //{
    //    double a_buf[] = {1, 2, 3, 4};
    //    double b_buf[] = {4, 5, 6, 7};
    //    double c_buf[] = {8, 9, 10, 11};
    //    
    //    NekMatrix<double> a(4, 4, a_buf, eDIAGONAL);
    //    NekMatrix<double> b(4, 4, b_buf, eDIAGONAL);
    //    NekMatrix<double> c(4, 4, c_buf, eDIAGONAL);
    //    
    //    NekMatrix<double> result = a*b + c;
    //    
    //    
    //    double expected_result_buf[] = {12, 0, 0, 0,
    //                                    0, 19, 0, 0,
    //                                    0, 0, 28, 0,
    //                                    0, 0, 0, 39};
    //    NekMatrix<double> expected_result(4, 4, expected_result_buf);
    //    BOOST_CHECK_EQUAL(expected_result, result);
    //}
    //
    BOOST_AUTO_TEST_CASE(TestMatrixMultiplyEqual01)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);
        
        NekMatrix<double> result = a*b*c;

        double expected_result_buf[] = {395, 584, 487, 720};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }

    BOOST_AUTO_TEST_CASE(TestMatrixMultiplyEqual02)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        double d_buf[] = {12, 13, 14, 15};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);
        NekMatrix<double> d(2, 2, d_buf);
        
        NekMatrix<double> result = a*b*c*d;

        double expected_result_buf[] = {11071, 16368, 12835, 18976};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }

    BOOST_AUTO_TEST_CASE(TestDgemmWithTransposedC)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);

        NekMatrix<double> test1 = 3.0*a*b + 4.0*c;
        NekMatrix<double> ct = Transpose(c);
        NekMatrix<double> test2 = 3.0*a*b + 4.0*ct;

        double expectedTest1Buf[] = {89, 120, 121, 164};
        NekMatrix<double> expectedTest1(2,2,expectedTest1Buf);

        double expectedTest2Buf[] = {89, 124, 117, 164};
        NekMatrix<double> expectedTest2(2,2,expectedTest2Buf);

        BOOST_CHECK_EQUAL(test1, expectedTest1);
        BOOST_CHECK_EQUAL(test2, expectedTest2);
    }

}


