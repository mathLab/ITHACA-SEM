///////////////////////////////////////////////////////////////////////////////
//
// File: TestTwoParameters.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>

namespace Nektar
{
    namespace TestTwoParameters
    {
        BOOST_AUTO_TEST_CASE(TestComplicated1)
        {
            //typedef Node<Node<double>, expt::MultiplyOp, Node<NekMatrix<double> > > Node1;
            //typedef Node<Node<NekMatrix<double> >, expt::AddOp, Node<NekMatrix<double> > > Node2;
            //typedef Node<Node<double>, expt::MultiplyOp, Node2> Node3;
            //typedef Node<Node1, expt::AddOp, Node3> Expression;
            //
            //typedef Expression::Indices Indices;
            //typedef Node3::Indices Node3Indices;
            //
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Indices, 0>::type::value == 0 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Indices, 1>::type::value == 1 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Indices, 2>::type::value == 2 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Indices, 3>::type::value == 3 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Indices, 4>::type::value == 4 ));
            //
            //BOOST_STATIC_ASSERT(( boost::mpl::size<Node3Indices>::type::value == 3 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Node3Indices, 0>::type::value == 0 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Node3Indices, 1>::type::value == 1 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Node3Indices, 2>::type::value == 2 ));

            //typedef RemoveUnecessaryTemporaries<Node3, Node3Indices>::TransformedIndicesType Node3TransformedIndices;
            //BOOST_STATIC_ASSERT(( boost::mpl::size<Node3TransformedIndices>::type::value == 3 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Node3TransformedIndices, 0>::type::value == 1 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Node3TransformedIndices, 1>::type::value == 2 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<Node3TransformedIndices, 2>::type::value == 0 ));

            //typedef RemoveUnecessaryTemporaries<Expression, Indices>::TransformedIndicesType TransformedIndicesType;
            //typedef RemoveUnecessaryTemporaries<Expression, Indices>::TransformedNodeType TransformedNodes;

            //BOOST_STATIC_ASSERT(( boost::mpl::size<TransformedIndicesType>::type::value == 5 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 0>::type::value == 0 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 1>::type::value == 1 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 2>::type::value == 3 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 3>::type::value == 4 ));
            //BOOST_STATIC_ASSERT(( boost::mpl::at_c<TransformedIndicesType, 4>::type::value == 2 ));


            ////BOOST_MPL_ASSERT(( boost::is_same<TransformedNodes, Expression> ));
            //
            //NekDouble gmat = 2.0;
            //double lapp00_buf[] = {1, 2, 3, 4};
            //NekMatrix<double> lap00(2, 2, lapp00_buf);
            //
            //double lapp01_buf[] = {5, 6, 7, 8};
            //NekMatrix<double> lap01(2, 2, lapp01_buf);
            //
            //double lapp11_buf[] = {9, 10, 11, 12};
            //NekMatrix<double> lap11(2, 2, lapp11_buf);

            //NekMatrix<double> lap = gmat*lap00 +
            //    gmat*(lap01 + Transpose(lap01))  +
            //    gmat*lap11;

        }

        BOOST_AUTO_TEST_CASE(TestComplicated2)
        {
            typedef NekMatrix<double> Matrix;
            Matrix gmat[9][9];
            for(unsigned int i = 0; i < 9; ++i)
            {
                for(unsigned int j = 0; j < 9; ++j)
                {
                    gmat[i][j] = Matrix(2,2);
                }
            }

            Matrix lap00(2,2);
            Matrix lap01(2,2);
            Matrix lap02(2,2);
            Matrix lap11(2,2);
            Matrix lap12(2,2);
            Matrix lap22(2,2);
            
            Matrix lap = (gmat[0][0]*gmat[0][0] + gmat[3][0]*gmat[3][0]
                                        + gmat[6][0]*gmat[6][0])*lap00
                               + (gmat[1][0]*gmat[1][0] + gmat[4][0]*gmat[4][0]
                                        + gmat[7][0]*gmat[7][0])*lap11
                               + (gmat[2][0]*gmat[2][0] + gmat[5][0]*gmat[5][0]
                                        + gmat[8][0]*gmat[8][0])*lap22
                               + (gmat[0][0]*gmat[1][0] + gmat[3][0]*gmat[4][0]
                                        + gmat[6][0]*gmat[7][0])
                                 *(lap01 + Transpose(lap01))
                               + (gmat[0][0]*gmat[2][0] + gmat[3][0]*gmat[5][0]
                                        + gmat[6][0]*gmat[8][0])
                                 *(lap02 + Transpose(lap02))
                               + (gmat[1][0]*gmat[2][0] + gmat[4][0]*gmat[5][0]
                                        + gmat[7][0]*gmat[8][0])
                                 *(lap12 + Transpose(lap12));
            
        }

    }
}
