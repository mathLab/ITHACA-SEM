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

namespace Nektar
{
    // This test passes if it compiles.
    BOOST_AUTO_TEST_CASE(TestAssociativeTree1)
    {
        typedef Node< Node<NekMatrix<double> >, AddOp, Node<NekMatrix<double> > > Exp1;
        BOOST_MPL_ASSERT(( TreeIsAssociative<Exp1> ));

        typedef Node< Node<NekMatrix<double> >, AddOp, Exp1> Exp2;
        BOOST_MPL_ASSERT(( TreeIsAssociative<Exp2> ));

        typedef Node< Node<NekMatrix<double> >, AddOp, Exp2> Exp3;
        BOOST_MPL_ASSERT(( TreeIsAssociative<Exp3> ));
    }

    BOOST_AUTO_TEST_CASE(TestAssociativeTree2)
    {
        typedef Node< Node<NekMatrix<double> >, MultiplyOp, Node<NekMatrix<double> > > Exp1;
        BOOST_MPL_ASSERT(( TreeIsAssociative<Exp1> ));

        typedef Node< Node<NekMatrix<double> >, AddOp, Exp1> Exp2;
        BOOST_MPL_ASSERT(( boost::mpl::not_<TreeIsAssociative<Exp2> > ));

        typedef Node< Node<NekMatrix<double> >, AddOp, Exp2> Exp3;
        BOOST_MPL_ASSERT(( boost::mpl::not_<TreeIsAssociative<Exp3> > ));
    }
    
}













