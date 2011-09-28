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
#include <boost/tuple/tuple.hpp>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{


//    BOOST_AUTO_TEST_CASE(TestVectorCreation)
//    {
//        boost::tuple<int, double> a = boost::make_tuple(1, 2.3);
//        int b = boost::get<0>(a);
//        double c = boost::get<1>(a);

//        const unsigned int n = 10000000;
//        const unsigned int reps = 50;

//        //const unsigned int n = 100;
//        //const unsigned int reps = 1;
//        NekVector<double> result(n);
//        NekVector<double> v0(n);
//        NekVector<double> v1(n);
//        NekVector<double> v2(n);
//        NekVector<double> v3(n);

//        BOOST_AUTO(t, boost::make_tuple(boost::cref(v0), boost::cref(v1), boost::cref(v2), boost::cref(v3)));

//        typedef Node<NekVector<double> > ConstantNode;
//        typedef Node<ConstantNode, expt::AddOp, ConstantNode> T0;
//        typedef Node<T0, expt::AddOp, ConstantNode> T1;
//        typedef Node<T1, expt::AddOp, ConstantNode> T2;

//        //boost::fusion::vector<const NekVector<double>&, const NekVector<double>&, const NekVector<double>&, const NekVector<double>&> f(v0, v1, v2, v3);
//        T2::VectorType f(v0, v1, v2, v3);
//        typedef T2::Indices Indices;
//        Timer timer;
//        //timer.Start();
//        //for(unsigned int j = 0; j < reps; ++j)
//        //{
//        //    for(unsigned int i = 0; i < n; ++i)
//        //    {
//        //        result[i] = boost::get<0>(t)[i] + boost::get<1>(t)[i] + boost::get<2>(t)[i] + boost::get<3>(t)[i];
//        //    }
//        //}
//        //timer.Stop();
//        //std::cout << "Tuple Addition: " <<  timer.TimePerTest(1) << std::endl;

//        timer.Start();
//        for(unsigned int j = 0; j < reps; ++j)
//        {
//            for(unsigned int i = 0; i < n; ++i)
//            {
//                //result[i] = boost::fusion::at_c<0>(f)[i];
//                result[i] = boost::fusion::at_c<0>(f)[i] + boost::fusion::at_c<1>(f)[i] + boost::fusion::at_c<2>(f)[i] + boost::fusion::at_c<3>(f)[i];
//                //result[i] = 1.2;
//            }
//        }
//        timer.Stop();
//        std::cout << "Fusion Addition: " << timer.TimePerTest(1) << std::endl;

//        timer.Start();
//        for(unsigned int j = 0; j < reps; ++j)
//        {
//            for(unsigned int i = 0; i < n; ++i)
//            {
//                result[i] = v0[i] + v1[i] + v2[i] + v3[i];
//            }
//        }
//        timer.Stop();
//        std::cout << "Hand Coded Addition: " << timer.TimePerTest(1) << std::endl;

//        timer.Start();
//        for(unsigned int j = 0; j < reps; ++j)
//        {
//            result = v0 + v1 + v2 + v3 + v0 + v1 + v2 + v3 +v0;
//        }
//        timer.Stop();
//        std::cout << "Full Expression template code: " << timer.TimePerTest(1) << std::endl;

//        timer.Start();
//        for(unsigned int j = 0; j < reps; ++j)
//        {
//            Unroll<Indices, 0, 4>::Execute(result, f);
//        }
//        timer.Stop();

//        std::cout << "Manual Call to Unroll: " << timer.TimePerTest(1) << std::endl;
//    }
}
