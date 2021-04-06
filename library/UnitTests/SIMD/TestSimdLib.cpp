///////////////////////////////////////////////////////////////////////////////
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
/// The above copyright notice and this permission notice shall be included
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

#include <LibUtilities/SimdLib/tinysimd.hpp>
#include <LibUtilities/SimdLib/io.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/core/ignore_unused.hpp>


#include <array>
#include <cmath>
#include <iostream>

// type in use
#if defined(__SSE2__) && defined(NEKTAR_ENABLE_SIMD_SSE2)
    #define USING_SSE2
#endif
#if defined(__AVX2__) && defined(NEKTAR_ENABLE_SIMD_AVX2)
    #define USING_AVX2
#endif
#if defined(__AVX512__) && defined(NEKTAR_ENABLE_SIMD_AVX512)
    #define USING_AVX512
#endif


namespace Nektar
{
namespace SimdLibTests
{
    using namespace tinysimd;

    BOOST_AUTO_TEST_CASE(SimdLib_width_alignment)
    {
        std::size_t width, alignment;

        #if !defined(USING_SSE2) && !defined(USING_AVX2) && !defined(USING_AVX512)
        std::cout << "scalar" << std::endl;
        // std::int32_t aka (usually) int
        width = simd<std::int32_t>::width;
        alignment = simd<std::int32_t>::alignment;
        BOOST_CHECK_EQUAL(width, 1);
        BOOST_CHECK_EQUAL(alignment, 4);
        // std::int64_t aka (usually) long
        width = simd<std::int64_t>::width;
        alignment = simd<std::int64_t>::alignment;
        BOOST_CHECK_EQUAL(width, 1);
        BOOST_CHECK_EQUAL(alignment, 8);
        // double
        width = simd<double>::width;
        alignment = simd<double>::alignment;
        BOOST_CHECK_EQUAL(width, 1);
        BOOST_CHECK_EQUAL(alignment, 8);
        #endif

        #if defined(USING_SSE2) && !defined(USING_AVX2) && !defined(USING_AVX512)
        std::cout << "sse2" << std::endl;
        // std::int32_t
        width = simd<std::int32_t>::width;
        alignment = simd<std::int32_t>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 16);
        #endif

        #if defined(USING_AVX2) && !defined(USING_AVX512)
        std::cout << "avx2" << std::endl;
        // int (sse2int4)
        width = simd<int>::width;
        alignment = simd<int>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 16);
        // long
        width = simd<long>::width;
        alignment = simd<long>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 32);
        // double
        width = simd<double>::width;
        alignment = simd<double>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 32);
        #endif

        #if defined(USING_AVX512)
        std::cout << "avx512" << std::endl;
        // int (avx2int8)
        width = simd<int>::width;
        alignment = simd<int>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 32);
        // long
        width = simd<long>::width;
        alignment = simd<long>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 64);
        // double
        width = simd<double>::width;
        alignment = simd<double>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 64);
        #endif

    }

    using vec_t = simd<double>;

    BOOST_AUTO_TEST_CASE(SimdLib_type_traits)
    {
        BOOST_CHECK_EQUAL(has_width<double>::value, false);
        BOOST_CHECK_EQUAL(has_width<vec_t>::value, true);

        BOOST_CHECK_EQUAL(has_alignment<double>::value, false);
        BOOST_CHECK_EQUAL(has_alignment<vec_t>::value, true);

        BOOST_CHECK_EQUAL(has_scalarType<double>::value, false);
        BOOST_CHECK_EQUAL(has_scalarType<vec_t>::value, true);

        BOOST_CHECK_EQUAL(is_vector<double>::value, false);
        BOOST_CHECK_EQUAL(is_vector<vec_t>::value, true);

        BOOST_CHECK_EQUAL(is_vector_floating_point<double>::value, false);
        BOOST_CHECK_EQUAL(is_vector_floating_point<simd<int>>::value, false);
        BOOST_CHECK_EQUAL(is_vector_floating_point<vec_t>::value, true);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_mem_size)
    {
        BOOST_CHECK_EQUAL(sizeof(vec_t), sizeof(double)*vec_t::width);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_ctors)
    {
        vec_t avec1;

        vec_t::scalarType ascalar = 0;
        vec_t avec2(ascalar);
        vec_t avec3{ascalar};
        vec_t avec4 = ascalar;

        vec_t avec5(avec2);
        vec_t avec6{avec4};

        vec_t avec7(avec2._data);
        vec_t avec8{avec2._data};

        vec_t::vectorType anative;
        vec_t avec9(anative);
        vec_t avec10{anative};


        boost::ignore_unused(avec1, avec3, avec5, avec6, avec7, avec8, avec9,
            avec10);
    }


    BOOST_AUTO_TEST_CASE(SimdLib_load)
    {
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        vec_t avec;
        avec.load(ascalararr.data());
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_implicit)
    {
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        vec_t avec;
        avec = *(reinterpret_cast<vec_t*>(ascalararr.data()));
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_aligned)
    {
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        vec_t avec;
        avec.load(ascalararr.data(), is_aligned);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_unaligned)
    {
        std::array<double, vec_t::width> ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        vec_t avec;
        avec.load(ascalararr.data(), is_not_aligned);
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(vec_t::alignment) std::array<double, vec_t::width> ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        avec.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store_aligned)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(vec_t::alignment) std::array<double, vec_t::width> ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        avec.store(ascalararr.data(), is_aligned);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store_unaligned)
    {
        double val = 4.0;
        vec_t avec(val);
        std::array<double, vec_t::width> ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        avec.store(ascalararr.data(), is_not_aligned);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_store_non_temporal)
    {
        double val = 4.0;
        vec_t avec(val);
        alignas(vec_t::alignment) std::array<double, vec_t::width> ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        avec.store(ascalararr.data(), is_not_reused);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_subscript_assign_read)
    {
        vec_t avec;
        std::array<double, vec_t::width> ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            ascalararr[i] = i;
        }


        for (size_t i = 0; i < vec_t::width; ++i)
        {
            avec[i] = ascalararr[i];
        }

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], avec[i]);
        }
    }

    #if defined(__AVX2__) && defined(NEKTAR_ENABLE_AVX2)
    BOOST_AUTO_TEST_CASE(SimdLib_load_any)
    {
        vec_t avec;
        double* p0, * p1, * p2, * p3;
        std::array<double, 16> ascalararr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

        p0 =  ascalararr.data();
        p1 =  ascalararr.data() + 3;
        p2 =  ascalararr.data() + 5;
        p3 =  ascalararr.data() + 6;

        avec.load(p0, p1, p2, p3);

        BOOST_CHECK_EQUAL(ascalararr[0], avec[0]);
        BOOST_CHECK_EQUAL(ascalararr[3], avec[1]);
        BOOST_CHECK_EQUAL(ascalararr[5], avec[2]);
        BOOST_CHECK_EQUAL(ascalararr[6], avec[3]);

    }
    #endif

    BOOST_AUTO_TEST_CASE(SimdLib_gather64)
    {
        vec_t avec;
        using index_t = simd<size_t>;
        index_t aindexvec;

        // create and fill index
        std::array<size_t, vec_t::width> aindex;
        aindex[0] = 0;
        if (vec_t::width > 2)
        {
            aindex[1] = 3;
            aindex[2] = 5;
            aindex[3] = 6;
        }
        if (vec_t::width > 4)
        {
            aindex[4] = 8;
            aindex[5] = 15;
            aindex[6] = 20;
        }

        // load index
        aindexvec.load(aindex.data(), is_not_aligned);

        // create and fill scalar array
        constexpr size_t scalarArraySize = 32;
        std::array<double, scalarArraySize> ascalararr;
        for (size_t i = 0; i < scalarArraySize; ++i)
        {
            ascalararr[i] = i;
        }

        avec.gather(ascalararr.data(), aindexvec);

        // check
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[aindex[i]], avec[i]);
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_scatter64)
    {
        vec_t avec;
        using index_t = simd<size_t>;
        index_t aindexvec;

        // create and fill index
        std::array<size_t, vec_t::width> aindex;
        aindex[0] = 1;
        if (vec_t::width > 2)
        {
            aindex[1] = 3;
            aindex[2] = 5;
            aindex[3] = 6;
        }
        if (vec_t::width > 4)
        {
            aindex[4] = 8;
            aindex[5] = 15;
            aindex[6] = 20;
            aindex[7] = 30;
        }

         // load index
        aindexvec.load(aindex.data(), is_not_aligned);

        // create scalar array
        constexpr size_t scalarArraySize = 32;
        std::array<double, scalarArraySize> ascalararr;

        // fill vector
        avec[0] = 10;
        if (vec_t::width > 2)
        {
            avec[1] =  9;
            avec[2] =  8;
            avec[3] =  7;
        }
        if (vec_t::width > 4)
        {
            avec[4] =  4;
            avec[5] =  3;
            avec[6] =  2;
            avec[7] =  1;
        }

        avec.scatter(ascalararr.data(), aindexvec);

        // check
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(avec[i], ascalararr[aindex[i]]);
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_add_unary)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        vec_t res(val1);
        vec_t avec(val2);
        res += avec;
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        res.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_sub_unary)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        vec_t res(val1);
        vec_t avec(val2);
        res -= avec;
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        res.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 - val2);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_mul_unary)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        vec_t res(val1);
        vec_t avec(val2);
        res *= avec;
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        res.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 * val2);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_div_unary)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        vec_t res(val1);
        vec_t avec(val2);
        res /= avec;
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        res.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 / val2);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_add_mul)
    {
        double val1 = -4.0;
        double val2 =  2.0;
        double val3 =  2.0;
        vec_t avec1(val1);
        vec_t avec2(val2);
        vec_t avec3(val3);
        vec_t res = avec1 + avec2 * avec3;
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        res.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2 + val3);
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_sqrt)
    {
        double val = 4.0;
        vec_t avec(val);
        vec_t asqrt = sqrt(avec);
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        asqrt.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::sqrt(val));
        }
    }


    BOOST_AUTO_TEST_CASE(SimdLib_abs)
    {
        double val = -4.0;
        vec_t avec(val);
        vec_t aabs = abs(avec);
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        aabs.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::abs(val));
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_log)
    {
        double val = 4.0;
        vec_t avec(val);
        vec_t alog = log(avec);
        alignas(vec_t::alignment) std::array<double, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        alog.store(ascalararr.data());

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(ascalararr[i], std::log(val));
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_greater)
    {
        double aval = 4.0;
        vec_t avec(aval);
        using mask_t = simd<bool>;
        mask_t amask;

        amask = avec > avec;
        // check
        alignas(vec_t::alignment) std::array<std::uint64_t, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        amask.store(ascalararr.data());
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            // type conversion make lvalue in rvalue, needed pre-c++17
            BOOST_CHECK_EQUAL(ascalararr[i], (std::uint64_t)mask_t::false_v);
        }

        double bval = 3.0;
        vec_t bvec(bval);

        amask = avec > bvec;
        // check
        amask.store(ascalararr.data());
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            // type conversion make lvalue in rvalue, needed pre-c++17
            BOOST_CHECK_EQUAL(ascalararr[i], (std::uint64_t)mask_t::true_v);
        }

        double cval = 5.0;
        vec_t cvec(cval);

        amask = avec > cvec;
        // check
        amask.store(ascalararr.data());
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            // type conversion make lvalue in rvalue, needed pre-c++17
            BOOST_CHECK_EQUAL(ascalararr[i], (std::uint64_t)mask_t::false_v);
        }


        if (vec_t::width == 4)
        {
            alignas(vec_t::alignment) std::array<double, 4>
                ascalararr2{{1.0,2.0,3.0,4.0}}; // double brace to deal with gcc 4.8.5 ...
            double dval = 2.0;
            vec_t dvec(dval);
            vec_t evec;
            evec.load(ascalararr2.data());

            simd<bool> amask;
            amask = dvec > evec;
            // check
            for (size_t i = 0; i < vec_t::width; ++i)
            {
                BOOST_CHECK_EQUAL(static_cast<bool>(amask[i]), dvec[i] > evec[i]);
            }

        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_logic_and)
    {
        double aval = 4.0;
        vec_t avec(aval);
        using mask_t = simd<bool>;
        mask_t amask;

        alignas(vec_t::alignment) std::array<std::uint64_t, vec_t::width>
            ascalararr{{}}; // double brace to deal with gcc 4.8.5 ...
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            ascalararr[i] = mask_t::true_v;
        }
        amask.load(ascalararr.data());

        // check
        BOOST_CHECK_EQUAL(amask && false, false);
        BOOST_CHECK_EQUAL(amask && true, true);

        for (size_t i = 0; i < vec_t::width; ++i)
        {
            ascalararr[i] = mask_t::false_v;
        }
        amask.load(ascalararr.data());

        // check
        BOOST_CHECK_EQUAL(amask && false, false);
        BOOST_CHECK_EQUAL(amask && true, false);

        if (vec_t::width > 1)
        {
            ascalararr[0] = mask_t::true_v;
            // check
            BOOST_CHECK_EQUAL(amask && false, false);
            BOOST_CHECK_EQUAL(amask && true, false);
        }
    }

    BOOST_AUTO_TEST_CASE(SimdLib_load_interleave_unload)
    {
        constexpr size_t nDof{5};
        // no padding in load_interleave deinterleave_store
        constexpr size_t nEle{vec_t::width * 5};
        constexpr size_t nDofBlock = nDof * vec_t::width;

        constexpr size_t size{nDof*nEle};
        std::array<double,size> dofScalarArr{{}};
        for (size_t i = 0; i < size; ++i)
        {
            dofScalarArr[i] = i;
        }

        // number of blocks
        size_t nBlock = nEle / vec_t::width;

        // aligned vector
        std::vector<vec_t, allocator<vec_t>> dofVectorArr(nDof);

        double* dataPtr = dofScalarArr.data();
        // loop over blocks vec_t::width elements at the time
        for (size_t b = 0; b < nBlock; ++b)
        {
            // load
            load_interleave(dataPtr, nDof, dofVectorArr);

            // manipulate each block
            for (size_t j = 0; j < nDof; ++j)
            {
                dofVectorArr[j] = dofVectorArr[j] + j;
            }

            // store
            deinterleave_store(dofVectorArr, nDof, dataPtr);
            dataPtr += nDofBlock;
        }

        // check
        for (size_t b = 0, i = 0; b < nBlock; ++b)
        {
            for (size_t j = 0; j < nDof; ++j, ++i)
            {
                BOOST_CHECK_EQUAL(dofScalarArr[i], i + j);
            }
        }

    }

    BOOST_AUTO_TEST_CASE(SimdLib_io)
    {
        vec_t avec(3.14);
        std::cout << avec << std::endl;
    }

}
}
