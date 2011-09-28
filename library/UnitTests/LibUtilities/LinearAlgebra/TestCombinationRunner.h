///////////////////////////////////////////////////////////////////////////////
//
// File: TestCombinationRunner.h
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

#ifndef NEKTAR_UNIT_TESTS_LIB_UTILITIES_LINEAR_ALGEBRA_TEST_COMBINATION_RUNNER_H
#define NEKTAR_UNIT_TESTS_LIB_UTILITIES_LINEAR_ALGEBRA_TEST_COMBINATION_RUNNER_H

#include <boost/shared_ptr.hpp>
#include <LibUtilities/BasicUtils/BoostUtil.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>


namespace Nektar
{
    class DoAddition
    {
        public:
            template<typename LhsType, typename RhsType>
            typename expt::AddOp::ResultType<LhsType, RhsType>::type
            operator()(const LhsType& lhs, const RhsType& rhs) const
            {
                return lhs + rhs;
            }
    };

    class DoSubtraction
    {
        public:
            template<typename LhsType, typename RhsType>
            typename expt::SubtractOp::ResultType<LhsType, RhsType>::type 
            operator()(const LhsType& lhs, const RhsType& rhs) const
            {
                return lhs - rhs;
            }
    };

    class DoMultiplication
    {
        public:
            template<typename LhsType, typename RhsType>
            typename expt::MultiplyOp::ResultType<LhsType, RhsType>::type 
            operator()(const LhsType& lhs, const RhsType& rhs) const
            {
                return lhs * rhs;
            }
    };

    class DoDivision
    {
        public:
            template<typename LhsType, typename RhsType>
            typename expt::DivideOp::ResultType<LhsType, RhsType>::type 
            operator()(const LhsType& lhs, const RhsType& rhs) const
            {
                return lhs / rhs;
            }
    };

    template<typename LhsScaledInnerMatrixType,
             typename LhsBlockInnerMatrixType,
             typename RhsScaledInnerMatrixType,
             typename RhsBlockInnerMatrixType,
             typename OpType>
    void RunAllTestCombinations(const NekMatrix<NekDouble, StandardMatrixTag>& l1,
                                const NekMatrix<NekMatrix<NekDouble, LhsScaledInnerMatrixType>, ScaledMatrixTag>& l2,
                                const NekMatrix<NekMatrix<NekDouble, LhsBlockInnerMatrixType>, BlockMatrixTag>& l3,
                                const NekMatrix<NekDouble, StandardMatrixTag>& r1,
                                const NekMatrix<NekMatrix<NekDouble, RhsScaledInnerMatrixType>, ScaledMatrixTag>& r2,
                                const NekMatrix<NekMatrix<NekDouble, RhsBlockInnerMatrixType>, BlockMatrixTag>& r3,
                                const NekMatrix<NekDouble, StandardMatrixTag>& result,
                                const OpType& f)
    {
        BOOST_CHECK_EQUAL(f(l1, r1), result);
        BOOST_CHECK_EQUAL(f(l1, r2), result);
        BOOST_CHECK_EQUAL(f(l1, r3), result);
        BOOST_CHECK_EQUAL(f(l2, r1), result);
        BOOST_CHECK_EQUAL(f(l2, r2), result);
        BOOST_CHECK_EQUAL(f(l2, r3), result);
        BOOST_CHECK_EQUAL(f(l3, r1), result);
        BOOST_CHECK_EQUAL(f(l3, r2), result);
        BOOST_CHECK_EQUAL(f(l3, r3), result);
    }

    //template<typename NumberType>
    //void GenerateFullMatrices(const NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag>& m1,
    //    NumberType scale, unsigned int blockRows, unsigned int blockColumns,
    //    boost::shared_ptr<NekMatrix<NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag> >& m2,
    //    boost::shared_ptr<NekMatrix<NekMatrix<NumberType>, FullMatrixTag, BlockMatrixTag> >& m3)
    //{
    //    NumberType* inner_values = new NumberType[m1.GetStorageSize()];
    //    std::transform(m1.begin(), m1.end(), inner_values, boost::bind(std::divides<NumberType>(), _1, scale));

    //    boost::shared_ptr<NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> > inner(
    //        new NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag>(m1.GetRows(), m1.GetColumns(), inner_values, m1.GetPolicySpecificDataHolderType())); 
    //    m2 = MakePtr(new NekMatrix<NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>(scale, inner));

    //    unsigned int numberOfRows = m1.GetRows()/blockRows;
    //    unsigned int numberOfColumns = m1.GetColumns()/blockColumns;
    //    m3 = MakePtr(new NekMatrix<NekMatrix<NumberType>, FullMatrixTag, BlockMatrixTag>(blockRows, blockColumns, numberOfRows, numberOfColumns));

    //    for(unsigned int blockRow = 0; blockRow < blockRows; ++blockRow)
    //    {
    //        for(unsigned int blockColumn = 0; blockColumn < blockColumns; ++blockColumn)
    //        {
    //            boost::shared_ptr<NekMatrix<NumberType> > block(new NekMatrix<NumberType>(numberOfRows, numberOfColumns));
    //            for(unsigned int i = 0; i < numberOfRows; ++i)
    //            {
    //                for(unsigned int j = 0; j < numberOfRows; ++j)
    //                {
    //                    (*block)(i,j) = m1(numberOfRows*blockRow + i, numberOfColumns*blockColumn + j);
    //                }
    //            }
    //            
    //            m3->SetBlock(blockRow, blockColumn, block);
    //        }
    //    }

    //    delete [] inner_values;


    //    m1 = MakePtr(new NekMatrix<NekDouble, FullMatrixTag, StandardMatrixTag>(4, 4, values));
    //    
    //    double inner_values[16];
    //    std::transform(values, values+16, inner_values, boost::bind(std::divides<NekDouble>(), _1, scale));

    //    boost::shared_ptr<NekMatrix<NekDouble> > inner(
    //        new NekMatrix<NekDouble>(4, 4, inner_values)); 
    //    m2 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, ScaledMatrixTag>(scale, inner));
    //    
    //    double block_1_values[] = {values[0], values[1], 
    //                        values[4], values[5]};
    //    double block_2_values[] = {values[2], values[3],
    //                        values[6], values[7]};
    //    double block_3_values[] = {values[8], values[9], 
    //                        values[12], values[13]};
    //    double block_4_values[] = {values[10], values[11],
    //                        values[14], values[15]};
    //    boost::shared_ptr<NekMatrix<NekDouble> > block1(new NekMatrix<NekDouble>(2, 2, block_1_values));
    //    boost::shared_ptr<NekMatrix<NekDouble> > block2(new NekMatrix<NekDouble>(2, 2, block_2_values));
    //    boost::shared_ptr<NekMatrix<NekDouble> > block3(new NekMatrix<NekDouble>(2, 2, block_3_values));
    //    boost::shared_ptr<NekMatrix<NekDouble> > block4(new NekMatrix<NekDouble>(2, 2, block_4_values));
    //    
    //    m3 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, BlockMatrixTag>(2, 2, 2, 2));
    //    m3->SetBlock(0,0, block1);
    //    m3->SetBlock(1,0, block2);
    //    m3->SetBlock(0,1, block3);
    //    m3->SetBlock(1,1, block4);
    //}

    //void GenerateUpperTriangularMatrices(NekDouble values[], NekDouble scale,
    //    boost::shared_ptr<NekMatrix<NekDouble, UpperTriangularMatrixTag, StandardMatrixTag> >& m1,
    //    boost::shared_ptr<NekMatrix<NekMatrix<NekDouble, UpperTriangularMatrixTag>, UpperTriangularMatrixTag, ScaledMatrixTag> >& m2,
    //    boost::shared_ptr<NekMatrix<NekMatrix<NekDouble>, UpperTriangularMatrixTag, BlockMatrixTag> >& m3)
    //{
    //    m1 = MakePtr(new NekMatrix<NekDouble, UpperTriangularMatrixTag, StandardMatrixTag>(4, 4, values));
    //    
    //    double inner_values[10];
    //    std::transform(values, values+10, inner_values, boost::bind(std::divides<NekDouble>(), _1, scale));

    //    boost::shared_ptr<NekMatrix<NekDouble, UpperTriangularMatrixTag> > inner(
    //        new NekMatrix<NekDouble, UpperTriangularMatrixTag>(4, 4, inner_values)); 
    //    m2 = MakePtr(new NekMatrix<NekMatrix<NekDouble, UpperTriangularMatrixTag>, UpperTriangularMatrixTag, ScaledMatrixTag>(scale, inner));
    //    
    //    double block_1_values[] = {values[0], 0.0,
    //                               values[1], values[2]};
    //    double block_2_values[] = {values[3], values[4],
    //                               values[6], values[7]};
    //    double block_4_values[] = {values[5], 0.0,
    //                               values[8], values[9]};
    //    boost::shared_ptr<NekMatrix<NekDouble> > block1(new NekMatrix<NekDouble>(2, 2, block_1_values));
    //    boost::shared_ptr<NekMatrix<NekDouble> > block2(new NekMatrix<NekDouble>(2, 2, block_2_values));
    //    boost::shared_ptr<NekMatrix<NekDouble> > block4(new NekMatrix<NekDouble>(2, 2, block_4_values));
    //    
    //    m3 = MakePtr(new NekMatrix<NekMatrix<NekDouble>, UpperTriangularMatrixTag, BlockMatrixTag>(2, 2, 2, 2));
    //    m3->SetBlock(0,0, block1);
    //    m3->SetBlock(0,1, block2);
    //    m3->SetBlock(1,1, block4);
    //}

    template<typename NumberType>
    void GenerateMatrices(const NekMatrix<NumberType, StandardMatrixTag>& m1,
        NumberType scale, unsigned int blockRows, unsigned int blockColumns,
        boost::shared_ptr<NekMatrix<NekMatrix<NumberType, StandardMatrixTag>, ScaledMatrixTag> >& m2,
        boost::shared_ptr<NekMatrix<NekMatrix<NumberType>, BlockMatrixTag> >& m3)
    {
        NumberType* inner_values = new NumberType[m1.GetStorageSize()];
        std::transform(m1.begin(), m1.end(), inner_values, boost::bind(std::divides<NumberType>(), _1, scale));
        MatrixStorage s = m1.GetType();
        
        boost::shared_ptr<NekMatrix<NumberType, StandardMatrixTag> > inner(
            new NekMatrix<NumberType, StandardMatrixTag>(m1.GetRows(), m1.GetColumns(), inner_values, s, m1.GetNumberOfSubDiagonals(), 
            m1.GetNumberOfSuperDiagonals())); 
        m2 = MakePtr(new NekMatrix<NekMatrix<NumberType, StandardMatrixTag>, ScaledMatrixTag>(scale, inner));

        unsigned int numberOfRows = m1.GetRows()/blockRows;
        unsigned int numberOfColumns = m1.GetColumns()/blockColumns;
        m3 = MakePtr(new NekMatrix<NekMatrix<NumberType>, BlockMatrixTag>(blockRows, blockColumns, numberOfRows, numberOfColumns));

        for(unsigned int blockRow = 0; blockRow < blockRows; ++blockRow)
        {
            for(unsigned int blockColumn = 0; blockColumn < blockColumns; ++blockColumn)
            {
                boost::shared_ptr<NekMatrix<NumberType> > block(new NekMatrix<NumberType>(numberOfRows, numberOfColumns));
                for(unsigned int i = 0; i < numberOfRows; ++i)
                {
                    for(unsigned int j = 0; j < numberOfRows; ++j)
                    {
                        (*block)(i,j) = m1(numberOfRows*blockRow + i, numberOfColumns*blockColumn + j);
                    }
                }
                
                m3->SetBlock(blockRow, blockColumn, block);
            }
        }

        delete [] inner_values;
    }
}

#endif //NEKTAR_UNIT_TESTS_LIB_UTILITIES_LINEAR_ALGEBRA_TEST_COMBINATION_RUNNER_H
