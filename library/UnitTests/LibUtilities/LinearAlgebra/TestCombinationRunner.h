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

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

#include <memory>

namespace Nektar
{
    template<typename LhsScaledInnerMatrixType,
             typename LhsBlockInnerMatrixType,
             typename RhsScaledInnerMatrixType,
             typename RhsBlockInnerMatrixType>
    void RunAllAddCombinations(const NekMatrix<NekDouble, StandardMatrixTag>& l1,
                               const NekMatrix<NekMatrix<NekDouble, LhsScaledInnerMatrixType>, ScaledMatrixTag>& l2,
                               const NekMatrix<NekMatrix<NekDouble, LhsBlockInnerMatrixType>, BlockMatrixTag>& l3,
                               const NekMatrix<NekDouble, StandardMatrixTag>& r1,
                               const NekMatrix<NekMatrix<NekDouble, RhsScaledInnerMatrixType>, ScaledMatrixTag>& r2,
                               const NekMatrix<NekMatrix<NekDouble, RhsBlockInnerMatrixType>, BlockMatrixTag>& r3,
                               const NekMatrix<NekDouble, StandardMatrixTag>& result)
    {
        BOOST_CHECK_EQUAL(l1 + r1, result);
        BOOST_CHECK_EQUAL(l1 + r2, result);
        BOOST_CHECK_EQUAL(l1 + r3, result);
        BOOST_CHECK_EQUAL(l2 + r1, result);
        BOOST_CHECK_EQUAL(l2 + r2, result);
        BOOST_CHECK_EQUAL(l2 + r3, result);
        BOOST_CHECK_EQUAL(l3 + r1, result);
        BOOST_CHECK_EQUAL(l3 + r2, result);
        BOOST_CHECK_EQUAL(l3 + r3, result);
    }

    template<typename LhsScaledInnerMatrixType,
             typename LhsBlockInnerMatrixType,
             typename RhsScaledInnerMatrixType,
             typename RhsBlockInnerMatrixType>
    void RunAllSubCombinations(const NekMatrix<NekDouble, StandardMatrixTag>& l1,
                               const NekMatrix<NekMatrix<NekDouble, LhsScaledInnerMatrixType>, ScaledMatrixTag>& l2,
                               const NekMatrix<NekMatrix<NekDouble, LhsBlockInnerMatrixType>, BlockMatrixTag>& l3,
                               const NekMatrix<NekDouble, StandardMatrixTag>& r1,
                               const NekMatrix<NekMatrix<NekDouble, RhsScaledInnerMatrixType>, ScaledMatrixTag>& r2,
                               const NekMatrix<NekMatrix<NekDouble, RhsBlockInnerMatrixType>, BlockMatrixTag>& r3,
                               const NekMatrix<NekDouble, StandardMatrixTag>& result)
    {
        BOOST_CHECK_EQUAL(l1 - r1, result);
        BOOST_CHECK_EQUAL(l1 - r2, result);
        BOOST_CHECK_EQUAL(l1 - r3, result);
        BOOST_CHECK_EQUAL(l2 - r1, result);
        BOOST_CHECK_EQUAL(l2 - r2, result);
        BOOST_CHECK_EQUAL(l2 - r3, result);
        BOOST_CHECK_EQUAL(l3 - r1, result);
        BOOST_CHECK_EQUAL(l3 - r2, result);
        BOOST_CHECK_EQUAL(l3 - r3, result);
    }

    template<typename NumberType>
    void GenerateMatrices(const NekMatrix<NumberType, StandardMatrixTag>& m1,
        NumberType scale, unsigned int blockRows, unsigned int blockColumns,
        std::shared_ptr<NekMatrix<NekMatrix<NumberType, StandardMatrixTag>, ScaledMatrixTag> >& m2,
        std::shared_ptr<NekMatrix<NekMatrix<NumberType>, BlockMatrixTag> >& m3)
    {
        NumberType* inner_values = new NumberType[m1.GetStorageSize()];
        std::transform(m1.begin(), m1.end(), inner_values, std::bind(std::divides<NumberType>(), std::placeholders::_1, scale));
        MatrixStorage s = m1.GetType();
        
        std::shared_ptr<NekMatrix<NumberType, StandardMatrixTag> > inner(
            new NekMatrix<NumberType, StandardMatrixTag>(m1.GetRows(), m1.GetColumns(), inner_values, s, m1.GetNumberOfSubDiagonals(), 
            m1.GetNumberOfSuperDiagonals())); 
        m2 = std::make_shared<NekMatrix<NekMatrix<NumberType, StandardMatrixTag>, ScaledMatrixTag>>(scale, inner);

        unsigned int numberOfRows = m1.GetRows()/blockRows;
        unsigned int numberOfColumns = m1.GetColumns()/blockColumns;
        m3 = std::make_shared<NekMatrix<NekMatrix<NumberType>, BlockMatrixTag>>(blockRows, blockColumns, numberOfRows, numberOfColumns);

        for(unsigned int blockRow = 0; blockRow < blockRows; ++blockRow)
        {
            for(unsigned int blockColumn = 0; blockColumn < blockColumns; ++blockColumn)
            {
                std::shared_ptr<NekMatrix<NumberType> > block(new NekMatrix<NumberType>(numberOfRows, numberOfColumns));
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
