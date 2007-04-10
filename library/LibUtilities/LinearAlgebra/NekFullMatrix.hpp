///////////////////////////////////////////////////////////////////////////////
//
// File: NekFullMatrix.hpp
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
// Description: Full matrix
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/call_traits.hpp>

#include <algorithm>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif


namespace Nektar
{
    template<typename DataType>
    class NekMatrixStoragePolicy<DataType, eFull> 
    {
        public: 
            static void Transpose(unsigned int& rows, unsigned int& columns,
                                  SharedArray<DataType>& data)
            {
                // TODO - Find a way to do this in-place.
                SharedArray<DataType> temp = MemoryManager::AllocateSharedArray<DataType>(data.GetSize());

                for(unsigned int row = 0; row < rows; ++row)
                {
                    for(unsigned int column = 0; column < columns; ++column)
                    {
                        unsigned int firstIndex = CalculateIndex(row, column, rows, columns);
                        unsigned int secondIndex = CalculateIndex(column, row, columns, rows);

                        temp[secondIndex] = data[firstIndex];
                    }
                }

                std::swap(rows, columns);
                std::swap(data, temp);
            }


            static typename boost::call_traits<DataType>::reference GetData(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns,
                                                                            SharedArray<DataType>& data)
            {
                return data[CalculateIndex(row, column, matrixRows, matrixColumns)];
            }

            static typename boost::call_traits<DataType>::const_reference GetConstData(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns,
                                                                            const SharedArray<DataType>& data)
            {
                return data[CalculateIndex(row, column, matrixRows, matrixColumns)];
            }

            static unsigned int CalculateIndex(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns)
            {
                return row*matrixColumns + column;
            }

            static unsigned int NumStorageElements(unsigned int rows, unsigned int cols)
            {
                return rows*cols;
            }

            static void Invert(unsigned int rows, unsigned int columns, SharedArray<DataType>& data)
            {
#ifdef NEKTAR_USING_LAPACK
                ASSERTL0(rows == columns, "Matrix Inversion only works for square arrays.");

                /// Incoming data is row major, make it column major for lapack calls.
                Transpose(rows, columns, data);

                int m = rows;
                int n = columns;
                int pivotSize = std::max(1, std::min(m, n));
                Nektar::NekIntSharedArray ipivot = Nektar::MemoryManager::AllocateSharedArray<int>(pivotSize);
                int info = 0;
                Lapack::Dgetrf(m, n, data.get(), m, ipivot.get(), info);

                if( info < 0 )
                {
                    std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
                            "th parameter had an illegal parameter for dgetrf";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                            boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
                    ASSERTL0(false, message.c_str());
                }

                unsigned int workSize = 64*n;
                Nektar::NekDoubleSharedArray work = Nektar::MemoryManager::AllocateSharedArray<NekDouble>(workSize);
                Lapack::Dgetri(n, data.get(), n, ipivot.get(), work.get(), workSize, info);

                if( info < 0 )
                {
                    std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
                            "th parameter had an illegal parameter for dgetri";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                            boost::lexical_cast<std::string>(info) + " is 0 from dgetri";
                    ASSERTL0(false, message.c_str());
                }

                // Put it back to row major form.
                Transpose(rows, columns, data);

#else
                // TODO
                BOOST_STATIC_ASSERT(0);
#endif //NEKTAR_USING_LAPACK

            }

            static void Factorize(unsigned int rows, unsigned int columns, SharedArray<DataType>& data)
            {

            }

        private:
            
    };

    template<typename DataType>
    class NekMatrixArithmeticPolicy<DataType, eFull>
    {
        public:
            template<MatrixBlockType BlockType, unsigned int space>
            static void PlusEqual(NekMatrix<DataType, eFull, BlockType, space>& lhs, 
                                  const NekMatrix<DataType, eFull, BlockType, space>& rhs)
            {
                ASSERTL0(lhs.GetRows() == rhs.GetRows() && lhs.GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");
                DataType* lhs_data = lhs.begin();
                const DataType* rhs_data = rhs.begin();

                for( ; lhs_data < lhs.end(); ++lhs_data, ++rhs_data )
                {
                    *lhs_data += *rhs_data;
                }
            }

            template<MatrixBlockType BlockType, unsigned int space>
            static void PlusEqual(NekMatrix<DataType, eFull, BlockType, space>& lhs, 
                                  const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
            {
                ASSERTL0(lhs.GetRows() == rhs.GetRows() && lhs.GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");

                for(unsigned int i = 0; i < lhs.GetRows(); ++i)
                {
                    lhs(i,i) += rhs(i,i);
                }
            }
    };

    template<typename DataType>
    class NekMatrixAssignmentPolicy<DataType, eFull>
    {
        public:
            template<MatrixBlockType BlockType, unsigned int space>
            static void Assign(NekMatrix<DataType, eFull, BlockType, space>& lhs, 
                               const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
            {
                lhs = NekMatrix<DataType, eFull, BlockType, space>(rhs.GetRows(), rhs.GetColumns(), DataType(0));
                for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                {
                    lhs(i,i) = rhs(i,i);
                }
            }
    };

//    template<typename DataType, unsigned int space>
//    class NekMatrix<DataType, eFull, eNormal, space>
//    {
//        public:
//            typedef NekMatrix<DataType, eFull, eNormal, space> ThisType;
//
//        public:
//            /// \brief Creates an empty, 0x0 matrix.  
//            NekMatrix() :
//                m_rows(0),
//                m_columns(0),
//                m_data()
//            {
//            }
//
//            /// \brief Create a matrix with the given size, initialized to the default value of DataType.
//            ///
//            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//            NekMatrix(unsigned int rows, unsigned int columns) :
//                m_rows(),
//                m_columns(),
//                m_data()
//            {
//                Initialize(rows, columns, DataType(0));
//            }
//
//            NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, PointerWrapper t = eCopy) :
//                m_rows(rows),
//                m_columns(columns),
//                m_data()
//            {
//                if( t == eCopy )
//                {
//                    m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
//                    std::copy(ptr, ptr + m_rows*m_columns, begin());
//                }
//                else
//                {
//                    m_data = SharedArray<DataType>(ptr, DeleteNothing<DataType>());
//                }
//            }
//             
//            NekMatrix(unsigned int rows, unsigned int columns, const DataType* ptr) :
//                m_rows(rows),
//                m_columns(columns),
//                m_data()
//            {
//                m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
//                std::copy(ptr, ptr + m_rows*m_columns, begin());
//            }
//                
//            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d) :
//                m_rows(),
//                m_columns(),
//                m_data()
//            {
//                Initialize(rows, columns, d);
//            }
//
//
//            NekMatrix(const NekMatrix<DataType, eFull, eNormal, space>& rhs) :
//                m_rows(rhs.m_rows),
//                m_columns(rhs.m_columns),
//                m_data(MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns))
//            {
//                std::copy(rhs.begin(), rhs.end(), begin());
//            }
//
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//
//            template<typename ExpressionPolicyType>
//            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
//                m_rows(rhs.GetMetadata().Rows),
//                m_columns(rhs.GetMetadata().Columns),
//                m_data(MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns))
//        {
//            BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eFull, eNormal, space> > ));
//            rhs.Apply(*this);
//        }
//#endif
//
//        NekMatrix<DataType, eFull, eNormal, space>& operator=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
//        {
//            NekMatrix<DataType, eFull, eNormal, space> temp(rhs);
//            Swap(temp);
//            return *this;
//        }
//
//
//        NekMatrix<DataType, eFull, eNormal, space>& operator=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
//        {
//            m_rows = rhs.GetRows();
//            m_columns = rhs.GetColumns();
//            m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
//            std::fill(begin(), end(), DataType(0));
//            
//            for(unsigned int i = 0; i < m_rows; ++i)
//            {
//                (*this)(i,i) = rhs(i,i);    
//            }
//            return *this;
//        }
//
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//        template<typename ExpressionPolicyType>
//        NekMatrix<DataType, eFull, eNormal, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
//        {
//            BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eFull, eNormal, space> > ));
//            m_rows = rhs.GetMetadata().Rows;
//            m_columns = rhs.GetMetadata().Columns;
//            m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
//            rhs.Apply(*this);
//            return *this;
//        }
//#endif
//
//        ~NekMatrix() 
//        {
//            m_rows = 0;
//            m_columns = 0;
//            m_data.reset();
//        }
//
//
//        unsigned int GetRows() const { return m_rows; }
//        unsigned int GetColumns() const { return m_columns; }
//
//        typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
//        {
//            ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
//            ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
//            return m_data[rowNumber*m_columns + colNumber];
//        }
//
//        typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
//        {
//            ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
//            ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
//            return m_data[rowNumber*m_columns + colNumber];
//        }
//
//        DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
//        {
//            return m_data.get() + rowNumber*m_columns + colNumber;
//        }
//
//        SharedArray<DataType> GetPtr()
//        {
//            return m_data;
//        }
//
//        // I can't figure out a way to return a constant pointer to the data 
//        // in a shared array.
//        //SharedArray<const DataType> GetPtr() const
//        //{
//        //    return SharedArray<const DataType>(m_data);
//        //}
//
//        typedef DataType* iterator;
//        typedef const DataType* const_iterator;
//
//        iterator begin() { return m_data.get(); }
//
//        iterator end() { return m_data.get() + m_columns*m_rows; }
//
//        const_iterator begin() const { return m_data.get(); }
//
//        const_iterator end() const { return m_data.get() + m_columns*m_rows; }
//
//        void Negate()
//        {
//            for(iterator iter = begin(); iter != end(); ++iter)
//            {
//                *iter = -*iter;
//            }
//        }
//        
//        void Transpose()
//        {
//            for(unsigned int row = 0; row < m_rows; ++row)
//            {
//                for(unsigned int column = row+1; column < m_columns; ++column)
//                {
//                    std::swap( (*this)(row, column), (*this)(column, row));
//                }
//            }
//        }
//
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//        expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eFull, eNormal, space> > >, expt::NegateOp> > operator-() const
//        {
//            return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eFull, eNormal, space> > >, expt::NegateOp> >(
//                    expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eFull, eNormal, space> > >(*this));
//        }
//#endif
//
//        NekMatrix<DataType, eFull, eNormal, space> operator+=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
//        {
//            ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");
//            DataType* lhs_data = begin();
//            const DataType* rhs_data = rhs.begin();
//
//            for( ; lhs_data < end(); ++lhs_data, ++rhs_data )
//            {
//                *lhs_data += *rhs_data;
//            }
//
//            return *this;
//        }
//
//        // This is wrong as well.  What if this is diagonal?
//        // enable if on the output.
//        NekMatrix<DataType, eFull, eNormal, space> operator+=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
//        {
//            ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");
//
//            for(unsigned int i = 0; i < rhs.GetRows(); ++i)
//            {
//                (*this)(i,i) += rhs(i, i);
//            }
//
//            return *this;
//        }
//
//        NekMatrix<DataType, eFull, eNormal, space> operator-=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
//        {
//            ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");
//            DataType* lhs_data = begin();
//            const DataType* rhs_data = rhs.begin();
//
//            for( ; lhs_data < end(); ++lhs_data, ++rhs_data )
//            {
//                *lhs_data -= *rhs_data;
//            }
//
//            return *this;
//        }
//
//        // This is wrong as well.  What if this is diagonal?
//        // enable if on the output.
//        NekMatrix<DataType, eFull, eNormal, space> operator-=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
//        {
//            ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");
//
//            for(unsigned int i = 0; i < rhs.GetRows(); ++i)
//            {
//                (*this)(i,i) -= rhs(i, i);
//            }
//
//            return *this;
//        }
//
//        NekMatrix<DataType, eFull, eNormal, space>& operator*=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
//        {
//            ASSERTL0(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");
//
//            NekMatrix<DataType, eFull, eNormal, space> result(GetRows(), rhs.GetColumns());
//
//            for(unsigned int i = 0; i < result.GetRows(); ++i)
//            {
//                for(unsigned int j = 0; j < result.GetColumns(); ++j)
//                {
//                    DataType t = DataType(0);
//
//                    // Set the result(i,j) element.
//                    for(unsigned int k = 0; k < GetColumns(); ++k)
//                    {
//                        t += (*this)(i,k)*rhs(k,j);
//                    }
//                    result(i,j) = t;
//                }
//            }
//
//            Swap(result);
//
//            return *this;
//        }
//
//        void Initialize(unsigned int rows, unsigned int columns)
//        {
//            m_rows = rows;
//            m_columns = columns;
//            m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
//        }
//        
//        void Initialize(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d)
//        {
//            Initialize(rows, columns);
//            std::fill(begin(), end(), d);
//        }
//
//    private:
//        void Swap(NekMatrix<DataType, eFull, eNormal, space>& rhs)
//        {
//            std::swap(m_rows, rhs.m_rows);
//            std::swap(m_columns, rhs.m_columns);
//            std::swap(m_data, rhs.m_data);
//        }
//
//        unsigned int m_rows;
//        unsigned int m_columns;
//        SharedArray<DataType> m_data;
//    };
    


}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP

/**
    $Log: NekFullMatrix.hpp,v $
    Revision 1.12  2007/04/09 03:16:53  bnelson
    *** empty log message ***

    Revision 1.11  2007/04/05 05:12:44  bnelson
    *** empty log message ***

    Revision 1.10  2007/04/04 02:26:41  bnelson
    *** empty log message ***

    Revision 1.9  2007/04/04 02:11:07  bnelson
    Added inversion

    Revision 1.8  2007/03/29 18:59:05  bnelson
    Refactoring in preparation for scaled matrices.  Fixed transpose problem.

    Revision 1.7  2007/02/15 06:56:54  bnelson
    *** empty log message ***

    Revision 1.6  2007/02/04 04:27:43  bnelson
    Updated linear systems to work with normal objects instead of only shared pointers.

    Revision 1.5  2007/01/29 01:30:20  bnelson
    Removed memory manager requirements.

    Revision 1.4  2007/01/23 03:12:50  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.3  2006/11/08 04:16:14  bnelson
    Added subtraction operators.

    Revision 1.2  2006/11/01 04:07:08  bnelson
    Changed block matrices to use the ConsistentObjectAccess object to store matrices or pointers to matrices so that the same pointer syntax works for both.

    Revision 1.1  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

 **/

