///////////////////////////////////////////////////////////////////////////////
//
// File: NekDiagonalMatrix.hpp
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
// Description: Diagonal matrix
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_DIAGONAL_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_DIAGONAL_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/Memory/DeleteNothing.hpp>

#include <boost/call_traits.hpp>

#include <algorithm>

namespace Nektar
{
    template<typename DataType, unsigned int space>
    class NekMatrix<DataType, eDiagonal, eNormal, space>
    {
        public:
            typedef NekMatrix<DataType, eDiagonal, eNormal, space> ThisType;


            /// Necessary to allow construction of default constructed matrices.
            friend class MemoryManager;
            static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eEnabled;

        public:

            NekMatrix() :
                m_numberOfElements(0),
                m_data()
            {
            }
            
            /// \brief Create an identity matrix.
            ///
            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
            explicit NekMatrix(unsigned int numberOfElements) :
                m_numberOfElements(numberOfElements),
                m_data(MemoryManager::AllocateSharedArray<DataType>(numberOfElements))
            {
                Initialize(DataType(1));
            }

            NekMatrix(unsigned int numberOfElements, const DataType* const ptr) :
                m_numberOfElements(numberOfElements),
                m_data(MemoryManager::AllocateSharedArray<DataType>(numberOfElements))
            {
                std::copy(ptr, ptr+m_numberOfElements, begin());
            }

            NekMatrix(unsigned int numberOfElements, DataType* ptr, MatrixDataHolderType t = eCopy) :
                m_numberOfElements(numberOfElements),
                m_data()
            {
                if( t == eCopy )
                {
                    m_data = MemoryManager::AllocateSharedArray<DataType>(numberOfElements);
                    std::copy(ptr, ptr+m_numberOfElements, begin());
                }
                else
                {
                    m_data = boost::shared_array<DataType>(ptr, DeleteNothing<DataType>());
                }
            }

            NekMatrix(unsigned int numberOfElements, typename boost::call_traits<DataType>::const_reference d) :
                m_numberOfElements(numberOfElements),
                m_data(MemoryManager::AllocateSharedArray<DataType>(numberOfElements))
            {
                Initialize(d);
            }


            NekMatrix(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs) :
                m_numberOfElements(rhs.m_numberOfElements),
                m_data(MemoryManager::AllocateSharedArray<DataType>(m_numberOfElements))
            {
                std::copy(rhs.begin(), rhs.end(), begin());
            }


            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_numberOfElements(rhs.GetMetadata().Rows),
                m_data(MemoryManager::AllocateSharedArray<DataType>(m_numberOfElements))
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eDiagonal, eNormal, space> > ));
                rhs.Apply(*this);
            }

            NekMatrix<DataType, eDiagonal, eNormal, space>& operator=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                NekMatrix<DataType, eDiagonal, eNormal, space> temp(rhs);
                Swap(temp);
                return *this;
            }

            template<typename ExpressionPolicyType>
            NekMatrix<DataType, eDiagonal, eNormal, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eDiagonal, eNormal, space> > ));
                
                m_numberOfElements = rhs.GetMetadata().Rows;
                m_data = MemoryManager::AllocateSharedArray<DataType>(m_numberOfElements);
                rhs.Apply(*this);
                return *this;
            }

            ~NekMatrix() {}

            unsigned int GetRows() const { return m_numberOfElements; }
            unsigned int GetColumns() const { return m_numberOfElements; }

            typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
            {
                if( rowNumber == colNumber )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return m_data[rowNumber];
                }
                else
                {
                    static DataType zeroElement(0);
                    return zeroElement;
                }
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
            {
                if( rowNumber == colNumber )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return m_data[rowNumber];
                }
                else
                {
                    static DataType zeroElement(0);
                    return zeroElement;
                }
            }

            DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
            {
                ASSERTL1(rowNumber == colNumber, "NekMatrix::GetPtr requires equal row an column numbers for diagonal matrices.");
                return m_data.get() + rowNumber;
            }

            boost::shared_ptr<DataType> GetPtr()
            {
                return m_data.get();
            }

            boost::shared_ptr<const DataType> GetPtr() const
            {
                return m_data.get();
            }

            typedef DataType* iterator;
            typedef const DataType* const_iterator;

            iterator begin()
            {
                return m_data.get();
            }

            iterator end()
            {
                return m_data.get() + m_numberOfElements;
            }

            const_iterator begin() const
            {
                return m_data.get();
            }

            const_iterator end() const
            {
                return m_data.get() + m_numberOfElements;
            }

            void Negate()
            {
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    m_data[i] = -m_data[i];
                }
            }
            
            void Transpose()
            {
                // A diagonal matrix is its own transpose.
            }

            expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eDiagonal, eNormal, space> > >, expt::NegateOp> > operator-() const
            {
                return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eDiagonal, eNormal, space> > >, expt::NegateOp> >(
                        expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eDiagonal, eNormal, space> > >(*this));
            }


//             NekMatrix<DataType, eFull, BlockType, space> operator+=(const NekMatrix<DataType, eFull, BlockType, space>& rhs)
//             {
//                 ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");
//                 DataType* lhs_data = begin();
//                 const DataType* rhs_data = rhs.begin();
// 
//                 for(unsigned int i = 0; i < m_rows*m_columns; ++i)
//                 {
//                     lhs_data[i] += rhs_data[i];
//                 }
// 
//                 return *this;
//             }

            // This is wrong as well.  What if this is diagonal?
            // enable if on the output.
            NekMatrix<DataType, eDiagonal, eNormal, space> operator+=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");

                for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                {
                    (*this)(i,i) += rhs(i,i);
                }

                return *this;
            }
            
            NekMatrix<DataType, eDiagonal, eNormal, space> operator-=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");

                for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                {
                    (*this)(i,i) -= rhs(i,i);
                }

                return *this;
            }

            // Full *= full = full
            // Full *= diagonal = full

            // diag *= diag = diag
            // diag *= full = full


//             template<NekMatrixForm rhsForm>
//             NekMatrix<DataType, eFull, BlockType, space>& operator*=(const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//             {
//                 ASSERTL0(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");
// 
//                 NekMatrix<DataType, eFull, BlockType, space> result(GetRows(), rhs.GetColumns());
// 
//                 for(unsigned int i = 0; i < result.GetRows(); ++i)
//                 {
//                     for(unsigned int j = 0; j < result.GetColumns(); ++j)
//                     {
//                         DataType t = DataType(0);
// 
//                         // Set the result(i,j) element.
//                         for(unsigned int k = 0; k < GetColumns(); ++k)
//                         {
//                             t += (*this)(i,k)*rhs(k,j);
//                         }
//                         result(i,j) = t;
//                     }
//                 }
// 
//                 Swap(result);
// 
//                 return *this;
//             }

            NekMatrix<DataType, eDiagonal, eNormal, space>& operator*=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator*=");

                for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                {
                    (*this)(i,i) *= rhs(i,i);
                }

                return *this;
            }
            
            void Initialize(typename boost::call_traits<DataType>::const_reference d)
            {
                std::fill(begin(), end(), d);
            }
            
        private:
            /// \brief Constructor used for block matrices.  Care must be used with this constructor
            ///         to initialize the matrix before use.
            //NekMatrix();
            
            //void Initialize(unsigned int rows, unsigned int columns);
            
            void Swap(NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                std::swap(m_numberOfElements, rhs.m_numberOfElements);
                std::swap(m_data, rhs.m_data);
            }

            unsigned int m_numberOfElements;
            boost::shared_array<DataType> m_data;

    };
            
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_DIAGONAL_MATRIX_HPP

/**
    $Log: NekDiagonalMatrix.hpp,v $
    Revision 1.2  2006/11/06 17:09:09  bnelson
    *** empty log message ***

    Revision 1.1  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

 **/
 
