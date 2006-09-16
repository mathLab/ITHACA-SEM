///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrix.hpp
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
// Description: Generic Matrix
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

#include <loki/Factory.h>
#include <loki/Singleton.h>
#include <loki/Typelist.h>

#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>
#include <boost/concept_check.hpp>
#include <boost/mpl/assert.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        enum NekMatrixForm
        {
            eFull,
            eDiagonal
//             eZero,
//             eSquareSymmetric,
//             eSquareSymmetricPositiveDefinite,
//             eSymmetricPositiveDefiniteBanded,
//             eSquareGeneral,
//             eSquareGeneralBanded
        };

        enum MatrixDataHolderType { eWrapper, eCopy };

        template<typename DataType, NekMatrixForm form, unsigned int space>
        class NekMatrixImpl;


        /// NekMatrix Class
        /// \param DataType The type of data to store in each element of the matrix.
        /// \param
        ///
        /// Some design decisions about the class.
        /// Only one factory object.  All impl objects need to implement the
        /// Initialize methods.  I could have made many constructors for the
        /// impl object constructors, but then I would have had a lot of messy
        /// Factory objects to created, one for each constructor.
        template<typename DataType, NekMatrixForm form = eFull, unsigned int space = 0>
        class NekMatrix
        {
            public:
                typedef NekMatrix<DataType, form, space> ThisType;

            public:
                NekMatrix(unsigned int rows, unsigned int columns) :
                    m_rows(rows),
                    m_columns(columns),
                    m_data(NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(rows, columns)),
                    m_dataIsDeletable(true)
                {
                }

                NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr) :
                    m_rows(rows),
                    m_columns(columns),
                    m_data(),
                    m_dataIsDeletable(true)
                {
                    m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(rows, columns);
                    std::copy(ptr, ptr+rows*columns, begin());
                }

                NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, MatrixDataHolderType t = eCopy) :
                    m_rows(rows),
                    m_columns(columns),
                    m_data(),
                    m_dataIsDeletable(false)
                {
                    if( t == eCopy )
                    {
                        m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(rows, columns);
                        std::copy(ptr, ptr+rows*columns, begin());
                        m_dataIsDeletable = true;
                    }
                    else
                    {
                        m_data = ptr;
                    }
                }

                NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d) :
                        m_rows(rows),
                        m_columns(columns),
                        m_data(NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(rows, columns)),
                        m_dataIsDeletable(true)
                {
                    std::fill(begin(), end(), d);
                }


                NekMatrix(const NekMatrix<DataType, form, space>& rhs) :
                    m_rows(rhs.m_rows),
                    m_columns(rhs.m_columns),
                    m_data(),
                    m_dataIsDeletable(rhs.m_dataIsDeletable)
                {
                    if( m_dataIsDeletable )
                    {
                        m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(m_rows, m_columns);
                        std::copy(rhs.begin(), rhs.end(), begin());
                    }
                    else
                    {
                        m_data = rhs.m_data;
                    }
                }

                template<NekMatrixForm rhsForm>
                NekMatrix(const NekMatrix<DataType, rhsForm, space>& rhs) :
                    m_rows(rhs.m_rows),
                    m_columns(rhs.m_columns),
                    m_data(),
                    m_dataIsDeletable(true)
                {
                
                    m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(m_rows, m_columns);
                    NekMatrixImpl<DataType, form, space>::CopyMatrixValues(m_data, rhs, m_rows, m_columns);
                }

                template<typename ExpressionPolicyType>
                NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                    m_rows(rhs.GetMetadata().Rows),
                    m_columns(rhs.GetMetadata().Columns),
                    m_data(),
                    m_dataIsDeletable(true)
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, form, space> > ));
                    m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(m_rows, m_columns);
                    rhs.Apply(*this);
                }

                NekMatrix<DataType, form, space>& operator=(const NekMatrix<DataType, form, space>& rhs)
                {
                    NekMatrix<DataType, form, space> temp(rhs);
                    Swap(temp);
                    return *this;
                }

                template<NekMatrixForm rhsForm>
                NekMatrix<DataType, form, space>& operator=(const NekMatrix<DataType, rhsForm, space>& rhs)
                {
                    m_rows = rhs.GetRows();
                    m_columns = rhs.GetColumns();
                    m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(m_rows, m_columns);
                    m_dataIsDeletable = true;
                    NekMatrixImpl<DataType, form, space>::CopyMatrixValues(m_data, rhs, m_rows, m_columns);
                    return *this;
                }

                template<typename ExpressionPolicyType>
                NekMatrix<DataType, form, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, form, space> > ));
                    m_rows = rhs.GetMetadata().Rows;
                    m_columns = rhs.GetMetadata().Columns;
                    m_data = NekMatrixImpl<DataType, form, space>::CreateMatrixStorage(m_rows, m_columns);
                    m_dataIsDeletable = true;
                    rhs.Apply(*this);
                    return *this;
                }

                ~NekMatrix()
                {
                    if( m_dataIsDeletable )
                    {
                        NekMatrixImpl<DataType, form, space>::DeallocateMatrixStorage(m_data, m_rows, m_columns);
                    }
                    m_data = NULL;
                }

                unsigned int GetRows() const
                {
                    return m_rows;
                }

                unsigned int GetColumns() const
                {
                    return m_columns;
                }

                inline typename boost::call_traits<DataType>::reference operator()
                    (unsigned int rowNumber, unsigned int colNumber)
                {
                    return NekMatrixImpl<DataType, form, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
                }

                inline typename boost::call_traits<DataType>::const_reference operator()
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    return NekMatrixImpl<DataType, form, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
                }


                typename boost::call_traits<DataType>::reference GetValue
                        (unsigned int rowNumber, unsigned int colNumber)
                {
                    return NekMatrixImpl<DataType, form, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
                }

                typename boost::call_traits<DataType>::const_reference GetValue
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    return NekMatrixImpl<DataType, form, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
                }

                void SetValue(unsigned int rowNumber, unsigned int colNumber, typename boost::call_traits<DataType>::param_type rhs)
                {
                    NekMatrixImpl<DataType, form, space>::SetValue(m_data, m_columns, rowNumber, colNumber, rhs);
                }

                DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
                {
                    return NekMatrixImpl<DataType, form, space>::GetPtr(m_data, m_columns, rowNumber, colNumber);
                }

                typedef DataType* iterator;
                typedef const DataType* const_iterator;

                iterator begin()
                {
                    return &m_data[0];
                }

                iterator end()
                {
                    return NekMatrixImpl<DataType, form, space>::end(m_data, m_rows, m_columns);
                }

                const_iterator begin() const
                {
                    return &m_data[0];
                }

                const_iterator end() const
                {
                    return NekMatrixImpl<DataType, form, space>::end(m_data, m_rows, m_columns);
                }

                NekMatrixForm GetForm() const
                {
                    return form;
                }

                void negate()
                {
                    for(iterator iter = begin(); iter != end(); ++iter)
                    {
                        *iter = -*iter;
                    }
                }

                expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, form, space> > >, expt::NegateOp> > operator-() const
                {
                    return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, form, space> > >, expt::NegateOp> >(
                            expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, form, space> > >(*this));
                }


                NekMatrix<DataType, eFull, space> operator+=(const NekMatrix<DataType, eFull, space>& rhs)
                {
                    ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");
                    DataType* lhs_data = begin();
                    const DataType* rhs_data = rhs.begin();

                    for(unsigned int i = 0; i < m_rows*m_columns; ++i)
                    {
                        lhs_data[i] += rhs_data[i];
                    }

                    return *this;
                }

                NekMatrix<DataType, eFull, space> operator+=(const NekMatrix<DataType, eDiagonal, space>& rhs)
                {
                    ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");

                    for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                    {
                        (*this)(i,i) += rhs(i,i);
                    }

                    return *this;
                }

            private:


                void Swap(NekMatrix<DataType, form, space>& rhs)
                {
                    std::swap(m_rows, rhs.m_rows);
                    std::swap(m_columns, rhs.m_columns);
                    std::swap(m_data, rhs.m_data);
                    std::swap(m_dataIsDeletable, rhs.m_dataIsDeletable);
                }

                unsigned int m_rows;
                unsigned int m_columns;
                DataType* m_data;
                bool m_dataIsDeletable;
        };

        template<typename DataType, unsigned int space>
        class NekMatrixImpl<DataType, eFull, space>
        {
            public:
                static void CopyMatrixValues(DataType* data, const NekMatrix<DataType, eDiagonal, space>& rhs, unsigned int rows, unsigned int columns)
                {
                    std::fill(data, end(data, rows, columns), DataType(0));
                    for(unsigned int i = 0; i < rows; ++i)
                    {
                        SetValue(data, columns, i, i, rhs.GetValue(i,i));
                    }
                }

                static inline DataType* CreateMatrixStorage(unsigned int rows, unsigned int columns)
                {
                    return Nektar::MemoryManager::AllocateArray<DataType>(rows*columns);
                }

                static inline void DeallocateMatrixStorage(DataType*& data, unsigned int rows, unsigned int columns)
                {
                    Nektar::MemoryManager::DeallocateArray<DataType>(data, rows*columns);
                }

                static inline typename NekMatrix<DataType, eFull, space>::iterator end(DataType* data, unsigned int rows, unsigned int columns)
                {
                    return &data[rows*columns];
                }

                static inline typename boost::call_traits<DataType>::reference GetValue(
                DataType* data, unsigned int matrixColumns, unsigned int rowNumber, unsigned int colNumber)
                {
                    return data[rowNumber*matrixColumns + colNumber];
                }

                static inline void SetValue(DataType* data, unsigned int matrixColumns, unsigned int rowNumber,
                        unsigned int colNumber, typename boost::call_traits<DataType>::param_type rhs)
                {
                    NekMatrixImpl<DataType, eFull, space>::GetValue(data, matrixColumns, rowNumber, colNumber) = rhs;
                }

                static inline DataType* GetPtr(DataType* data, unsigned int matrixColumns, unsigned int rowNumber, unsigned int colNumber)
                {
                    return &data[rowNumber*matrixColumns + colNumber];
                }
        };

        template<typename DataType, unsigned int space>
        class NekMatrixImpl<DataType, eDiagonal, space>
        {
            public:
                static inline DataType* CreateMatrixStorage(unsigned int rows, unsigned int columns)
                {
                    ASSERTL0(rows == columns, "Digaonal matrices must be square.");
                    return Nektar::MemoryManager::AllocateArray<DataType>(rows);
                }

                static inline void DeallocateMatrixStorage(DataType*& data, unsigned int rows, unsigned int /*columns*/)
                {
                    Nektar::MemoryManager::DeallocateArray<DataType>(data, rows);
                }

                static inline typename NekMatrix<DataType, eDiagonal, space>::iterator end(DataType* data, unsigned int rows, unsigned int /*columns*/)
                {
                    return &data[rows];
                }

                static inline typename boost::call_traits<DataType>::reference GetValue(
                DataType* data, unsigned int /*matrixColumns*/, unsigned int rowNumber, unsigned int /*colNumber*/)
                {
                    return data[rowNumber];
                }

                static inline void SetValue(DataType* data, unsigned int matrixColumns, unsigned int rowNumber,
                        unsigned int colNumber, typename boost::call_traits<DataType>::param_type rhs)
                {
                    NekMatrixImpl<DataType, eFull, space>::GetValue(data, matrixColumns, rowNumber, colNumber) = rhs;
                }

                static inline DataType* GetPtr(DataType* data, unsigned int /*matrixColumns*/, unsigned int rowNumber, unsigned int /*colNumber*/)
                {
                    return &data[rowNumber];
                }
        };


        template<typename DataType, unsigned int space>
        const NekMatrix<DataType, eFull, space>& convertToFull(const NekMatrix<DataType, eFull, space>& m)
        {
            return m;
        }

        template<typename DataType, unsigned int space>
        NekMatrix<DataType, eFull, space> convertToFull(const NekMatrix<DataType, eDiagonal, space>& m)
        {
            NekMatrix<DataType, eFull, space> result(m.GetRows(), m.GetColumns(), DataType());
            for(unsigned int i = 0; i < m.GetRows(); ++i)
            {
                result(i,i) = m(i,i);
            }

            return result;
        }
    }

    // All of the expression interfaces for NekMatrix should go here.
    namespace expt
    {
        template<typename DataType, LibUtilities::NekMatrixForm form, unsigned int space>
        class ExpressionTraits<LibUtilities::NekMatrix<DataType, form, space> >
        {
            public:
                typedef LibUtilities::NekMatrixMetadata MetadataType;
        };

        // Binary expression specializations for NekMatrix.
        template<typename DataType, Nektar::LibUtilities::NekMatrixForm lhsForm, Nektar::LibUtilities::NekMatrixForm rhsForm, unsigned int space>
        class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, lhsForm, space>, Nektar::LibUtilities::NekMatrix<DataType, rhsForm, space> >
        {
            public:
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> AdditionResultType;
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> SubtractionResultType;
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> DivisionResultType;
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> MultiplicationResultType;
        };

        //// Binary expression specializations for NekMatrix.
        //template<typename DataType, Nektar::LibUtilities::NekMatrixForm rhsForm, unsigned int space>
        //class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space>, Nektar::LibUtilities::NekMatrix<DataType, rhsForm, space> >
        //{
        //    public:
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> AdditionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> SubtractionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> DivisionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> MultiplicationResultType;
        //};

        //// Binary expression specializations for NekMatrix.
        //template<typename DataType, Nektar::LibUtilities::NekMatrixForm lhsForm, unsigned int space>
        //class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, lhsForm, space>, Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> >
        //{
        //    public:
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> AdditionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> SubtractionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> DivisionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> MultiplicationResultType;
        //};

        // Binary expression specializations for NekMatrix.
        template<typename DataType, unsigned int space>
        class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eDiagonal, space>, Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eDiagonal, space> >
        {
            public:
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eDiagonal, space> AdditionResultType;
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eDiagonal, space> SubtractionResultType;
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eDiagonal, space> DivisionResultType;
                typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eDiagonal, space> MultiplicationResultType;
        };
        
        
        template<typename DataType, Nektar::LibUtilities::NekMatrixForm lhsForm, unsigned int vectorDim, unsigned int space>
        class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, lhsForm, space>, Nektar::LibUtilities::NekVector<DataType, vectorDim, space> >
        {
            public:
                typedef Nektar::LibUtilities::NekVector<DataType, vectorDim, space> MultiplicationResultType;
        };

    }

    namespace LibUtilities
    {
        template<typename DataType, NekMatrixForm form, unsigned int space>
        void negate(NekMatrix<DataType, form, space>& rhs)
        {
            rhs.negate();
        }
        
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, unsigned int space>
        expt::Expression<expt::BinaryExpressionPolicy<
            expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, space> > >,
            expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, space> > >,
            expt::AddOp > > operator+(
            const NekMatrix<DataType, lhsForm, space>& lhs,
            const NekMatrix<DataType, rhsForm, space>& rhs)
        {
            typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, space> > > LhsExpressionType;
            typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, space> > > RhsExpressionType;

            return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::AddOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
        }
#else
        template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, unsigned int space>
        NekMatrix<DataType, eFull, space> operator+(
            const NekMatrix<DataType, lhsForm, space>& lhs,
            const NekMatrix<DataType, rhsForm, space>& rhs)
        {
            NekMatrix<DataType, NekMatrixOperationResult<lhsForm, rhsForm>::AdditionResultType, space> result(lhs);
            result += rhs;
            return result;
        }

        template<typename DataType, unsigned int space>
        NekMatrix<DataType, eFull, space> operator+(
            const NekMatrix<DataType, eFull, space>& lhs,
            const NekMatrix<DataType, eDiagonal, space>& rhs)
        {
            NekMatrix<DataType, eFull, space> result = lhs;
            result += rhs;
            return result;
        }

        template<typename DataType, unsigned int space>
        NekMatrix<DataType, eFull, space> operator+(
            const NekMatrix<DataType, eDiagonal, space>& lhs,
            const NekMatrix<DataType, eFull, space>& rhs)
        {
            NekMatrix<DataType, eFull, space> result = rhs;
            result += lhs;
            return result;
        }
#endif

//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//        template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, unsigned int space>
//        Expression<BinaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, space> > >, Expression<ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, space> > >, MultiplyOp > > operator*(
//                const NekMatrix<DataType, lhsForm, space>& lhs,
//                const NekMatrix<DataType, rhsForm, space>& rhs)
//        {
//            typedef Expression<ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, space> > > LhsExpressionType;
//            typedef Expression<ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, space> > > RhsExpressionType;
//
//            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//        }
//
//        template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, unsigned int space>
//        Expression<BinaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, space> > >, Expression<ConstantExpressionPolicy<NekVector<DataType, vectorDim, space> > >, MultiplyOp > > operator*(
//                const NekMatrix<DataType, lhsForm, space>& lhs,
//                const NekVector<DataType, vectorDim, space>& rhs)
//        {
//            typedef Expression<ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, space> > > LhsExpressionType;
//            typedef Expression<ConstantExpressionPolicy<NekVector<DataType, vectorDim, space> > > RhsExpressionType;
//
//            return Expression<BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//        }
//#else
        template<typename DataType, NekMatrixForm form, unsigned int space>
        NekMatrix<DataType, form, space> operator*(
        const NekMatrix<DataType, form, space>& lhs,
        const NekMatrix<DataType, form, space>& rhs)
        {
            ASSERTL0(lhs.GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");

            NekMatrix<DataType, form, space> result(lhs.GetRows(), rhs.GetColumns());

            for(unsigned int i = 0; i < result.GetRows(); ++i)
            {
                for(unsigned int j = 0; j < result.GetColumns(); ++j)
                {
                    DataType t = DataType(0);

                    // Set the result(i,j) element.
                    for(unsigned int k = 0; k < lhs.GetColumns(); ++k)
                    {
                        t += lhs(i,k)*rhs(k,j);
                    }
                    result(i,j) = t;
                }
            }

            return result;
        }

        template<typename DataType, NekMatrixForm form, unsigned int space>
        NekVector<DataType, 0, space> operator*(
        const NekMatrix<DataType, form, space>& lhs,
        const NekVector<DataType, 0, space>& rhs)
        {
            ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");

            NekVector<DataType, 0, space> result(rhs.GetDimension(), DataType(0));

            for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
            {
                DataType t = DataType(0);
                for(unsigned int j = 0; j < rhs.GetDimension(); ++j)
                {
                    t += lhs(i,j)*rhs(j);
                }
                result(i) = t;
            }

            return result;
        }
//#endif


        template<typename DataType, NekMatrixForm form, unsigned int space>
        bool operator==(const NekMatrix<DataType, form, space>& lhs,
                        const NekMatrix<DataType, form, space>& rhs)
        {
            if( lhs.GetRows() != rhs.GetRows() )
            {
                return false;
            }

            if( lhs.GetColumns() != rhs.GetColumns() )
            {
                return false;
            }

            typename NekMatrix<DataType, form, space>::const_iterator lhs_iter = lhs.begin();
            typename NekMatrix<DataType, form, space>::const_iterator rhs_iter = rhs.begin();

            for( ; lhs_iter != lhs.end(); ++lhs_iter, ++rhs_iter )
            {
                if( *lhs_iter != *rhs_iter )
                {
                    return false;
                }
            }

            return true;
        }

        template<typename DataType, NekMatrixForm form, unsigned int space>
        std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, form, space>& rhs)
        {
            for(unsigned int i = 0; i < rhs.GetRows(); ++i)
            {
                os << "[";
                for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
                {
                    os << rhs(i,j);
                    if( j != rhs.GetColumns() - 1 )
                    {
                        os << ", ";
                    }
                }
                os << "]";
                if( i != rhs.GetRows()-1 )
                {
                    os << std::endl;
                }
            }
            return os;
        }
    }

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

/**
    $Log: NekMatrix.hpp,v $
    Revision 1.8  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.7  2006/09/11 03:26:26  bnelson
    Updated to use new policy based expression templates.

    Revision 1.6  2006/08/28 02:40:21  bnelson
    *** empty log message ***

    Revision 1.5  2006/08/25 03:05:16  bnelson
    Fixed gcc compile errors.

    Revision 1.4  2006/08/25 01:28:53  bnelson
    Changed the way specialized matrices are handled.

    Added support for a matrix wrapper around a pre-allocated array of data.

    Added the NekMemoryManager for allocations.

    Revision 1.3  2006/08/14 02:29:49  bnelson
    Updated points, vectors, and matrix classes to work with ElVis.  Added a variety of methods to all of these classes.

    Revision 1.2  2006/06/01 13:44:28  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:12:41  kirby
    *** empty log message ***

    Revision 1.11  2006/05/31 23:24:21  bnelson
    Moved the diagonal matrix to NekMatrix.hpp.

    Updated method names for the coding standard.

    Revision 1.10  2006/05/31 04:20:16  bnelson
    Changed matrix implementation so the form is a template parameter.

    Revision 1.9  2006/05/29 04:32:18  bnelson
    Changed the data holder to boost::shared_array.

    Revision 1.8  2006/05/29 03:45:04  bnelson
    Updated operator+= to be more efficient using iterators.

    Revision 1.7  2006/05/29 03:40:12  bnelson
    Changed the implementation from individual NekMatrixImpl objects for each type of matrix to an enumeration to store the type.

    Revision 1.6  2006/05/25 03:02:40  bnelson
    Added Matrix/Vector multiplication.

    Revision 1.5  2006/05/18 04:21:06  bnelson
    Removed the ability to specify the rows and columns in the template parameter list.  If this behavior is desired we'll need to create a fixed size array class.

    Added multiplication to the arrays.

    Revision 1.4  2006/05/15 05:06:55  bnelson
    Removed use of MemoryManager pending review of some memory problems.

    Revision 1.3  2006/05/15 04:13:36  bnelson
    no message

    Revision 1.2  2006/05/14 21:32:03  bnelson
 *** empty log message ***

    Revision 1.1  2006/05/04 18:57:43  kirby
 *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


 **/


