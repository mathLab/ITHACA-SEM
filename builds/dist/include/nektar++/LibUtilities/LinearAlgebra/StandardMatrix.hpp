///////////////////////////////////////////////////////////////////////////////
//
// File: StandardMatrix.hpp
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
// Description: Interface classes for matrices
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP

#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>
#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>
#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/ExpressionEvaluator.hpp>

namespace Nektar
{
    /// \brief Standard Matrix
    /// \param DataType The type stored in each element.
    ///
    /// Matrices are stored in column major order to make it easier to interoperate with 
    /// Blas and Lapack.
    template<typename DataType>
    class NekMatrix<DataType, StandardMatrixTag> : public Matrix<DataType>
    {
        public:
            typedef Matrix<DataType> BaseType;
            typedef NekMatrix<DataType, StandardMatrixTag> ThisType;
            typedef DataType NumberType;
            typedef const DataType& ConstGetValueType;
            typedef DataType GetValueType;
          
        public:
                            
            template<typename T, typename MatrixType>
            class iterator_impl
            {
                public:
                    typedef T value_type;

                    template<typename Z>
                    struct TagType
                    {
                        typedef std::forward_iterator_tag type;
                    };

                    template<typename Z>
                    struct TagType<const Z>
                    {
                        typedef std::input_iterator_tag type;
                    };

                    typedef typename TagType<T>::type iterator_category;
                    typedef unsigned int difference_type;
                    typedef typename boost::call_traits<value_type>::reference reference;
                    typedef typename boost::call_traits<value_type>::const_reference const_reference;
                    typedef typename boost::remove_reference<value_type>::type* pointer;

                public:
                    iterator_impl(pointer d, pointer e, bool isEnd = false) :
                        m_data(d),
                        m_end(e),
                        m_curRow(std::numeric_limits<unsigned int>::max()),
                        m_curColumn(std::numeric_limits<unsigned int>::max()),
                        m_matrix(NULL),
                        m_curIndex(std::numeric_limits<unsigned int>::max()),
                        m_transpose('N')
                    {
                        if( isEnd )
                        {
                            m_data = m_end;
                        }
                    }

                    iterator_impl(MatrixType* m, char transpose, bool isEnd = false) :
                        m_data(NULL),
                        m_end(NULL),
                        m_curRow(0),
                        m_curColumn(0),
                        m_matrix(m),
                        m_curIndex(0),
                        m_transpose(transpose)
                    {
                        if( isEnd )
                        {
                            m_curRow = std::numeric_limits<unsigned int>::max();
                            m_curColumn = std::numeric_limits<unsigned int>::max();
                            m_curIndex = std::numeric_limits<unsigned int>::max();
                        }
                    }

                    iterator_impl(const iterator_impl<T, MatrixType>& rhs) :
                        m_data(rhs.m_data),
                        m_end(rhs.m_end),
                        m_curRow(rhs.m_curRow),
                        m_curColumn(rhs.m_curColumn),
                        m_matrix(rhs.m_matrix),
                        m_curIndex(rhs.m_curIndex),
                        m_transpose(rhs.m_transpose)
                    {
                    }

                    iterator_impl<T, MatrixType>& operator=(const iterator_impl<T, MatrixType>& rhs)
                    {
                        m_data = rhs.m_data;
                        m_end = rhs.m_end;
                        m_curRow = rhs.m_curRow;
                        m_curColumn = rhs.m_curColumn;
                        m_matrix = rhs.m_matrix;
                        m_curIndex = rhs.m_curIndex;
                        m_transpose = rhs.m_transpose;
                        return *this;
                    }

                    reference operator*()
                    {
                        if( m_data )
                        {
                            ASSERTL1(m_data < m_end, "Attempt to dereference matrix iterator after its end.");
                            return *m_data;
                        }
                        else
                        {
                            return m_matrix->GetPtr()[m_curIndex];
                        }
                    }

                    const_reference operator*() const
                    {
                        if( m_data )
                        {
                            ASSERTL1(m_data < m_end, "Attempt to dereference matrix iterator after its end.");
                            return *m_data;
                        }
                        else
                        {
                            return m_matrix->GetPtr()[m_curIndex];
                        }
                    }

                    /// \brief Prefix increment operator.
                    iterator_impl<T, MatrixType>& operator++()
                    {
                        if( m_data )
                        {
                            ++m_data;
                        }
                        else
                        {
                            boost::tie(m_curRow, m_curColumn) =
                                m_matrix->Advance(m_curRow, m_curColumn, m_transpose);
                            if( m_curRow == std::numeric_limits<unsigned int>::max() )
                            {
                                m_curIndex = m_curRow;
                            }
                            else
                            {
                                m_curIndex = m_matrix->CalculateIndex(m_curRow, m_curColumn, m_transpose);
                            }
                        }
                        return *this;
                    }

                    /// \postfix increment operator.
                    iterator_impl<T, MatrixType> operator++(int)
                    {
                        iterator_impl<T, MatrixType> result = *this;
                        ++(*this);
                        return result;
                    }

                    bool operator==(const iterator_impl<T, MatrixType>& rhs)
                    {
                        return m_data == rhs.m_data &&
                               m_end == rhs.m_end &&
                               m_curRow == rhs.m_curRow &&
                               m_curColumn == rhs.m_curColumn &&
                               m_matrix == rhs.m_matrix &&
                               m_curIndex == rhs.m_curIndex &&
                               m_transpose == rhs.m_transpose;
                    }

                    bool operator!=(const iterator_impl<T, MatrixType>& rhs)
                    {
                        return !(*this == rhs);
                    }

                private:
                    // Used when the matrix is not transposed
                    T* m_data;
                    T* m_end;

                    // Used when the matrix is transposed.
                    unsigned int m_curRow;
                    unsigned int m_curColumn;
                    MatrixType* m_matrix;
                    unsigned int m_curIndex;
                    char m_transpose;
            };

        public:
            /// \brief Creates an empty matrix.
            LIB_UTILITIES_EXPORT NekMatrix() ;
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            LIB_UTILITIES_EXPORT NekMatrix(unsigned int rows, unsigned int columns, MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max());

            LIB_UTILITIES_EXPORT NekMatrix(unsigned int rows, unsigned int columns, 
                      MatrixStorage policy,
                      unsigned int subDiagonals,
                      unsigned int superDiagonals,
                      unsigned int capacity);
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief initValue The value used to initialize each element.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            LIB_UTILITIES_EXPORT NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max());
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief data An array of data use to initialize the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            LIB_UTILITIES_EXPORT NekMatrix(unsigned int rows, unsigned int columns, const DataType* data,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max());
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief d An array of data used to initialize the matrix.  Values from d are copied into the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            LIB_UTILITIES_EXPORT NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max());
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief d An array of data used to initialize the matrix.  If wrapperType is eCopy, then
            ///          each element is copied from d to the matrix.  If wrapperType is eWrapper, then 
            ///          the matrix uses d directly as its matrix data and no copies are made.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            LIB_UTILITIES_EXPORT NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, DataType>& d, PointerWrapper wrapperType = eCopy,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max());
            
            LIB_UTILITIES_EXPORT NekMatrix(const ThisType& rhs);
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                template<typename L, typename Op, typename R>
                NekMatrix(const expt::Node<L, Op, R>& rhs) :
                    BaseType(0, 0),
                    m_wrapperType(eCopy),
                    m_storagePolicy(eFULL),
                    m_numberOfSuperDiagonals(std::numeric_limits<unsigned int>::max()),
                    m_numberOfSubDiagonals(std::numeric_limits<unsigned int>::max()),
                    m_tempSpace()
                {
                    boost::tuple<unsigned int, unsigned int, unsigned int> sizes  =
                            MatrixSize<expt::Node<L, Op, R>, typename expt::Node<L, Op, R>::Indices, 0>::GetRequiredSize(rhs.GetData());
                            
                    unsigned int rows = sizes.get<0>();
                    unsigned int columns = sizes.get<1>();
                    unsigned int bufferCapacity = sizes.get<2>();
                    
                    m_data = Array<OneD, DataType>(bufferCapacity);
                    this->Resize(rows, columns);

                    expt::ExpressionEvaluator::EvaluateWithoutAliasingCheck(rhs, *this);
                    this->RemoveExcessCapacity();
                }

            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            
            LIB_UTILITIES_EXPORT MatrixStorage GetType() const;
           
            LIB_UTILITIES_EXPORT ThisType& operator=(const ThisType& rhs);

            template<typename InnerMatrixType>
            ThisType& operator=(const NekMatrix<InnerMatrixType, ScaledMatrixTag>& rhs)
            {
                BaseType::operator=(rhs);
                m_storagePolicy = rhs.GetType();
                m_numberOfSubDiagonals = rhs.GetNumberOfSubDiagonals();
                m_numberOfSuperDiagonals = rhs.GetNumberOfSuperDiagonals();
                
                ResizeDataArrayIfNeeded();
                
                unsigned int requiredStorageSize = GetRequiredStorageSize();
                DataType scale = rhs.Scale();
                
                DataType* lhs_array = m_data.data();
                const DataType* rhs_array = rhs.GetRawPtr();
                
                for(unsigned int i = 0; i < requiredStorageSize; ++i)
                {
                    lhs_array[i] = scale*rhs_array[i];
                }
                
                return *this;

            }
            

            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                template<typename L, typename Op, typename R>
                ThisType& operator=(const expt::Node<L, Op, R>& rhs)
                {
					if( this->GetWrapperType() == eWrapper || 
                        expt::ExpressionEvaluator::ContainsAlias(rhs, *this) )
					{
						// Assignment into a wrapped matrix can't change the matrix size (the matrix 
						// can be wrapping a portion of a larger array that can't be resized).
						// If the matrix is wrapped, we'll evaluate the expression into  temporary 
						// and then assign the temporary, relying on the operator= for the error 
						// checking relating to size.

                        // If the expression is aliased, we need to do the same thing, as 
                        // we can't change the size of the accumulator before evaluation.
						ThisType temp(rhs);
						*this = temp;
					}
					else
					{
                        // If the matrix is not wrapped, then we are free to resize as necessary.

						boost::tuple<unsigned int, unsigned int, unsigned int> sizes = 
                        MatrixSize<expt::Node<L, Op, R>, typename expt::Node<L, Op, R>::Indices, 0>::GetRequiredSize(rhs.GetData());
						unsigned int rows = sizes.get<0>();
						unsigned int columns = sizes.get<1>();
						unsigned int bufferSize = sizes.get<2>();

                        if( this->GetRows() != rows ||
							this->GetColumns() != columns )
						{
							this->Resize(rows, columns);	
						}

                        // We don't resize the temp workspace since we don't know if we will 
                        // need it.  It will be sized when needed during expression evaluation.
                        this->ResizeDataArrayIfNeeded(bufferSize);

                        this->SetTransposeFlag('N');

                        expt::ExpressionEvaluator::Evaluate(rhs, *this);
                        this->RemoveExcessCapacity();
					}

                    return *this;        
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
            
            /// \brief Returns the element value at the given row and column.
            LIB_UTILITIES_EXPORT ConstGetValueType operator()(unsigned int row, unsigned int column) const;

            class Proxy
            {
                public:
                    Proxy() : m_value(defaultReturnValue) {}                
                    explicit Proxy(DataType& value) : m_value(value) {}
                    Proxy(const Proxy& rhs) : m_value(rhs.m_value) {}
                    Proxy& operator=(const Proxy& rhs)
                    {
                        m_value = rhs.m_value;
                        return *this;
                    }
                    
                    DataType& operator*() { return m_value; }
                    const DataType& operator*() const { return m_value; }
                    operator DataType&() { return m_value; }
                    operator const DataType&() const { return m_value; }
                    void operator=(const DataType& newValue)
                    {
                        if( &m_value != &defaultReturnValue )
                        {
                            m_value = newValue;
                        }
                    }
                    
                private:
                    DataType& m_value;
                    static DataType defaultReturnValue;
            };

            /// \brief Returns the element value at the given row and column.  
            /// \brief row The element's row.
            /// \brief column The element's column.
            /// \brief transpose If transpose = 'N', then the return value is element [row, column].
            ///                  If transpose = 'T', then the return value is element [column, row].
            LIB_UTILITIES_EXPORT ConstGetValueType operator()(unsigned int row, unsigned int column, char transpose) const;
            
            LIB_UTILITIES_EXPORT unsigned int GetRequiredStorageSize() const;
            
            LIB_UTILITIES_EXPORT unsigned int CalculateIndex(unsigned int row, unsigned int col, const char transpose) const;
            
            /// \brief Returns the element value at the given row and column.
            LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference GetValue(unsigned int row, unsigned int column) const;

            /// \brief Returns the element value at the given row and column.  
            /// \brief row The element's row.
            /// \brief column The element's column.
            /// \brief transpose If transpose = 'N', then the return value is element [row, column].
            ///                  If transpose = 'T', then the return value is element [column, row].
            LIB_UTILITIES_EXPORT ConstGetValueType GetValue(unsigned int row, unsigned int column, char transpose) const;
            
            LIB_UTILITIES_EXPORT const Array<OneD, const DataType>& GetPtr() const;
            
            /// \brief Returns the scaling used by this matrix.
            ///
            /// Since this matrix is not a scaled matrix, this method always returns 1.
            LIB_UTILITIES_EXPORT DataType Scale() const;

            LIB_UTILITIES_EXPORT const DataType* GetRawPtr() const;
            
            typedef iterator_impl<const DataType, const ThisType> const_iterator;
            typedef iterator_impl<DataType, ThisType> iterator;
            
            LIB_UTILITIES_EXPORT iterator begin();

            LIB_UTILITIES_EXPORT iterator begin(char transpose);

            LIB_UTILITIES_EXPORT iterator end();
            LIB_UTILITIES_EXPORT iterator end(char transpose);

            LIB_UTILITIES_EXPORT const_iterator begin() const;

            LIB_UTILITIES_EXPORT const_iterator begin(char transpose) const;

            LIB_UTILITIES_EXPORT const_iterator end() const;

            LIB_UTILITIES_EXPORT const_iterator end(char transpose) const;
            
            
            LIB_UTILITIES_EXPORT unsigned int GetStorageSize() const;
            
             // Banded matrices only.
            // Get the specified number of sub diagonals, or calculate the 
            // number if this object has not been initialized.
            LIB_UTILITIES_EXPORT unsigned int GetNumberOfSubDiagonals() const;

            LIB_UTILITIES_EXPORT unsigned int GetNumberOfSuperDiagonals() const;
            
            LIB_UTILITIES_EXPORT unsigned int CalculateNumberOfRows() const;

            /// \brief Returns true if the this matrix and rhs are equivalent.
            ///
            /// Two matrices are equivalent if they have the same size and each element
            /// is the same.
            LIB_UTILITIES_EXPORT bool operator==(const NekMatrix<DataType, StandardMatrixTag>& rhs) const;

            LIB_UTILITIES_EXPORT PointerWrapper GetWrapperType() const;


            LIB_UTILITIES_EXPORT char GetTransposeFlag() const ;
            
            
            LIB_UTILITIES_EXPORT boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(unsigned int curRow, unsigned int curColumn) const;
            
            LIB_UTILITIES_EXPORT boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(unsigned int curRow, unsigned int curColumn, char transpose) const;
         
            LIB_UTILITIES_EXPORT static ThisType CreateWrapper(const ThisType& rhs);
            
            LIB_UTILITIES_EXPORT static boost::shared_ptr<ThisType> CreateWrapper(const boost::shared_ptr<ThisType>& rhs);

            LIB_UTILITIES_EXPORT void SetSize(unsigned int rows, unsigned int cols);

            LIB_UTILITIES_EXPORT Proxy operator()(unsigned int row, unsigned int column);
            
                        
            LIB_UTILITIES_EXPORT Proxy operator()(unsigned int row, unsigned int column, char transpose);
            
            LIB_UTILITIES_EXPORT void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d);

            LIB_UTILITIES_EXPORT void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d, char transpose);

            LIB_UTILITIES_EXPORT Array<OneD, DataType>& GetPtr();

            LIB_UTILITIES_EXPORT DataType* GetRawPtr();

            

            LIB_UTILITIES_EXPORT NekDouble AbsMaxtoMinEigenValueRatio(void);

            LIB_UTILITIES_EXPORT void EigenSolve(Array<OneD, NekDouble> &EigValReal,
                            Array<OneD, NekDouble> &EigValImag,
                            Array<OneD, NekDouble> &EigVecs =NullNekDouble1DArray);

            LIB_UTILITIES_EXPORT void Invert();
                
            LIB_UTILITIES_EXPORT Array<OneD, DataType>& GetTempSpace();
            
            LIB_UTILITIES_EXPORT void SwapTempAndDataBuffers();
                        
            LIB_UTILITIES_EXPORT ThisType& operator*=(const NumberType& s);

            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            expt::Node<expt::Node<NekMatrix<DataType, StandardMatrixTag> >, expt::NegateOp, void > operator-() const
            {
                expt::Node<NekMatrix<DataType, StandardMatrixTag> > leafNode(*this);
                return expt::Node<expt::Node<NekMatrix<DataType, StandardMatrixTag> >, expt::NegateOp, void >(leafNode);
            }
            #else
            LIB_UTILITIES_EXPORT NekMatrix<DataType, StandardMatrixTag> operator-() const;
            #endif

        protected:
            LIB_UTILITIES_EXPORT NekMatrix(const ThisType& rhs, PointerWrapper wrapperType);
            
            LIB_UTILITIES_EXPORT Array<OneD, DataType>& GetData();
            
            LIB_UTILITIES_EXPORT void RemoveExcessCapacity();
            
            LIB_UTILITIES_EXPORT void ResizeDataArrayIfNeeded(unsigned int requiredStorageSize);
            
            LIB_UTILITIES_EXPORT void ResizeDataArrayIfNeeded();
            

        private:
            LIB_UTILITIES_EXPORT void PerformCopyConstruction(const ThisType& rhs);
            
            LIB_UTILITIES_EXPORT virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const;
            
            LIB_UTILITIES_EXPORT virtual unsigned int v_GetStorageSize() const ;
            
            LIB_UTILITIES_EXPORT virtual MatrixStorage v_GetStorageType() const;
            
            // We need to rethink class structure a little.  This shouldn't be necessary.
            LIB_UTILITIES_EXPORT virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d);

            Array<OneD, DataType> m_data;
            PointerWrapper m_wrapperType;
            MatrixStorage m_storagePolicy;
            
            // Only used by banded matrices.
            unsigned int m_numberOfSuperDiagonals;
            unsigned int m_numberOfSubDiagonals;

            // Used during chained matrix multiplication as workspace.
            Array<OneD, DataType> m_tempSpace;
    };
   
    
    template<typename DataType>
    DataType NekMatrix<DataType, StandardMatrixTag>::Proxy::defaultReturnValue;

    template<typename DataType>
    LIB_UTILITIES_EXPORT NekMatrix<DataType, StandardMatrixTag>
    Transpose(NekMatrix<DataType, StandardMatrixTag>& rhs);

    template<typename DataType>
    LIB_UTILITIES_EXPORT void NegateInPlace(NekMatrix<DataType, StandardMatrixTag>& v);
}
    
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
namespace expt
{
    template<typename DataType>
    struct IsAlias<Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag>, Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag> >
    {
        static bool Apply(const Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag>& lhs, const Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag>& rhs)
        {
            return lhs.GetPtr().Overlaps(rhs.GetPtr());
        }
    };

    template<typename DataType, typename NodeType, typename Indices, unsigned int StartIndex>
    struct CreateFromTree<Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag>, NodeType, Indices, StartIndex>
    {
        template<typename ArgVectorType>
        static Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag> Apply(const ArgVectorType& tree)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> sizes = 
                Nektar::MatrixSize<NodeType, Indices, StartIndex>::GetRequiredSize(tree);

            unsigned int rows = sizes.get<0>();
            unsigned int columns = sizes.get<1>();
            unsigned int bufferSize = sizes.get<2>();

            return Nektar::NekMatrix<DataType, Nektar::StandardMatrixTag>(rows, columns,
                Nektar::eFULL, std::numeric_limits<unsigned int>::max(),
                std::numeric_limits<unsigned int>::max(), bufferSize);
        }
    };
}
#endif


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP
