///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixBase.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP

#include <LibUtilities/LinearAlgebra/MatrixType.h>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <boost/lexical_cast.hpp>

#ifdef max
#undef max
#endif

namespace Nektar
{    
    template<typename DataType>
    class ConstMatrix
    {
        public:
            typedef ConstMatrix<DataType> ThisType;
            virtual ~ConstMatrix() {}
            
        public:
            typename boost::call_traits<DataType>::value_type operator()(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetColumns()) +
                    std::string(" columns"));
                return v_GetValue(row, column);
            }
            
            unsigned int GetStorageSize() const
            {
                return v_GetStorageSize();
            }
            
            MatrixStorage GetStorageType() const 
            {
                return v_GetStorageType();
            }
            
            unsigned int GetRows() const
            {
                return GetTransposedRows(GetTransposeFlag());
            }

            unsigned int GetTransposedRows(char transpose) const 
            {
                if( transpose == 'N' )
                {
                    return m_size[0]; 
                }
                else
                {
                    return m_size[1];
                }
            }
            
            unsigned int GetColumns() const
            {
                return GetTransposedColumns(GetTransposeFlag());
            }

            unsigned int GetTransposedColumns(char transpose) const 
            { 
                if( transpose == 'N' )
                {
                    return m_size[1]; 
                }
                else
                {
                    return m_size[0];
                }
            }

            const unsigned int* GetSize() const { return m_size; }
            
            void Transpose() 
            {
                if( m_transpose == 'N' )
                {
                    m_transpose = 'T';
                }
                else
                {
                    m_transpose = 'N';
                }
                v_Transpose();
            }

            char GetTransposeFlag() const 
            {
                return v_GetTransposeFlag();
            }       
            
            static unsigned int CalculateIndex(MatrixStorage type, 
                unsigned int row, unsigned int col, 
                unsigned int numRows, unsigned int numColumns, const char transpose =  'N',
                unsigned int numSubDiags = 0, unsigned int numSuperDiags = 0) 
            {
                if(transpose == 'T' )
                {
                    std::swap(row, col);
                }
                switch(type)
                {
                    case eFULL:
                        return FullMatrixFuncs::CalculateIndex(numRows, numColumns, row, col);
                        break;
                    case eDIAGONAL:
                        return DiagonalMatrixFuncs::CalculateIndex(row, col);
                        break;
                    case eUPPER_TRIANGULAR:
                        return UpperTriangularMatrixFuncs::CalculateIndex(row, col);                        
                        break;
                    case eLOWER_TRIANGULAR:
                        return LowerTriangularMatrixFuncs::CalculateIndex(numColumns, row, col);
                        break;
                    case eSYMMETRIC:
                    case ePOSITIVE_DEFINITE_SYMMETRIC:
                        return SymmetricMatrixFuncs::CalculateIndex(row, col);
                        break;
                    case eBANDED:
                        return BandedMatrixFuncs::CalculateIndex(numRows, numColumns,
                            row, col, numSubDiags, numSuperDiags);
                        break;
                    case eSYMMETRIC_BANDED:
                    case ePOSITIVE_DEFINITE_SYMMETRIC_BANDED:
                        {
                            ASSERTL1(numSubDiags==numSuperDiags,
                                     std::string("Number of sub- and superdiagonals should ") + 
                                     std::string("be equal for a symmetric banded matrix"));
                            return SymmetricBandedMatrixFuncs::CalculateIndex(row, col, 
                                numSuperDiags);
                        }
                        break;
                    case eUPPER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Not yet implemented.");
                        break;
                    case eLOWER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Not yet implemented.");
                        break;
                        
                    default:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                }
                
                return std::numeric_limits<unsigned int>::max();
            }
            
            static unsigned int GetRequiredStorageSize(MatrixStorage type, unsigned int rows, 
                unsigned int columns, unsigned int subDiags = 0, unsigned int superDiags = 0)
            {
                switch(type)
                {
                    case eFULL:
                        return FullMatrixFuncs::GetRequiredStorageSize(rows, columns);
                        break;
                    case eDIAGONAL:
                        return DiagonalMatrixFuncs::GetRequiredStorageSize(rows, columns);
                        break;
                    case eUPPER_TRIANGULAR:
                        return UpperTriangularMatrixFuncs::GetRequiredStorageSize(rows, columns);
                        break;
                    case eLOWER_TRIANGULAR:
                        return LowerTriangularMatrixFuncs::GetRequiredStorageSize(rows, columns);
                        break;
                    case eSYMMETRIC:
                    case ePOSITIVE_DEFINITE_SYMMETRIC:
                        return SymmetricMatrixFuncs::GetRequiredStorageSize(rows, columns);
                        break;
                    case eBANDED:
                        return BandedMatrixFuncs::GetRequiredStorageSize(rows, columns,
                            subDiags, superDiags);
                        break;
                    case eSYMMETRIC_BANDED:
                    case ePOSITIVE_DEFINITE_SYMMETRIC_BANDED:
                        {
                            ASSERTL1(subDiags==superDiags,
                                     std::string("Number of sub- and superdiagonals should ") + 
                                     std::string("be equal for a symmetric banded matrix"));
                            return SymmetricBandedMatrixFuncs::GetRequiredStorageSize(rows, columns,
                                superDiags);
                        }
                        break;
                    case eUPPER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eLOWER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                        
                    default:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                }

                return 0;
            }

        protected:
            
            // All constructors are private to enforce the abstract nature of ConstMatrix without
            // resorting to pure virtual functions.
            ConstMatrix(unsigned int rows, unsigned int columns) :
                m_size(),
                m_transpose('N')
            {
                m_size[0] = rows;
                m_size[1] = columns;
            }
            
            ConstMatrix(const ThisType& rhs) :
                m_size(),
                m_transpose(rhs.m_transpose)
            {
                m_size[0] = rhs.m_size[0];
                m_size[1] = rhs.m_size[1];
            }

            ThisType& operator=(const ThisType& rhs)
            {
                m_size[0] = rhs.m_size[0];
                m_size[1] = rhs.m_size[1];
                m_transpose = rhs.m_transpose;
                return *this;
            }
            
            /// \brief Resets the rows and columns in the array.
            /// This method does not update the data storage to match the new row and column counts.
            void Resize(unsigned int rows, unsigned int columns)
            {
                m_size[0] = rows;
                m_size[1] = columns;
            }

            void SetTransposeFlag(char newValue)
            {
                m_transpose = newValue; 
            }

            char GetRawTransposeFlag() const { return m_transpose; }

        private:
            
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const = 0;            
            virtual unsigned int v_GetStorageSize() const = 0;            
            virtual MatrixStorage v_GetStorageType() const = 0;
            virtual void v_Transpose() {};
            virtual char v_GetTransposeFlag() const { return m_transpose; }
            unsigned int m_size[2];
            char m_transpose;
    };
    
    template<typename DataType>
    class Matrix : public ConstMatrix<DataType>
    {
        public:
            typedef Matrix<DataType> ThisType;
            typedef ConstMatrix<DataType> BaseType;
            
        public:
            virtual ~Matrix() {}
            
            void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                v_SetValue(row, column, d);
            }
            
        protected:            
            // All constructors are private to enforce the abstract nature of ConstMatrix without
            // resorting to pure virtual functions.
            Matrix(unsigned int rows, unsigned int columns) :
                BaseType(rows, columns)
            {
            }
            
            Matrix(const ThisType& rhs) :
                BaseType(rhs)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
            
            ThisType& operator=(const ConstMatrix<DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
            
        private:
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d) = 0;
    };
    
    
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_BASE_HPP
