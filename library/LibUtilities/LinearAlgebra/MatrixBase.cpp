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

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>

namespace Nektar
{

    template<typename DataType>
    ConstMatrix<DataType>::~ConstMatrix() {}

    template<typename DataType>
    typename boost::call_traits<DataType>::value_type ConstMatrix<DataType>::operator()(unsigned int row, unsigned int column) const
    {
        ASSERTL2(row < GetRows(), std::string("Row ") + std::to_string(row) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(GetRows()) +
            std::string(" rows"));
        ASSERTL2(column < GetColumns(), std::string("Column ") + std::to_string(column) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(GetColumns()) +
            std::string(" columns"));
        return v_GetValue(row, column);
    }

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::GetStorageSize() const
    {
        return v_GetStorageSize();
    }

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::GetRows() const
    {
        return GetTransposedRows(GetTransposeFlag());
    }

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::GetTransposedRows(char transpose) const
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

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::GetColumns() const
    {
        return GetTransposedColumns(GetTransposeFlag());
    }

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::GetTransposedColumns(char transpose) const
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

    template<typename DataType>
    const unsigned int* ConstMatrix<DataType>::GetSize() const { return m_size; }

    template<typename DataType>
    void ConstMatrix<DataType>::Transpose()
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

    template<typename DataType>
    char ConstMatrix<DataType>::GetTransposeFlag() const
    {
        return v_GetTransposeFlag();
    }

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::CalculateIndex(MatrixStorage type,
        unsigned int row, unsigned int col,
        unsigned int numRows, unsigned int numColumns, const char transpose,
        unsigned int numSubDiags, unsigned int numSuperDiags )
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

    template<typename DataType>
    unsigned int ConstMatrix<DataType>::GetRequiredStorageSize(MatrixStorage type, unsigned int rows,
        unsigned int columns, unsigned int subDiags, unsigned int superDiags)
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

    template<typename DataType>
    ConstMatrix<DataType>::ConstMatrix(unsigned int rows, unsigned int columns,
                                       MatrixStorage policy) :
        m_size(),
        m_transpose('N'),
        m_storageType(policy)
    {
        m_size[0] = rows;
        m_size[1] = columns;
    }

    template<typename DataType>
    ConstMatrix<DataType>::ConstMatrix(const ConstMatrix<DataType>& rhs) :
        m_size(),
        m_transpose(rhs.m_transpose),
        m_storageType(rhs.m_storageType)
    {
        m_size[0] = rhs.m_size[0];
        m_size[1] = rhs.m_size[1];
    }

    template<typename DataType>
    ConstMatrix<DataType>& ConstMatrix<DataType>::operator=(const ConstMatrix<DataType>& rhs)
    {
        m_size[0] = rhs.m_size[0];
        m_size[1] = rhs.m_size[1];
        m_transpose = rhs.m_transpose;
        m_storageType = rhs.m_storageType;
        return *this;
    }


    template<typename DataType>
    void ConstMatrix<DataType>::Resize(unsigned int rows, unsigned int columns)
    {
        m_size[0] = rows;
        m_size[1] = columns;
    }

    template<typename DataType>
    void ConstMatrix<DataType>::SetTransposeFlag(char newValue)
    {
        m_transpose = newValue;
    }

    template<typename DataType>
    void ConstMatrix<DataType>::v_Transpose() {}

    template<typename DataType>
    char ConstMatrix<DataType>::v_GetTransposeFlag() const { return m_transpose; }


    template<typename DataType>
    Matrix<DataType>::~Matrix() {}

    template<typename DataType>
    void Matrix<DataType>::SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
    {
        ASSERTL2(row < this->GetRows(), std::string("Row ") + std::to_string(row) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetRows()) +
            std::string(" rows"));
        ASSERTL2(column < this->GetColumns(), std::string("Column ") + std::to_string(column) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetColumns()) +
            std::string(" columns"));
        v_SetValue(row, column, d);
    }

    template<typename DataType>
    Matrix<DataType>::Matrix(unsigned int rows, unsigned int columns,
                             MatrixStorage policy) :
        ConstMatrix<DataType>(rows, columns,policy)
    {
    }

    template<typename DataType>
    Matrix<DataType>::Matrix(const Matrix<DataType>& rhs) :
        ConstMatrix<DataType>(rhs)
    {
    }

    template<typename DataType>
    Matrix<DataType>& Matrix<DataType>::operator=(const Matrix<DataType>& rhs)
    {
        ConstMatrix<DataType>::operator=(rhs);
        return *this;
    }

    template<typename DataType>
    Matrix<DataType>& Matrix<DataType>::operator=(const ConstMatrix<DataType>& rhs)
    {
        ConstMatrix<DataType>::operator=(rhs);
        return *this;
    }


    template LIB_UTILITIES_EXPORT class ConstMatrix<NekDouble>;
    template LIB_UTILITIES_EXPORT class Matrix<NekDouble>;

    template LIB_UTILITIES_EXPORT class ConstMatrix<int>;
    template LIB_UTILITIES_EXPORT class Matrix<int>;

    template LIB_UTILITIES_EXPORT class ConstMatrix<unsigned int>;
    template LIB_UTILITIES_EXPORT class Matrix<unsigned int>;

    template LIB_UTILITIES_EXPORT class ConstMatrix<float>;
    template LIB_UTILITIES_EXPORT class Matrix<float>;
}
