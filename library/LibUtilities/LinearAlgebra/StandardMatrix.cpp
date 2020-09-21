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

#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>


namespace Nektar
{


    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix() :
        Matrix<DataType>(0, 0, eFULL),
        m_data(),
        m_wrapperType(eCopy),
        m_numberOfSuperDiagonals(std::numeric_limits<unsigned int>::max()),
        m_numberOfSubDiagonals(std::numeric_limits<unsigned int>::max()),
        m_tempSpace()
    {
        m_data = Array<OneD, DataType>(GetRequiredStorageSize());
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(unsigned int rows, unsigned int columns, MatrixStorage policy,
              unsigned int subDiagonals,
              unsigned int superDiagonals) :
        Matrix<DataType>(rows, columns, policy),
        m_data(),
        m_wrapperType(eCopy),
        m_numberOfSuperDiagonals(superDiagonals),
        m_numberOfSubDiagonals(subDiagonals),
        m_tempSpace()
    {
        m_data = Array<OneD, DataType>(GetRequiredStorageSize());
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(unsigned int rows, unsigned int columns,
              MatrixStorage policy,
              unsigned int subDiagonals,
              unsigned int superDiagonals,
              unsigned int capacity) :
        Matrix<DataType>(rows, columns, policy),
        m_data(),
        m_wrapperType(eCopy),
        m_numberOfSuperDiagonals(superDiagonals),
        m_numberOfSubDiagonals(subDiagonals),
        m_tempSpace()
    {
        unsigned int requiredStorage = this->GetRequiredStorageSize();
        unsigned int actualSize = std::max(requiredStorage, capacity);
        m_data = Array<OneD, DataType>(actualSize);
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue,
              MatrixStorage policy,
              unsigned int subDiagonals,
              unsigned int superDiagonals) :
        Matrix<DataType>(rows, columns, policy),
        m_data(),
        m_wrapperType(eCopy),
        m_numberOfSuperDiagonals(superDiagonals),
        m_numberOfSubDiagonals(subDiagonals),
        m_tempSpace()
    {
        m_data = Array<OneD, DataType>(GetRequiredStorageSize());
        std::fill(m_data.begin(), m_data.end(), initValue);
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(unsigned int rows, unsigned int columns, const DataType* data,
              MatrixStorage policy,
              unsigned int subDiagonals,
              unsigned int superDiagonals) :
        Matrix<DataType>(rows, columns, policy),
        m_data(),
        m_wrapperType(eCopy),
        m_numberOfSuperDiagonals(superDiagonals),
        m_numberOfSubDiagonals(subDiagonals),
        m_tempSpace()
    {
        unsigned int size = GetRequiredStorageSize();
        m_data = Array<OneD, DataType>(size);
        std::copy(data, data + size, m_data.begin());
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d,
              MatrixStorage policy,
              unsigned int subDiagonals,
              unsigned int superDiagonals) :
        Matrix<DataType>(rows, columns, policy),
        m_data(),
        m_wrapperType(eCopy),
        m_numberOfSuperDiagonals(superDiagonals),
        m_numberOfSubDiagonals(subDiagonals),
        m_tempSpace()
    {
        m_data = Array<OneD, DataType>(GetRequiredStorageSize());
        CopyArrayN(d, m_data, m_data.size());
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, DataType>& d, PointerWrapper wrapperType,
              MatrixStorage policy,
              unsigned int subDiagonals,
              unsigned int superDiagonals) :
        Matrix<DataType>(rows, columns, policy),
        m_data(),
        m_wrapperType(wrapperType),
        m_numberOfSuperDiagonals(superDiagonals),
        m_numberOfSubDiagonals(subDiagonals),
        m_tempSpace()
    {
        if( wrapperType == eWrapper )
        {
             m_data = Array<OneD, DataType>(d, eVECTOR_WRAPPER);
        }
        else
        {
            m_data = Array<OneD, DataType>(GetRequiredStorageSize());
            CopyArrayN(d, m_data, m_data.size());
        }
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(const NekMatrix<DataType, StandardMatrixTag>& rhs) :
        Matrix<DataType>(rhs),
        m_data(),
        m_wrapperType(rhs.m_wrapperType),
        m_numberOfSuperDiagonals(rhs.m_numberOfSuperDiagonals),
        m_numberOfSubDiagonals(rhs.m_numberOfSubDiagonals),
        m_tempSpace()
    {
        PerformCopyConstruction(rhs);
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>& NekMatrix<DataType, StandardMatrixTag>::operator=(const NekMatrix<DataType, StandardMatrixTag>& rhs)
    {
        if( this == &rhs )
        {
            return *this;
        }

        Matrix<DataType>::operator=(rhs);
        m_numberOfSubDiagonals = rhs.m_numberOfSubDiagonals;
        m_numberOfSuperDiagonals = rhs.m_numberOfSuperDiagonals;

        ResizeDataArrayIfNeeded();

        unsigned int requiredStorageSize = GetRequiredStorageSize();
        std::copy(rhs.m_data.data(), rhs.m_data.data() + requiredStorageSize, m_data.data());

        return *this;
    }

    /// Fill matrix with scalar
    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>& 
        NekMatrix<DataType, StandardMatrixTag>::operator=(const DataType & rhs)
    {
        unsigned int requiredStorageSize = GetRequiredStorageSize();
        
        DataType* lhs_array = m_data.data();
        
        for (unsigned int i = 0; i < requiredStorageSize; ++i)
        {
            lhs_array[i] = rhs;
        }
        
        return *this;
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::ConstGetValueType NekMatrix<DataType, StandardMatrixTag>::operator()(unsigned int row, unsigned int column) const
    {
        ASSERTL2(row < this->GetRows(), std::string("Row ") + std::to_string(row) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetRows()) +
            std::string(" rows"));
        ASSERTL2(column < this->GetColumns(), std::string("Column ") + std::to_string(column) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetColumns()) +
            std::string(" columns"));

        return this->GetValue(row, column, this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::ConstGetValueType NekMatrix<DataType, StandardMatrixTag>::operator()(unsigned int row, unsigned int column, char transpose) const
    {
        return this->GetValue(row, column, transpose);
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::GetRequiredStorageSize() const
    {
        return Matrix<DataType>::GetRequiredStorageSize(this->GetStorageType(),
            this->GetRows(), this->GetColumns(),
            this->GetNumberOfSubDiagonals(), this->GetNumberOfSuperDiagonals());
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::CalculateIndex(unsigned int row, unsigned int col, const char transpose) const
    {
        unsigned int numRows = this->GetSize()[0];
        unsigned int numColumns = this->GetSize()[1];
        return Matrix<DataType>::CalculateIndex(this->GetStorageType(),
            row, col, numRows, numColumns, transpose,
            m_numberOfSubDiagonals, m_numberOfSuperDiagonals);
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::const_reference NekMatrix<DataType, StandardMatrixTag>::GetValue(unsigned int row, unsigned int column) const
    {
        ASSERTL2(row < this->GetRows(), std::string("Row ") + std::to_string(row) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetRows()) +
            std::string(" rows"));
        ASSERTL2(column < this->GetColumns(), std::string("Column ") + std::to_string(column) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetColumns()) +
            std::string(" columns"));

        return GetValue(row, column, this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::ConstGetValueType NekMatrix<DataType, StandardMatrixTag>::GetValue(unsigned int row, unsigned int column, char transpose) const
    {
        static DataType defaultReturnValue;
        unsigned int index = CalculateIndex(row, column, transpose);
        if( index != std::numeric_limits<unsigned int>::max() )
        {
            return m_data[index];
        }
        else
        {
            return defaultReturnValue;
        }
    }

    template<typename DataType>
    const Array<OneD, const DataType>& NekMatrix<DataType, StandardMatrixTag>::GetPtr() const
    {
        return m_data;
    }

    template<typename DataType>
    DataType NekMatrix<DataType, StandardMatrixTag>::Scale() const
    {
        return DataType(1);
    }

    template<typename DataType>
    const DataType* NekMatrix<DataType, StandardMatrixTag>::GetRawPtr() const
    {
        return m_data.data();
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::const_iterator NekMatrix<DataType, StandardMatrixTag>::begin() const
    {
        return begin(this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::const_iterator NekMatrix<DataType, StandardMatrixTag>::begin(char transpose) const
    {
        if( transpose == 'N' )
        {
            return const_iterator(m_data.data(), m_data.data() + m_data.size());
        }
        else
        {
            return const_iterator(this, transpose);
        }
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::const_iterator NekMatrix<DataType, StandardMatrixTag>::end() const
    {
        return end(this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::const_iterator NekMatrix<DataType, StandardMatrixTag>::end(char transpose) const
    {
        if( transpose == 'N' )
        {
            return const_iterator(m_data.data(), m_data.data() + m_data.size(), true);
        }
        else
        {
            return const_iterator(this, transpose, true);
        }
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::GetStorageSize() const
    {
        return m_data.size();
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::GetNumberOfSubDiagonals() const
    {
        if( m_numberOfSubDiagonals != std::numeric_limits<unsigned int>::max() )
        {
            return m_numberOfSubDiagonals;
        }
        else if( this->GetRows() > 0 )
        {
            return this->GetRows()-1;
        }
        else
        {
            return 0;
        }
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::GetNumberOfSuperDiagonals() const
    {
        if( m_numberOfSuperDiagonals != std::numeric_limits<unsigned int>::max() )
        {
            return m_numberOfSuperDiagonals;
        }
        else if( this->GetRows() > 0 )
        {
            return this->GetRows()-1;
        }
        else
        {
            return 0;
        }
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::CalculateNumberOfRows() const
    {
        return GetNumberOfSubDiagonals() + GetNumberOfSuperDiagonals() + 1;
    }

    template<typename DataType>
    bool NekMatrix<DataType, StandardMatrixTag>::operator==(const NekMatrix<DataType, StandardMatrixTag>& rhs) const
    {
        if( this->GetRows() != rhs.GetRows() ||
            this->GetColumns() != rhs.GetColumns() )
        {
            return false;
        }

        if( this->GetTransposeFlag() == rhs.GetTransposeFlag() )
        {
            return std::equal(begin(), end(), rhs.begin());
        }
        else
        {
            for(unsigned int i = 0; i < this->GetRows(); ++i)
            {
                for(unsigned int j = 0; j < this->GetColumns(); ++j)
                {
                    if( (*this)(i,j) != rhs(i,j) )
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    template<typename DataType>
    PointerWrapper NekMatrix<DataType, StandardMatrixTag>::GetWrapperType() const { return m_wrapperType; }

    template<typename DataType>
    std::tuple<unsigned int, unsigned int>
    NekMatrix<DataType, StandardMatrixTag>::Advance(unsigned int curRow, unsigned int curColumn) const
    {
        return Advance(curRow, curColumn, this->GetTransposeFlag());
    }

    template<typename DataType>
    std::tuple<unsigned int, unsigned int>
    NekMatrix<DataType, StandardMatrixTag>::Advance(unsigned int curRow, unsigned int curColumn, char transpose) const
    {
        unsigned int numRows = this->GetTransposedRows(transpose);
        unsigned int numColumns = this->GetTransposedColumns(transpose);

        switch(this->GetStorageType())
        {
            case eFULL:
                return FullMatrixFuncs::Advance(
                    numRows, numColumns, curRow, curColumn);
                break;
            case eDIAGONAL:
                return DiagonalMatrixFuncs::Advance(
                    numRows, numColumns, curRow, curColumn);
                break;
            case eUPPER_TRIANGULAR:
                return UpperTriangularMatrixFuncs::Advance(
                    numRows, numColumns, curRow, curColumn);
                break;

            case eLOWER_TRIANGULAR:
                return LowerTriangularMatrixFuncs::Advance(
                    numRows, numColumns, curRow, curColumn);
                break;

            case eSYMMETRIC:
            case ePOSITIVE_DEFINITE_SYMMETRIC:
                return SymmetricMatrixFuncs::Advance(
                    numRows, numColumns, curRow, curColumn);
                break;
            case eBANDED:
                return BandedMatrixFuncs::Advance(
                    numRows, numColumns, curRow, curColumn);
                break;
            case eSYMMETRIC_BANDED:
            case ePOSITIVE_DEFINITE_SYMMETRIC_BANDED:
                NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
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
        return std::tuple<unsigned int, unsigned int>(curRow, curColumn);
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag> NekMatrix<DataType, StandardMatrixTag>::CreateWrapper(const NekMatrix<DataType, StandardMatrixTag>& rhs)
    {
        return NekMatrix<DataType, StandardMatrixTag>(rhs, eWrapper);
    }

    template<typename DataType>
    std::shared_ptr<NekMatrix<DataType, StandardMatrixTag> > NekMatrix<DataType, StandardMatrixTag>::CreateWrapper(const std::shared_ptr<NekMatrix<DataType, StandardMatrixTag> >& rhs)
    {
        return std::shared_ptr<NekMatrix<DataType, StandardMatrixTag> >(new NekMatrix<DataType, StandardMatrixTag>(*rhs, eWrapper));
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>::NekMatrix(const NekMatrix<DataType, StandardMatrixTag>& rhs, PointerWrapper wrapperType)  :
        BaseType(rhs),
        m_data(),
        m_wrapperType(wrapperType),
        m_numberOfSuperDiagonals(rhs.m_numberOfSuperDiagonals),
        m_numberOfSubDiagonals(rhs.m_numberOfSubDiagonals)
    {
        PerformCopyConstruction(rhs);
    }


    template<typename DataType>
    Array<OneD, DataType>& NekMatrix<DataType, StandardMatrixTag>::GetData() { return m_data; }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::RemoveExcessCapacity()
    {
        if( m_wrapperType == eCopy )
        {
            unsigned int requiredStorageSize = GetRequiredStorageSize();
            if( m_data.size() > requiredStorageSize )
            {
                Array<OneD, DataType> newArray(requiredStorageSize);
                CopyArrayN(m_data, newArray, requiredStorageSize);
                m_data = newArray;
            }
        }
        else if( m_wrapperType == eWrapper )
        {
            ASSERTL0(true, "Can't call RemoveExcessCapacity on a wrapped matrix.");
        }
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::ResizeDataArrayIfNeeded(unsigned int requiredStorageSize)
    {
        if( m_wrapperType == eCopy  )
        {
            if( m_data.size() < requiredStorageSize )
            {
                Array<OneD, DataType> newData(requiredStorageSize);
                std::copy(m_data.data(), m_data.data() + m_data.size(), newData.data());
                m_data = newData;
            }
        }
        else if( m_wrapperType == eWrapper )
        {
            // If the current matrix is wrapped, then just copy over the top,
            // but the sizes of the two matrices must be the same.
            ASSERTL0(m_data.size() >= requiredStorageSize, "Wrapped NekMatrices must have the same dimension in operator=");
        }
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::ResizeDataArrayIfNeeded()
    {
        unsigned int requiredStorageSize = GetRequiredStorageSize();
        ResizeDataArrayIfNeeded(requiredStorageSize);
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::PerformCopyConstruction(const ThisType& rhs)
    {
        if( m_wrapperType == eWrapper )
        {
            m_data = rhs.m_data;
        }
        else
        {
            m_data = Array<OneD, DataType>(GetRequiredStorageSize());
            CopyArrayN(rhs.m_data, m_data, m_data.size());
        }
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::value_type NekMatrix<DataType, StandardMatrixTag>::v_GetValue(unsigned int row, unsigned int column) const
    {
        return NekMatrix<DataType, StandardMatrixTag>::operator()(row, column);
    }

    template<typename DataType>
    unsigned int NekMatrix<DataType, StandardMatrixTag>::v_GetStorageSize() const
    {
        return NekMatrix<DataType, StandardMatrixTag>::GetStorageSize();
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
    {
        return NekMatrix<DataType, StandardMatrixTag>::SetValue(row, column, d);
    }



    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::SetSize(unsigned int rows, unsigned int cols)
    {
        this->Resize(rows, cols);

        // Some places in Nektar++ access the matrix data array and
        // use size() to see how big it is.  When using
        // expression templates, the data array's capacity is often larger
        // than the actual number of elements, so this statement is
        // required to report the correct number of elements.
        this->GetData().ChangeSize(this->GetRequiredStorageSize());
        ASSERTL0(this->GetRequiredStorageSize() <= this->GetData().size(), "Can't resize matrices if there is not enough capacity.");
    }


    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::Proxy NekMatrix<DataType, StandardMatrixTag>::operator()(unsigned int row, unsigned int column)
    {
        ASSERTL2(row < this->GetRows(), std::string("Row ") + std::to_string(row) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetRows()) +
            std::string(" rows"));
        ASSERTL2(column < this->GetColumns(), std::string("Column ") + std::to_string(column) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetColumns()) +
            std::string(" columns"));

        return (*this)(row, column, this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::Proxy NekMatrix<DataType, StandardMatrixTag>::operator()(unsigned int row, unsigned int column, char transpose)
    {
        unsigned int index = this->CalculateIndex(row, column, transpose);
        if( index != std::numeric_limits<unsigned int>::max() )
        {
            return Proxy(this->GetData()[index]);
        }
        else
        {
            return Proxy();
        }

    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
    {
        ASSERTL2(row < this->GetRows(), std::string("Row ") + std::to_string(row) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetRows()) +
            std::string(" rows"));
        ASSERTL2(column < this->GetColumns(), std::string("Column ") + std::to_string(column) +
            std::string(" requested in a matrix with a maximum of ") + std::to_string(this->GetColumns()) +
            std::string(" columns"));
        SetValue(row, column, d, this->GetTransposeFlag());
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d, char transpose)
    {
        unsigned int index = this->CalculateIndex(row, column, transpose);
        if( index != std::numeric_limits<unsigned int>::max() )
        {
            this->GetData()[index] = d;
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "Can't assign values into zeroed elements of a special array.");
        }
    }

    template<typename DataType>
    Array<OneD, DataType>& NekMatrix<DataType, StandardMatrixTag>::GetPtr()
    {
        return this->GetData();
    }

    template<typename DataType>
    DataType* NekMatrix<DataType, StandardMatrixTag>::GetRawPtr()
    {
        return this->GetData().data();
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::iterator NekMatrix<DataType, StandardMatrixTag>::begin()
    {
        return begin(this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::iterator NekMatrix<DataType, StandardMatrixTag>::begin(char transpose)
    {
        if( transpose == 'N' )
        {
            return iterator(this->GetData().data(), this->GetData().data() + this->GetData().size());
        }
        else
        {
            return iterator(this, transpose);
        }
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::iterator NekMatrix<DataType, StandardMatrixTag>::end()
    {
        return end(this->GetTransposeFlag());
    }

    template<typename DataType>
    typename NekMatrix<DataType, StandardMatrixTag>::iterator NekMatrix<DataType, StandardMatrixTag>::end(char transpose)
    {
        if( transpose == 'N' )
        {
            return iterator(this->GetData().data(), this->GetData().data() + this->GetData().size(), true);
        }
        else
        {
            return iterator(this, transpose, true);
        }
    }

    template<typename DataType>
    NekDouble NekMatrix<DataType, StandardMatrixTag>::AbsMaxtoMinEigenValueRatio(void)
    {
        NekDouble returnval;
        int nvals = this->GetColumns();
        Array<OneD, NekDouble> EigValReal(nvals);
        Array<OneD, NekDouble> EigValImag(nvals);

        EigenSolve(EigValReal,EigValImag);

        Vmath::Vmul(nvals,EigValReal,1,EigValReal,1,EigValReal,1);
        Vmath::Vmul(nvals,EigValImag,1,EigValImag,1,EigValImag,1);
        Vmath::Vadd(nvals,EigValReal,1,EigValImag,1,EigValReal,1);

        returnval = sqrt(Vmath::Vmax(nvals,EigValReal,1)/Vmath::Vmin(nvals,EigValReal,1));

        return returnval;
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::EigenSolve(Array<OneD, NekDouble> &EigValReal,
                    Array<OneD, NekDouble> &EigValImag,
                    Array<OneD, NekDouble> &EigVecs)
    {
        ASSERTL0(this->GetRows()==this->GetColumns(), "Only square matrices can be called");

        switch(this->GetType())
        {
        case eFULL:
            FullMatrixFuncs::EigenSolve(this->GetRows(),
                                        this->GetData(), EigValReal,
                                        EigValImag, EigVecs);
            break;
        case eDIAGONAL:
            Vmath::Vcopy(this->GetRows(),&(this->GetData())[0],1, &EigValReal[0],1);
            Vmath::Zero(this->GetRows(),&EigValImag[0],1);
            break;
        case eUPPER_TRIANGULAR:
        case eLOWER_TRIANGULAR:
        case eSYMMETRIC:
                case eBANDED:
        case eSYMMETRIC_BANDED:
        case eUPPER_TRIANGULAR_BANDED:
        case eLOWER_TRIANGULAR_BANDED:
            NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
            break;

        default:
            NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
        }
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::Invert()
    {
        ASSERTL0(this->GetRows()==this->GetColumns(), "Only square matrices can be inverted.");
        ASSERTL0(this->GetTransposeFlag()=='N', "Only untransposed matrices may be inverted.");

        switch(this->GetType())
        {
            case eFULL:
                FullMatrixFuncs::Invert(this->GetRows(), this->GetColumns(),
                    this->GetData(), this->GetTransposeFlag());
                break;
            case eDIAGONAL:
                DiagonalMatrixFuncs::Invert(this->GetRows(), this->GetColumns(),
                    this->GetData());
                break;
            case eSYMMETRIC:
                SymmetricMatrixFuncs::Invert(this->GetRows(), this->GetColumns(),
                    this->GetData());
                break;
            case eUPPER_TRIANGULAR:
            case eLOWER_TRIANGULAR:
            case eBANDED:
                NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                break;
            case eSYMMETRIC_BANDED:
                NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
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
    }

    template<typename DataType>
    Array<OneD, DataType>& NekMatrix<DataType, StandardMatrixTag>::GetTempSpace()
    {
        if( m_tempSpace.capacity() == 0 )
        {
            m_tempSpace = Array<OneD, DataType>(this->GetData().capacity());
        }
        return m_tempSpace;
    }

    template<typename DataType>
    void NekMatrix<DataType, StandardMatrixTag>::SwapTempAndDataBuffers()
    {
        std::swap(m_tempSpace, this->GetData());
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag> NekMatrix<DataType, StandardMatrixTag>::operator-() const
    {
        NekMatrix<DataType, StandardMatrixTag> result(*this);
        NegateInPlace(result);
        return result;
    }

    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>& NekMatrix<DataType, StandardMatrixTag>::operator*=(const DataType& s)
    {
        for(unsigned int i = 0; i < this->GetPtr().size(); ++i)
        {
            this->GetPtr()[i] *= s;
        }
        return *this;
    }


    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>
    Transpose(NekMatrix<DataType, StandardMatrixTag>& rhs)
    {
        NekMatrix<DataType, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColumns(),
            rhs.GetPtr(), eWrapper, rhs.GetType(), rhs.GetNumberOfSubDiagonals(),
            rhs.GetNumberOfSuperDiagonals());
        result.Transpose();
        return result;
    }

    template LIB_UTILITIES_EXPORT NekMatrix<NekDouble, StandardMatrixTag> Transpose(NekMatrix<NekDouble, StandardMatrixTag>& rhs);

    template LIB_UTILITIES_EXPORT class NekMatrix<NekDouble, StandardMatrixTag>;

    template<typename DataType>
    void NegateInPlace(NekMatrix<DataType, StandardMatrixTag>& m)
    {
        for(unsigned int i = 0; i < m.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < m.GetColumns(); ++j)
            {
                m(i,j) *= -1.0;
            }
        }
    }

    template LIB_UTILITIES_EXPORT void NegateInPlace(NekMatrix<double, StandardMatrixTag>& v);

}


