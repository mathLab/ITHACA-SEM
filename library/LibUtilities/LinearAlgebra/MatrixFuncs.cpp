///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixFuncs.cpp
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
// Description: Matrix functions that depend on storage policy.  Putting 
// methods in these separate classes makes it easier to use them from 
// normal, scaled, and block matrices.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>

namespace Nektar
{

    unsigned int BandedMatrixFuncs::GetRequiredStorageSize(unsigned int totalRows, unsigned int totalColumns,
                                               unsigned int subDiags, unsigned int superDiags)
    {
        return CalculateNumberOfRows(totalRows, subDiags, superDiags)*totalColumns;
    }

    unsigned int BandedMatrixFuncs::CalculateNumberOfDiags(unsigned int totalRows, unsigned int diags)
    {
        if( diags != std::numeric_limits<unsigned int>::max() )
        {
            return diags;
        }
        else if( totalRows > 0 )
        {
            return totalRows-1;
        }
        else
        {
            return 0;
        }
    }
                
    unsigned int BandedMatrixFuncs::CalculateNumberOfRows(unsigned int totalRows, unsigned int subDiags, unsigned int superDiags)
    {
        return CalculateNumberOfDiags(totalRows, subDiags) + CalculateNumberOfDiags(totalRows, superDiags) + 1;
    }

    unsigned int BandedMatrixFuncs::CalculateIndex(unsigned int totalRows, 
                                                        unsigned int totalColumns,
                                                        unsigned int row, unsigned int column,
                                                        unsigned int sub, unsigned int super)
    {
        if( (column <= row && (row - column) <= CalculateNumberOfDiags(totalRows, sub)) ||
            (column > row && (column - row) <= CalculateNumberOfDiags(totalRows, super)) )
        {
            unsigned int elementRow = CalculateNumberOfDiags(totalRows, super)+row-column;
            unsigned int elementColumn = column;

            return elementRow + elementColumn*CalculateNumberOfRows(totalRows, sub, super);
        }
        else
        {
            return std::numeric_limits<unsigned int>::max();
        }
    }


    boost::tuples::tuple<unsigned int, unsigned int> 
    BandedMatrixFuncs::Advance(const unsigned int totalRows, const unsigned int totalColumns,
            const unsigned int curRow, const unsigned int curColumn)
    {
        unsigned int nextRow = curRow;
        unsigned int nextColumn = curColumn;

        //if( transpose == 'N' )
        //{
        //    if( (column <= row && (row - column) <= data.GetNumberOfSubDiagonals(totalRows)) ||
        //    (column > row && (column - row) <= data.GetNumberOfSuperDiagonals(totalRows)) )
        //    {
        //        unsigned int arrayColumns = totalColumns;

        //        unsigned int elementRow = data.GetNumberOfSuperDiagonals(totalRows)+row-column;
        //        unsigned int elementColumn = column;

        //        return elementRow + elementColumn*CalculateNumberOfRows(totalRows, data);
        //    }
        //    else
        //    {
        //        return boost::optional<unsigned int>();
        //    }
        //}
        //else
        //{
        //}
        return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
    }

  
    unsigned int FullMatrixFuncs::GetRequiredStorageSize(unsigned int rows, unsigned int columns)
    {
        return rows*columns;
    }
    
    unsigned int FullMatrixFuncs::CalculateIndex(unsigned int totalRows, unsigned int totalColumns, unsigned int curRow, unsigned int curColumn)
    {
        return curColumn*totalRows + curRow;
    }

    
    boost::tuples::tuple<unsigned int, unsigned int> 
    FullMatrixFuncs::Advance(const unsigned int totalRows, const unsigned int totalColumns,
            const unsigned int curRow, const unsigned int curColumn)
    {
        unsigned int nextRow = curRow;
        unsigned int nextColumn = curColumn;

        if( nextRow < totalRows )
        {
            ++nextRow;
        }

        if( nextRow >= totalRows )
        {
            nextRow = 0;
            ++nextColumn;
        }

        if( nextColumn >= totalColumns )
        {
            nextRow = std::numeric_limits<unsigned int>::max();
            nextColumn = std::numeric_limits<unsigned int>::max();
        }

        return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
    }
        
        
    unsigned int TriangularMatrixFuncs::GetRequiredStorageSize(unsigned int rows, unsigned int columns)
    {
        ASSERTL0(rows==columns, "Triangular matrices must be square.");
        return rows*(rows+1)/2;
    }
    
    unsigned int UpperTriangularMatrixFuncs::CalculateIndex(unsigned int curRow, unsigned int curColumn)
    {
        if( curRow <= curColumn )
        {
            return curRow + curColumn*(curColumn+1)/2;
        }
        else
        {
            return std::numeric_limits<unsigned int>::max();
        }
    }

    boost::tuples::tuple<unsigned int, unsigned int> 
    UpperTriangularMatrixFuncs::Advance(const unsigned int totalRows, const unsigned int totalColumns,
            const unsigned int curRow, const unsigned int curColumn)
    {
        ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
        ASSERTL1(curRow < totalRows, "Attemping to iterate through an element on row " +
            boost::lexical_cast<std::string>(curRow) + " of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
        ASSERTL1(curColumn < totalColumns, "Attemping to iterate through an element on row " +
            boost::lexical_cast<std::string>(curColumn) + " of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");
        ASSERTL1(curRow <= curColumn, "Attemping to iterate through element (" +
            boost::lexical_cast<std::string>(curRow) + ", " +
            boost::lexical_cast<std::string>(curColumn) + ") of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " upper triangular matrix.");

        unsigned int nextRow = curRow;
        unsigned int nextColumn = curColumn;

        if( nextRow <= nextColumn )
        {
            ++nextRow;
        }

        if( nextRow > nextColumn )
        {
            ++nextColumn;
            nextRow = 0;
        }
        
        if( nextColumn >= totalColumns )
        {
            nextRow = std::numeric_limits<unsigned int>::max();
            nextColumn = std::numeric_limits<unsigned int>::max();
        }

        return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
    }

                 
    unsigned int LowerTriangularMatrixFuncs::CalculateIndex(unsigned int totalColumns, unsigned int curRow, unsigned int curColumn)
    {
        if( curRow >= curColumn )
        {
            return curRow + (2*totalColumns - curColumn - 1)*(curColumn)/2;
        }
        else
        {
            return std::numeric_limits<unsigned int>::max();
        }
    }

    boost::tuples::tuple<unsigned int, unsigned int> 
    LowerTriangularMatrixFuncs::Advance(const unsigned int totalRows, const unsigned int totalColumns,
            const unsigned int curRow, const unsigned int curColumn,
            char transpose)
    {
        ASSERTL1(totalRows == totalColumns, "Triangular matrices must be square.");
        ASSERTL1(curRow < totalRows, "Attemping to iterate through an element on row " +
            boost::lexical_cast<std::string>(curRow) + " of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");
        ASSERTL1(curColumn < totalColumns, "Attemping to iterate through an element on row " +
            boost::lexical_cast<std::string>(curColumn) + " of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");
        ASSERTL1(curRow >= curColumn, "Attemping to iterate through element (" +
            boost::lexical_cast<std::string>(curRow) + ", " +
            boost::lexical_cast<std::string>(curColumn) + ") of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " lower triangular matrix.");

        if( transpose == 'T' )
        {
            return UpperTriangularMatrixFuncs::Advance(
                totalColumns, totalRows, curColumn, curRow);
        }

        unsigned int nextRow = curRow;
        unsigned int nextColumn = curColumn;

        if( nextRow < totalRows )
        {
            ++nextRow;
        }

        if( nextRow >= totalRows )
        {
            ++nextColumn;
            nextRow = nextColumn;
        }

        if( nextColumn >= totalColumns )
        {
            nextRow = std::numeric_limits<unsigned int>::max();
            nextColumn = std::numeric_limits<unsigned int>::max();
        }

        return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
    }
    

    unsigned int SymmetricMatrixFuncs::CalculateIndex(unsigned int curRow, unsigned int curColumn)
    {
        if( curRow <= curColumn )
        {
            return curRow + curColumn*(curColumn+1)/2;
        }
        else
        {
            return curColumn + curRow*(curRow + 1)/2;
        }
    }

    boost::tuples::tuple<unsigned int, unsigned int> 
    SymmetricMatrixFuncs::Advance(const unsigned int totalRows, const unsigned int totalColumns,
            const unsigned int curRow, const unsigned int curColumn)
    {
        ASSERTL1(totalRows == totalColumns, "Symmetric matrices must be square.");
        ASSERTL1(curRow < totalRows, "Attemping to iterate through an element on row " +
            boost::lexical_cast<std::string>(curRow) + " of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");
        ASSERTL1(curColumn < totalColumns, "Attemping to iterate through an element on row " +
            boost::lexical_cast<std::string>(curColumn) + " of a (" +
            boost::lexical_cast<std::string>(totalRows) + ", " +
            boost::lexical_cast<std::string>(totalColumns) + " symmetric matrix.");

        unsigned int nextRow = curRow;
        unsigned int nextColumn = curColumn;

        if( nextRow < totalRows )
        {
            ++nextRow;
        }

        if( nextRow >= totalRows )
        {
            nextRow = 0;
            ++nextColumn;
        }

        if( nextColumn >= totalColumns )
        {
            nextRow = std::numeric_limits<unsigned int>::max();
            nextColumn = std::numeric_limits<unsigned int>::max();
        }

        return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
    }

    boost::tuples::tuple<unsigned int, unsigned int> 
    DiagonalMatrixFuncs::Advance(const unsigned int totalRows, const unsigned int totalColumns,
            const unsigned int curRow, const unsigned int curColumn)
    {
        ASSERTL0(curRow == curColumn, "Iteration of a diagonal matrix is only valid along the diagonal.");

        unsigned int nextRow = curRow;
        unsigned int nextColumn = curColumn;

        if( nextRow < totalRows )
        {
            ++nextRow;
            ++nextColumn;
        }

        if( nextRow >= totalRows )
        {
            nextRow = std::numeric_limits<unsigned int>::max();
            nextColumn = std::numeric_limits<unsigned int>::max();
        }

        return boost::tuples::tuple<unsigned int, unsigned int>(nextRow, nextColumn);
    }
    
    unsigned int DiagonalMatrixFuncs::GetRequiredStorageSize(unsigned int rows, unsigned int columns)
    {
        ASSERTL0(rows==columns, "Diagonal matrices must be square.");
        return rows;
    }
    
    unsigned int DiagonalMatrixFuncs::CalculateIndex(unsigned int row, unsigned int col)
    {
        if( row == col )
        {
            return row;
        }
        else
        {
            return std::numeric_limits<unsigned int>::max();
        }
    }

    unsigned int TriangularBandedMatrixFuncs::GetRequiredStorageSize(unsigned int rows, unsigned int columns,
                                                                     unsigned int nSubSuperDiags)
    {
        ASSERTL0(rows==columns, "Triangular matrices must be square.");
        return (nSubSuperDiags+1)*columns;
    }

    unsigned int SymmetricBandedMatrixFuncs::CalculateIndex(unsigned int curRow, unsigned int curColumn, 
                                                            unsigned int nSuperDiags)
    {
        if( curRow <= curColumn )
        {
            if( (curColumn - curRow) <= nSuperDiags )
            {
                unsigned int elementRow = nSuperDiags - (curColumn - curRow);
                unsigned int elementColumn = curColumn;
                
                return elementRow + elementColumn*(nSuperDiags+1);
            }
            else
            {
                return std::numeric_limits<unsigned int>::max();
            }
        }
        else
        {
            return CalculateIndex(curColumn,curRow,nSuperDiags);
        }

    }
}

