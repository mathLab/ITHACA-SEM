///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrixMetadata.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_METADATA_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_METADATA_HPP

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>


namespace Nektar
{

    // Interface for expression templates.
    class NekMatrixMetadata
    {
        public:
            template<typename MatrixType>
            explicit NekMatrixMetadata(const MatrixType& matrix) :
                Rows(matrix.GetRows()),
                Columns(matrix.GetColumns())
            {
            }

            NekMatrixMetadata(const NekMatrixMetadata& rhs) :
                Rows(rhs.Rows),
                Columns(rhs.Columns)
            {
            }

            static NekMatrixMetadata CreateForNegation(const NekMatrixMetadata& rhs)
            {
                return NekMatrixMetadata(rhs);
            }

            static NekMatrixMetadata CreateForAddition(const NekMatrixMetadata& lhs, const NekMatrixMetadata& rhs)
            {
                ASSERTL1(lhs.Rows == rhs.Rows && lhs.Columns == rhs.Columns, "Matrix dimensions must agree in operator+");
                return NekMatrixMetadata(lhs);
            }

            static NekMatrixMetadata CreateForMultiplication(const NekMatrixMetadata& lhs, const NekMatrixMetadata& rhs)
            {
                ASSERTL1(lhs.Columns == rhs.Rows, "Matrix dimensions must agree in operator*");
                NekMatrixMetadata result;
                result.Rows = lhs.Rows;
                result.Columns = rhs.Columns;
                return result;
            }

            NekMatrixMetadata& operator=(const NekMatrixMetadata& rhs)
            {
                Rows = rhs.Rows;
                Columns = rhs.Columns;
                return *this;
            }

            unsigned int Rows;
            unsigned int Columns;

        private:
            NekMatrixMetadata() :
                Rows(0),
                Columns(0)
            {
            }
    };
    
    class NekBlockMatrixMetadata
    {
        public:
            template<typename MatrixType>
            explicit NekBlockMatrixMetadata(const MatrixType& matrix) :
                Rows(matrix.GetRows()),
                Columns(matrix.GetColumns()),
                BlockRows(matrix.GetBlockRows()),
                BlockColumns(matrix.GetBlockColumns())
            {
            }

            NekBlockMatrixMetadata(const NekBlockMatrixMetadata& rhs) :
                Rows(rhs.Rows),
                Columns(rhs.Columns),
                BlockRows(rhs.BlockRows),
                BlockColumns(rhs.BlockColumns)
            {
            }

            static NekBlockMatrixMetadata CreateForNegation(const NekBlockMatrixMetadata& rhs)
            {
                return NekBlockMatrixMetadata(rhs);
            }

            static NekBlockMatrixMetadata CreateForAddition(const NekBlockMatrixMetadata& lhs, const NekBlockMatrixMetadata& rhs)
            {
                ASSERTL1(lhs.Rows == rhs.Rows && lhs.Columns == rhs.Columns, "Matrix dimensions must agree in operator+");
                ASSERTL1(lhs.BlockRows == rhs.BlockRows && lhs.BlockColumns == rhs.BlockColumns, "Matrix block dimensions must agree in operator+");
                return NekBlockMatrixMetadata(lhs);
            }

//             static NekBlockMatrixMetadata CreateForMultiplication(const NekBlockMatrixMetadata& lhs, const NekBlockMatrixMetadata& rhs)
//             {
//                 ASSERTL1(lhs.Columns == rhs.Rows, "Matrix dimensions must agree in operator*");
//                 NekBlockMatrixMetadata result;
//                 result.Rows = lhs.Rows;
//                 result.Columns = rhs.Columns;
//                 return result;
//             }

            NekBlockMatrixMetadata& operator=(const NekBlockMatrixMetadata& rhs)
            {
                Rows = rhs.Rows;
                Columns = rhs.Columns;
                BlockRows = rhs.BlockRows;
                BlockColumns = rhs.BlockColumns;
                return *this;
            }

            unsigned int Rows;
            unsigned int Columns;
            std::vector<unsigned int> BlockRows;
            std::vector<unsigned int> BlockColumns;

        private:
            NekBlockMatrixMetadata() :
                Rows(0),
                Columns(0),
                BlockRows(),
                BlockColumns()
            {
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_METADATA_HPP

/**
    $Log: NekMatrixMetadata.hpp,v $
    Revision 1.3  2006/09/30 15:18:37  bnelson
    no message

    Revision 1.2  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.1  2006/09/11 03:26:27  bnelson
    Updated to use new policy based expression templates.

 **/


