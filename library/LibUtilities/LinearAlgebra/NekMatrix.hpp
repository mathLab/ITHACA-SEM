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
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/MatrixOperations.hpp>

namespace Nektar
{
    template<typename DataType, typename FormType>
    std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, FormType>& rhs)
    {
        int oswidth = 9;
        int osprecision = 6;

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            os << "[";
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                os.width(oswidth);
                os.precision(osprecision);
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

    template<typename DataType, typename FormType>
    std::ostream& operator>>(std::ostream& os, const NekMatrix<DataType, FormType>& rhs)
    {
        NekDouble tol = 1e-12;
        os <<  "[" << std::endl;

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                if((NekDouble)rhs(i,j) > tol)
                {
                    os << '+';
                }
                else if((NekDouble)rhs(i,j) < -tol)
                {
                    os << '*';
                }
                else
                {
                    os << '-';
                }
            }
            if( i != rhs.GetRows()-1 )
            {
                os << std::endl;
            }
        }
        os <<  "]" << std::endl;
        return os;
    }
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

