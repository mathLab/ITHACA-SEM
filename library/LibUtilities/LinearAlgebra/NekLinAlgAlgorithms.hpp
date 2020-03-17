///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinAlgAlgorithms.hpp
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
// Description: Linear Algebra Algorithms
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LIN_ALG_ALGORITHMS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LIN_ALG_ALGORITHMS_HPP

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <vector>

namespace Nektar
{
    /// \brief Calculates the orthogonal, normalized vectors from
    /// linearly independent input.
    /// \return A list of ortogonalized vectors corresponding to the
    ///         input vectors.  If the size of the output doesn't
    ///         match the size of the input then an error occurred.
    /// This algorithm is taken from "Parallel Scientific Computing in
    /// C++ and MPI", Karniadakis and Kirby, page 55.
    template<typename DataType>
    std::vector<NekVector<DataType> > 
    GramSchmidtOrthogonalization(const std::vector<NekVector<DataType> >& x)
    {
        typedef NekVector<DataType> VectorType;
        typedef NekMatrix<DataType> MatrixType;
        
        //typename dim = x[0].GetDimension(); 
        unsigned int dim = x[0].GetDimension();
        std::vector<VectorType> q(x.size(), VectorType());
        
        // This matrix holds the r_ij values.  Using the matrix object
        // is a convenience since it provides a 2D access to a table
        // of values.
        MatrixType r(dim, dim);
        r(0,0) = x[0].L2Norm();

        if( r(0,0) == DataType(0) )
        {
            return q;
        }

        q[0] = x[0]/r(0,0);

        for(unsigned int j = 1; j < x.size(); ++j)
        {
            for(unsigned int i = 0; i <= j-1; ++i)
            {
                r(i,j) = q[i].Dot(x[j]);
            }
            
            VectorType y = x[j];
            for(unsigned int i = 0; i <= j-1; ++i)
            {
                y = y - r(i,j)*q[i];
            }
            
            r(j,j) = y.L2Norm();
            if( r(j,j) == DataType(0) )
            {
                return q;
            }
            
            q[j] = y/r(j,j);
        }
        
        return q;
    }
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LIN_ALG_ALGORITHMS_HPP

