///////////////////////////////////////////////////////////////////////////////
//
// File: lapack.cpp
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

#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <algorithm>
#include <string>
#include <boost/lexical_cast.hpp>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace Lapack
{
    void dgetrs(char trans, int matrixRows, int matrixColumns, const double* A, double* x)
    {
        // The factorization step replaces the input matrix, which we can't allow 
        // (because A is a NekMatrix object and the user will not expect it to be changed).
        // So we need to copy.
        
        Nektar::NekMatrix<double, Nektar::eFull> factoredMatrix(matrixRows, matrixColumns, A);
        factoredMatrix.Transpose();
        
        // Step 1 - Obtain the LU factorization of the matrix.
        // dgetrf expects a general m by n matrix.
        int m = matrixRows;
        int n = matrixColumns;
        
        // Pivot information
        int pivotSize = std::max(1, std::min(m, n));
        int info = 0;
        Nektar::Int1DSharedArray ipivot = Nektar::MemoryManager::AllocateSharedPtr<Nektar::Int1DSharedArrayBase>(boost::extents[pivotSize]);
        Lapack::Dgetrf(m, n, factoredMatrix.GetPtr().get(), m, ipivot->data(), info);

        if( info < 0 )
        {
            std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
                    "th parameter had an illegal parameter for dgetrf";
            ASSERTL0(false, message.c_str());
        }
        else if( info > 0 )
        {
            std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                    boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
            ASSERTL0(false, message.c_str());
        }
        // Now that we have a factored matrix, solve it.
        //void    dgetrs(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);
        int nrhs = 1; // ONly 1 right hand side.

        Lapack::Dgetrs(trans, matrixRows, 1, factoredMatrix.GetPtr().get(), matrixRows, ipivot->data(), x, matrixRows, info);
        
        if( info < 0 )
        {
            std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
                    "th parameter had an illegal parameter for dgetrs";
            ASSERTL0(false, message.c_str());
        }
    }
}

