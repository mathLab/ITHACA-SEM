///////////////////////////////////////////////////////////////////////////////
//
// File GlobalMatrix.cpp
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
// Description: GlobalMatrix definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/SparseBlas.hpp>

#include <MultiRegions/GlobalMatrix.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalMatrix
         * This matrix is essentially a wrapper around a DNekSparseMat.
         */

        /**
         * Allocates a new DNekSparseMat object from the given specification.
         * @param   rows        Number of rows in matrix.
         * @param   columns     Number of columns in matrix.
         * @param   cooMat      ?
         */
        GlobalMatrix::GlobalMatrix(unsigned int rows, 
                                   unsigned int columns,
                                   const COOMatType &cooMat):
            m_matrix(MemoryManager<DNekSparseMat>
                                        ::AllocateSharedPtr(rows,columns,
                                                            cooMat))
        {
        }


        /**
         * Performs a matrix-vector multiply using the Sparse BLAS routine
         * DCSRMV.
         * @param   in          Input vector.
         * @param   out         Output vector.
         */
        void GlobalMatrix::Multiply(const Array<OneD,const NekDouble> &in, 
                                          Array<OneD,      NekDouble> &out)
        {
            SparseBlas::Dcsrmv(*m_matrix,in,out);
        }
	
    } //end of namespace
} //end of namespace

