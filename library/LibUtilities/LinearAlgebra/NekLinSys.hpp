///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinSys.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_HPP

#include <LibUtilities/BasicUtils/Lapack.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/Memory/DeleteNothing.hpp>

#include <boost/shared_ptr.hpp> 
  
namespace Nektar
{
    template<typename MatrixType, typename VectorType>
    struct LinearSystemSolver;

    template<typename DataType, unsigned int space, unsigned int vectorDim>
    struct LinearSystemSolver<NekMatrix<DataType, eDiagonal, eNormal, space>, NekVector<DataType, vectorDim, space> >
    {
        typedef NekVector<DataType, vectorDim, space> VectorType;
        typedef NekMatrix<DataType, eDiagonal, eNormal, space> MatrixType;

        static void Solve(const boost::shared_ptr<MatrixType>& A, const VectorType& b, VectorType& x)
        {
            ASSERTL0(A->GetColumns() == b.GetRows(), "ERROR: NekLinSys::Solve matrix columns must equal vector rows");

            for(unsigned int i = 0; i < A->GetColumns(); ++i)
            {
                x[i] = b[i]/(*A)(i,i);
            }
        }
    };

    template<typename DataType, NekMatrixForm form, unsigned int space, MatrixBlockType BlockType, unsigned int vectorDim>
    struct LinearSystemSolver<NekMatrix<DataType, form, BlockType, space>, NekVector<DataType, vectorDim, space> >
    {
        typedef NekMatrix<DataType, form, BlockType, space> MatrixType;
        typedef NekVector<DataType, vectorDim, space> VectorType;

        static void Solve(const boost::shared_ptr<MatrixType>& A, const VectorType& b, VectorType& x)
        {
            Lapack::dgetrs(A->GetRows(),A->GetColumns(),A->GetPtr().get(),x.GetPtr());
        }
    };

    template<typename MatrixType>
    class LinearSystem;

    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    class LinearSystem<NekMatrix<DataType, form, BlockType, space> >
    {
        public:
            typedef NekMatrix<DataType, form, BlockType, space> MatrixType;
            typedef LinearSystem<MatrixType> ThisType;

        public:
            explicit LinearSystem(const boost::shared_ptr<MatrixType>& theA) :
                A(theA)
            {
            }

            LinearSystem(const ThisType& rhs) :
                A(rhs.A)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                ThisType temp(rhs);
                swap(temp);
                return *this;
            }

            ~LinearSystem() {}


            template<unsigned int dim>
            NekVector<DataType, dim, space> Solve(boost::shared_ptr<NekVector<DataType, dim, space> > b) const
            {
                typedef NekVector<DataType, dim, space> VectorType;
                VectorType x(*b);
                LinearSystemSolver<MatrixType, VectorType>::Solve(A, *b, x);
                return x;
            }

            template<unsigned int dim>
            NekVector<DataType, dim, space> Solve(const NekVector<DataType, dim, space>& b) const
            {
                typedef NekVector<DataType, dim, space> VectorType;
                VectorType x(b);
                LinearSystemSolver<MatrixType, VectorType>::Solve(A, b, x);
                return x;
            }

            //template<typename VectorType>
            //VectorType Solve(const boost::shared_ptr<const VectorType>& b,
            //                 const boost::shared_ptr<VectorType>& x) const
            //{
            //    return LinearSystemSolver<MatrixType, VectorType>::Solve(A, *b, *x);
            //}

            //template<typename VectorType>
            //VectorType Solve(const boost::shared_ptr<const VectorType>& b,
            //                 const VectorType& x) const
            //{
            //    return LinearSystemSolver<MatrixType, VectorType>::Solve(A, *b, x);
            //}

            //template<typename VectorType>
            //VectorType Solve(const VectorType& b,
            //                 const boost::shared_ptr<VectorType>& x) const
            //{
            //    return LinearSystemSolver<MatrixType, VectorType>::Solve(A, b, *x);
            //}

            //template<typename VectorType>
            //VectorType Solve(const VectorType& b,
            //                 const VectorType& x) const
            //{
            //    return LinearSystemSolver<MatrixType, VectorType>::Solve(A, b, *x);
            //}

        private:
            void swap(ThisType& rhs)
            {
                std::swap(A, rhs.A);
            }

            boost::shared_ptr<MatrixType> A;
    };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_HPP
