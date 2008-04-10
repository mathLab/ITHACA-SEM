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

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/Memory/DeleteNothing.hpp>
#include <LibUtilities/LinearAlgebra/MatrixType.h>
#include <LibUtilities/BasicUtils/ConsistentObjectAccess.hpp>
#include <LibUtilities/BasicUtils/RawType.hpp>

#include <boost/shared_ptr.hpp> 
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace Nektar
{
    template<typename DataType>
    struct IsSharedPointer : public boost::false_type {};

    template<typename DataType>
    struct IsSharedPointer<boost::shared_ptr<DataType> > : public boost::true_type {};

    template<typename MatrixType>
    struct LinearSystemSolver;

    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, DiagonalMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, DiagonalMatrixTag, StandardMatrixTag> MatrixType;

        template<typename BVectorType, typename XVectorType>
        static void Solve(const MatrixType& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x)
        {
            ASSERTL0(A.GetColumns() == b.GetRows(), "ERROR: NekLinSys::Solve matrix columns must equal vector rows");

            for(unsigned int i = 0; i < A.GetColumns(); ++i)
            {
                x[i] = b[i]*A(i,i);
            }
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const MatrixType& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x)
        {
            Solve(A,ipivot,b,x);
        }
    };

    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, FullMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, FullMatrixTag, StandardMatrixTag> MatrixType;

        template<typename BVectorType, typename XVectorType>
        static void Solve(const MatrixType& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x)
        {
            x = b;
            int info = 0;
            Lapack::Dgetrs('N',A.GetRows(),1,A.GetPtr().get(),A.GetRows(),(int *)ipivot.get(),x.GetRawPtr(),A.GetRows(),info);
            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dgetrs";
                ASSERTL0(false, message.c_str());
            }
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const MatrixType& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x)
        {
            x = b;
            int info = 0;
            Lapack::Dgetrs('T',A.GetRows(),1,A.GetPtr().get(),A.GetRows(),(int *)ipivot.get(),x.GetRawPtr(),A.GetRows(),info);

            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dgetrs";
                ASSERTL0(false, message.c_str());
            }
        }
    };

    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, SymmetricMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, SymmetricMatrixTag, StandardMatrixTag> MatrixType;

        template<typename BVectorType, typename XVectorType>
        static void Solve(const MatrixType& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x)
        {
            x = b;
            int info = 0;
            Lapack::Dsptrs('U', A.GetRows(), 1, A.GetRawPtr(), ipivot.get(), x.GetRawPtr(), x.GetRows(), info);
            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dsptrs";
                ASSERTL0(false, message.c_str());
            }
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const MatrixType& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x)
        {
            return Solve(A, ipivot, b, x);
        }
    };
    
    template<typename MatrixType>
    class LinearSystem;

    template<typename DataType, typename StorageType, typename Type>
    class LinearSystem<NekMatrix<DataType, StorageType, Type> >
    {
        public:
            typedef NekMatrix<DataType, StorageType, Type> MatrixType;
            typedef LinearSystem<MatrixType> ThisType;

        public:
            explicit LinearSystem(const boost::shared_ptr<MatrixType> &theA) :
                A(*theA) 
            {
                // At some point we should fix this.  We should upate the copy of 
                // A to be transposd for this to work.
                ASSERTL0(theA->GetTransposeFlag() == 'N', "LinearSystem requires a non-transposed matrix.");
                FactorMatrix(A);
            }
            
            explicit LinearSystem(const MatrixType& theA) :
                A(theA)
            {
                // At some point we should fix this.  We should upate the copy of 
                // A to be transposd for this to work.
                ASSERTL0(theA.GetTransposeFlag() == 'N', "LinearSystem requires a non-transposed matrix.");
                FactorMatrix(A);
            }

            LinearSystem(const ThisType& rhs) :
                A(rhs.A),
                m_ipivot(rhs.m_ipivot)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                ThisType temp(rhs);
                swap(temp);
                return *this;
            }

            ~LinearSystem() {}

            // In the following calls to Solve, VectorType must be a NekVector.
            // Anything else won't compile.        
            template<typename VectorType>
            typename RemoveVectorConst<typename RawType<VectorType>::type>::type Solve(const VectorType& b)
            {
                typename RemoveVectorConst<typename RawType<VectorType>::type>::type x(ConsistentObjectAccess<VectorType>::const_reference(b).GetRows());
                LinearSystemSolver<MatrixType>::Solve(A, m_ipivot, ConsistentObjectAccess<VectorType>::const_reference(b), x);
                return x;
            }    

            template<typename BType, typename XType>
            void Solve(const BType& b, XType& x) const
            {
                LinearSystemSolver<MatrixType>::Solve(A, m_ipivot,
                    ConsistentObjectAccess<BType>::const_reference(b), 
                    ConsistentObjectAccess<XType>::reference(x));
            }

            // Transpose variant of solve
            template<typename VectorType>
            typename RemoveVectorConst<typename RawType<VectorType>::type>::type SolveTranspose(const VectorType& b)
            {
                typename RemoveVectorConst<typename RawType<VectorType>::type>::type x(ConsistentObjectAccess<VectorType>::const_reference(b).GetRows());
                LinearSystemSolver<MatrixType>::SolveTranspose(A, m_ipivot, ConsistentObjectAccess<VectorType>::const_reference(b), x);
                return x;
            }    

            template<typename BType, typename XType>
            void SolveTranspose(const BType& b, XType& x) const
            {
                LinearSystemSolver<MatrixType>::SolveTranspose(A, m_ipivot,
                    ConsistentObjectAccess<BType>::const_reference(b), 
                    ConsistentObjectAccess<XType>::reference(x));
            }
            
            

            unsigned int GetRows() const { return A.GetRows(); }
            unsigned int GetColumns() const { return A.GetColumns(); }
            
        private:
            void FactorMatrix(NekMatrix<DataType,DiagonalMatrixTag, StandardMatrixTag> &theA)
            {
                for(unsigned int i = 0; i < A.GetColumns(); ++i)
                {
                    theA.SetValue(i, i, 1.0/theA(i,i));
                }

            }

            void FactorMatrix(NekMatrix<DataType,FullMatrixTag,StandardMatrixTag> &theA)
            {
                int m = theA.GetRows();
                int n = theA.GetColumns();
                int pivotSize = std::max(1, std::min(m, n));
                int info = 0;
                m_ipivot = Array<OneD, int>(pivotSize);

                Lapack::Dgetrf(m, n, theA.GetPtr().get(), m, m_ipivot.get(), info);

                if( info < 0 )
                {
                    std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dgetrf";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +   boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
                    ASSERTL0(false, message.c_str());
                }            
            }

            void FactorMatrix(NekMatrix<DataType, SymmetricMatrixTag, StandardMatrixTag>& theA)
            {
                int info = 0;
                
                Lapack::Dsptrf('U', theA.GetRows(), theA.GetRawPtr(), m_ipivot.get(), info);
                
                if( info < 0 )
                {
                    std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dsptrf";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +   boost::lexical_cast<std::string>(info) + " is 0 from dsptrf";
                    ASSERTL0(false, message.c_str());
                }    
            }
            
            void FactorMatrix(NekMatrix<DataType, BandedMatrixTag, StandardMatrixTag>& theA)
            {
            }
            
            void swap(ThisType& rhs)
            {
                std::swap(A, rhs.A);
                std::swap(m_ipivot,rhs.m_ipivot);
            }

            MatrixType A;
            Array<OneD, int> m_ipivot;        
        };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_HPP
