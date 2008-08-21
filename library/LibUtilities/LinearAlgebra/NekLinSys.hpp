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
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <iostream>

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
        typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;
        
        template<typename BVectorType, typename XVectorType>
        static void Solve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            for(unsigned int i = 0; i < A.num_elements(); ++i)
            {
                x[i] = b[i]*A[i];
            }
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            Solve(A,ipivot,b,x,n,policyData);
        }
    };

    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, FullMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, FullMatrixTag, StandardMatrixTag> MatrixType;
        typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;
        
        template<typename BVectorType, typename XVectorType>
        static void Solve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            x = b;
            int info = 0;
            Lapack::Dgetrs('N',n,1,A.get(),n,(int *)ipivot.get(),x.GetRawPtr(),n,info);
            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dgetrs";
                ASSERTL0(false, message.c_str());
            }
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const Array<OneD, double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            x = b;
            int info = 0;
            Lapack::Dgetrs('T',n,1,A.get(),n,(int *)ipivot.get(),x.GetRawPtr(), n,info);

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
        typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;
        
        template<typename BVectorType, typename XVectorType>
        static void Solve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            x = b;
            int info = 0;
            Lapack::Dsptrs('U', n, 1, A.get(), ipivot.get(), x.GetRawPtr(), x.GetRows(), info);
            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dsptrs";
                ASSERTL0(false, message.c_str());
            }
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            return Solve(A, ipivot, b, x, n, policyData);
        }
    };
    
    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, BandedMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, BandedMatrixTag, StandardMatrixTag> MatrixType;
        typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;

        template<typename BVectorType, typename XVectorType>
        static void Solve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            PerformSolve(A, ipivot, b, x, 'N', n, policyData);
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            PerformSolve(A, ipivot, b, x, 'T', n, policyData);
        }
        
        template<typename BVectorType, typename XVectorType>
        static void PerformSolve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, char trans, unsigned int n, const PolicySpecificDataType& policyData)
        {
            x = b;
            int KL = policyData.GetNumberOfSubDiagonals(n);
            int KU = policyData.GetNumberOfSuperDiagonals(n);
            int info = 0;
            
            Lapack::Dgbtrs(trans, n, KL, KU, 1, A.get(), 2*KL+KU+1, ipivot.get(), x.GetRawPtr(), n, info);
            
            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dgbtrs";
                ASSERTL0(false, message.c_str());
            }
        }
    };

    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag> MatrixType;
        typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;

        template<typename BVectorType, typename XVectorType>
        static void Solve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            PerformSolve(A, ipivot, b, x, 'N', n, policyData);
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            PerformSolve(A, ipivot, b, x, 'T', n, policyData);
        }

        template<typename BVectorType, typename XVectorType>
        static void PerformSolve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, char trans, unsigned int n, const PolicySpecificDataType& policyData)
        {
            x = b;
            int info = 0;
            Lapack::Dtptrs('L', trans, 'N', n, 1, A.get(), x.GetRawPtr(), n, info);

            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dtrtrs";
                ASSERTL0(false, message.c_str());
            }
            else if( info > 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th diagonal element of A is 0 for dtrtrs";
                ASSERTL0(false, message.c_str());
            }
        }
    };

    template<typename DataType>
    struct LinearSystemSolver<NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag> >
    {
        typedef NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag> MatrixType;
        typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;

        template<typename BVectorType, typename XVectorType>
        static void Solve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            PerformSolve(A, ipivot, b, x, 'N', n, policyData);
        }

        template<typename BVectorType, typename XVectorType>
        static void SolveTranspose(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, unsigned int n, const PolicySpecificDataType& policyData)
        {
            PerformSolve(A, ipivot, b, x, 'T', n, policyData);
        }

        template<typename BVectorType, typename XVectorType>
        static void PerformSolve(const Array<OneD, const double>& A, const Array<OneD, const int>& ipivot, const BVectorType& b, XVectorType& x, char trans, unsigned int n, const PolicySpecificDataType& policyData)
        {
            x = b;
            int info = 0;
            Lapack::Dtptrs('U', trans, 'N', n, 1, A.get(), x.GetRawPtr(), n, info);

            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dtrtrs";
                ASSERTL0(false, message.c_str());
            }
            else if( info > 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th diagonal element of A is 0 for dtrtrs";
                ASSERTL0(false, message.c_str());
            }
        }
    };

    template<typename StorageType, typename Type>
    class LinearSystem<NekMatrix<double, StorageType, Type> >
    {
        public:
            typedef NekMatrix<double, StorageType, Type> MatrixType;
            typedef LinearSystem<MatrixType> ThisType;
            typedef typename MatrixType::PolicySpecificDataHolderType PolicySpecificDataType;
            
        public:
            explicit LinearSystem(const boost::shared_ptr<MatrixType> &theA, PointerWrapper wrapperType = eCopy) :
                n(theA->GetRows()),
                A(theA->GetPtr(), eVECTOR_WRAPPER),
                m_ipivot(),
                m_policySpecificData(theA->GetPolicySpecificDataHolderType())
            {
                // At some point we should fix this.  We should upate the copy of 
                // A to be transposd for this to work.
                ASSERTL0(theA->GetTransposeFlag() == 'N', "LinearSystem requires a non-transposed matrix.");
                ASSERTL0( (wrapperType == eWrapper && !boost::is_same<StorageType, BandedMatrixTag>::value) || wrapperType == eCopy , "Banded matrices can't be wrapped");
                
                if( wrapperType == eCopy )
                {
                    A = Array<OneD, double>(theA->GetPtr().num_elements());
                    CopyArray(theA->GetPtr(), A);
                }
                
                FactorMatrix(*theA);
            }
            
            explicit LinearSystem(const MatrixType& theA, PointerWrapper wrapperType = eCopy) :
                n(theA.GetRows()),
                A(theA.GetPtr(), eVECTOR_WRAPPER),
                m_ipivot(),
                m_policySpecificData(theA.GetPolicySpecificDataHolderType())
            {
                // At some point we should fix this.  We should upate the copy of 
                // A to be transposd for this to work.
                ASSERTL0(theA.GetTransposeFlag() == 'N', "LinearSystem requires a non-transposed matrix.");
                ASSERTL0( (wrapperType == eWrapper && !boost::is_same<StorageType, BandedMatrixTag>::value) || wrapperType == eCopy, "Banded matrices can't be wrapped" );
                
                if( wrapperType == eCopy )
                {
                    A = Array<OneD, double>(theA.GetPtr().num_elements());
                    CopyArray(theA.GetPtr(), A);
                }

                FactorMatrix(theA);
            }

            LinearSystem(const ThisType& rhs) :
                n(rhs.n),
                A(rhs.A),
                m_ipivot(rhs.m_ipivot),
                m_policySpecificData(rhs.m_policySpecificData)
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
                LinearSystemSolver<MatrixType>::Solve(A, m_ipivot, ConsistentObjectAccess<VectorType>::const_reference(b), x, n, m_policySpecificData);
                return x;
            }    

            template<typename BType, typename XType>
            void Solve(const BType& b, XType& x) const
            {
                LinearSystemSolver<MatrixType>::Solve(A, m_ipivot,
                    ConsistentObjectAccess<BType>::const_reference(b), 
                    ConsistentObjectAccess<XType>::reference(x), n, m_policySpecificData);
            }

            // Transpose variant of solve
            template<typename VectorType>
            typename RemoveVectorConst<typename RawType<VectorType>::type>::type SolveTranspose(const VectorType& b)
            {
                typename RemoveVectorConst<typename RawType<VectorType>::type>::type x(ConsistentObjectAccess<VectorType>::const_reference(b).GetRows());
                LinearSystemSolver<MatrixType>::SolveTranspose(A, m_ipivot, ConsistentObjectAccess<VectorType>::const_reference(b), x, n, m_policySpecificData);
                return x;
            }    

            template<typename BType, typename XType>
            void SolveTranspose(const BType& b, XType& x) const
            {
                LinearSystemSolver<MatrixType>::SolveTranspose(A, m_ipivot,
                    ConsistentObjectAccess<BType>::const_reference(b), 
                    ConsistentObjectAccess<XType>::reference(x), n, m_policySpecificData);
            }
            
            

            unsigned int GetRows() const { return n; }
            unsigned int GetColumns() const { return n; }
            
        private:
            void FactorMatrix(const NekMatrix<double, DiagonalMatrixTag, StandardMatrixTag> &theA)
            {
                for(unsigned int i = 0; i < theA.GetColumns(); ++i)
                {
                    A[i] = 1.0/theA(i,i);
                }

            }

            template<typename TriangularType>
            void FactorMatrix(const NekMatrix<double, TriangularType, StandardMatrixTag>& theA,
                              typename boost::enable_if
                                <
                                    boost::mpl::or_
                                    <
                                        boost::is_same<TriangularType, UpperTriangularMatrixTag>,
                                        boost::is_same<TriangularType, LowerTriangularMatrixTag>
                                    >
                                >::type* p = 0)
            {
            }

            void FactorMatrix(const NekMatrix<double,FullMatrixTag,StandardMatrixTag> &theA)
            {
                int m = theA.GetRows();
                int n = theA.GetColumns();
                
                int pivotSize = std::max(1, std::min(m, n));
                int info = 0;
                m_ipivot = Array<OneD, int>(pivotSize);

                Lapack::Dgetrf(m, n, A.get(), m, m_ipivot.get(), info);

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

            void FactorMatrix(const NekMatrix<double, SymmetricMatrixTag, StandardMatrixTag>& theA)
            {
                int info = 0;
                int pivotSize = theA.GetRows();
                m_ipivot = Array<OneD, int>(pivotSize);
                
                Lapack::Dsptrf('U', theA.GetRows(), A.get(), m_ipivot.get(), info);
                
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

            void FactorMatrix(const NekMatrix<double, BandedMatrixTag, StandardMatrixTag>& theA)
            {
                int M = theA.GetRows();
                int N = theA.GetColumns();
                int KL = m_policySpecificData.GetNumberOfSubDiagonals(M);
                int KU = m_policySpecificData.GetNumberOfSuperDiagonals(M);
                
                // The array we pass in to dgbtrf must have enough space for KL
                // subdiagonals and KL+KU superdiagonals (see lapack users guide,
                // in the section discussing band storage.
                PolicySpecificDataType t(KL, KL+KU);
                unsigned int requiredStorageSize = MatrixStoragePolicy<double, BandedMatrixTag>::GetRequiredStorageSize(
                    M, N, t);
                
                unsigned int rawRows = KL+KU+1;
                A = Array<OneD, double>(requiredStorageSize);

                // Put the extra elements up front.
                for(unsigned int i = 0; i < theA.GetColumns(); ++i)
                {
                    std::copy(theA.GetRawPtr() + i*rawRows, theA.GetRawPtr() + (i+1)*rawRows,
                        A.get() + (i+1)*KL + i*rawRows);
                }
                       
                int info = 0;
                int pivotSize = theA.GetRows();
                m_ipivot = Array<OneD, int>(pivotSize);
                
                Lapack::Dgbtrf(M, N, KL, KU, A.get(), 2*KL+KU+1, m_ipivot.get(), info);
              
                if( info < 0 )
                {
                    std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + "th parameter had an illegal parameter for dgbtrf";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +   boost::lexical_cast<std::string>(info) + " is 0 from dgbtrf";
                    ASSERTL0(false, message.c_str());
                }    
            }

            
            void swap(ThisType& rhs)
            {
                std::swap(n, rhs.n);
                std::swap(A, rhs.A);
                std::swap(m_ipivot,rhs.m_ipivot);
                std::swap(m_policySpecificData, rhs.m_policySpecificData);
            }

            unsigned int n;
            Array<OneD, double> A;
            Array<OneD, int> m_ipivot;  
            PolicySpecificDataType m_policySpecificData;      
        };
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_HPP
