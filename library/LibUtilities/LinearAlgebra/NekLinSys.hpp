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

#include <LibUtilities/LinearAlgebra/lapack.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

#include <boost/shared_ptr.hpp> 
  
#ifdef NEKTAR_USING_LAPACK
namespace Nektar
{

    template<typename DataType>
    class NekObjectAccessor
    {
        public:
	explicit NekObjectAccessor(DataType& obj) :
                m_object(obj)
            {
            }
            
            NekObjectAccessor(const NekObjectAccessor<DataType>& rhs) :
                m_object(rhs.m_object)
            {
            }
                  
                
            DataType* operator->()
            {
                return &m_object;
            }
           
            
        private:
            NekObjectAccessor<DataType>& operator=(const NekObjectAccessor<DataType>& rhs);
            DataType& m_object;
    };
    
    template<typename DataType>
    class NekObjectAccessor<DataType*>
    {
        public:
            explicit NekObjectAccessor(DataType& obj) :
                m_object(obj)
            {
            }
            
            NekObjectAccessor(const NekObjectAccessor<DataType>& rhs) :
                m_object(rhs.m_object)
            {
            }
                  
                
            DataType* operator->()
            {
                return m_object;
            }
           
            
        private:
            NekObjectAccessor<DataType>& operator=(const NekObjectAccessor<DataType>& rhs);
            DataType* m_object;
    };
            
    // Linear system object.
    // A linear system is one 

    template<typename MatrixType, typename VectorType>
    class LinearSystem
    {
        public:

        private:

    };

    // Ax=b
    template<typename DataType, unsigned int space, unsigned int vectorDim>
    class LinearSystem<NekMatrix<DataType, eDiagonal, eNormal, space>, NekVector<DataType, vectorDim, space> >
    {
        public:
            typedef LinearSystem<NekMatrix<DataType, eDiagonal, eNormal, space>, NekVector<DataType, vectorDim, space> > ThisType;
            typedef NekVector<DataType, vectorDim, space> ResultType; 
        public:
            LinearSystem(const boost::shared_ptr<NekMatrix<DataType, eDiagonal, eNormal, space> >& theA,
                         const boost::shared_ptr<NekVector<DataType, vectorDim, space> >& theB) :
                A(theA),
                b(theB)
            {
            }

            LinearSystem(const ThisType& rhs) :
                A(rhs.A),
                b(rhs.b)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                ThisType temp(rhs);
                swap(temp);
                return *this;
            }

            ~LinearSystem() {}


            ResultType Solve() const
            {
                ASSERTL0(A->GetColumns() == b->GetRows(), "ERROR: NekLinSys::Solve matrix columns must equal vector rows");
                ResultType result(*b);
        
                for(unsigned int i = 0; i < A->GetColumns(); ++i)
                {
                    result[i] = (*b)[i]/(*A)(i,i);
                }
        
                return result;
            }

        private:
            void swap(ThisType& rhs)
            {
                std::swap(A, rhs.A);
                std::swap(b, rhs.b);
            }

            boost::shared_ptr<NekMatrix<DataType, eDiagonal, eNormal, space> > A;
            boost::shared_ptr<NekVector<DataType, vectorDim, space> > b;

    };
    
    template<typename DataType, NekMatrixForm form, unsigned int space, MatrixBlockType BlockType, unsigned int vectorDim>
            class LinearSystem<NekMatrix<DataType, form, BlockType, space>, NekVector<DataType, vectorDim, space> >
    {
        public:
            typedef LinearSystem<NekMatrix<DataType, form, BlockType, space>, NekVector<DataType, vectorDim, space> > ThisType;
            typedef NekVector<DataType, vectorDim, space> ResultType; 
            
        public:
            LinearSystem(const boost::shared_ptr<NekMatrix<DataType, form, BlockType, space> >& theA,
                         const boost::shared_ptr<NekVector<DataType, vectorDim, space> >& theB) :
                A(theA),
                b(theB)
            {
            }

            LinearSystem(const ThisType& rhs) :
                A(rhs.A),
                b(rhs.b)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                ThisType temp(rhs);
                swap(temp);
                return *this;
            }

            ~LinearSystem() {}


            ResultType Solve() const
            {
#ifdef NEKTAR_USING_LAPACK
                ResultType result(*b);
                dgetrs(A->GetRows(), A->GetColumns(), A->GetPtr().get(), result.GetPtr());
                return result;
#else
#error NekLinSys::Solve not yet supported with without lapack support.
#endif //NEKTAR_USING_LAPACK
            }

        private:
            void swap(ThisType& rhs)
            {
                std::swap(A, rhs.A);
                std::swap(b, rhs.b);
            }

            boost::shared_ptr<NekMatrix<DataType, form, BlockType, space> > A;
            boost::shared_ptr<NekVector<DataType, vectorDim, space> > b;

    };
}
#endif //NEKTAR_USING_LAPACK
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_HPP
