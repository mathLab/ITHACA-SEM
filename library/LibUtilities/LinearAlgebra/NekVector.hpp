///////////////////////////////////////////////////////////////////////////////
//
// File: NekVector.hpp
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
// Description: Generic N-Dimensional Vector.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP
#define NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorObject.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixObject.hpp>

#include <functional>
#include <algorithm>
#include <math.h>

#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_array.hpp>


namespace Nektar
{

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    namespace expt
    {
        template<typename DataType, unsigned int dim, unsigned int space>
        class ExpressionTraits<NekVector<DataType, dim, space> >
        {
            public:
                typedef NekVectorMetadata MetadataType;
        };
    }
#endif

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>
    operator+(const NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim, space>& rhs)
    {
        NekVector<DataType, dim, space> result(lhs);
        result += rhs;
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>
    operator-(const NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim, space>& rhs)
    {
        NekVector<DataType, dim, space> result(lhs);
        result -= rhs;
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>
    operator*(const NekVector<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::const_reference rhs)
    {
        NekVector<DataType, dim, space> result(lhs);
        result *= rhs;
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>
    operator*(typename boost::call_traits<DataType>::const_reference lhs, const NekVector<DataType, dim, space>& rhs)
    {
        NekVector<DataType, dim, space> result(rhs);
        result *= lhs;
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>
    operator/(const NekVector<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::const_reference rhs)
    {
        NekVector<DataType, dim, space> result(lhs);
        result /= rhs;
        return result;
    }


    template<typename DataType, unsigned int dim, unsigned int space>
    std::ostream& operator<<(std::ostream& os, const NekVector<DataType, dim, space>& rhs)
    {
        os << rhs.AsString();
        return os;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space> createVectorFromPoints(const NekPoint<DataType, dim, space>& source,
                                                           const NekPoint<DataType, dim, space>& dest)
    {
        NekVector<DataType, dim, space> result;
        for(unsigned int i = 0; i < dim; ++i)
        {
            result[i] = dest[i]-source[i];
        }
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekPoint<DataType, dim, space> findPointAlongVector(const NekVector<DataType, dim, space>& lhs,
                                                        typename boost::call_traits<DataType>::const_reference t)
    {
        NekPoint<DataType, dim, space> result;
        for(unsigned int i = 0; i < dim; ++i)
        {
            result[i] = lhs[i]*t;
        }

        return result;
    }


    template<typename DataType, Nektar::NekMatrixForm lhsForm, MatrixBlockType BlockType, unsigned int space, unsigned int dim>
    NekVector<DataType, dim, space> operator*(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
                                              const NekVector<DataType, dim, space>& rhs)
    {
        NekVector<DataType, dim, space> result(lhs.GetRows(), DataType(0));
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result[i] += lhs(i,j)*rhs(j);
            }
        }

        return result;
    }

//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//    template<typename DataType, Nektar::NekMatrixForm lhsForm, MatrixBlockType BlockType, unsigned int space, unsigned int dim>
//    typename expt::BinaryExpressionType<NekMatrix<DataType, lhsForm, BlockType, space>, 
//                                        expt::MultiplyOp,
//                                        NekVector<DataType, dim, space> >::Type
//    operator*(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
//              const NekVector<DataType, dim, space>& rhs)
//    {
//        return expt::CreateBinaryExpression<expt::MultiplyOp>(lhs, rhs);
//    }
//
//    template<typename DataType, Nektar::NekMatrixForm lhsForm, MatrixBlockType BlockType, unsigned int space, unsigned int dim>
//    class MultiplicationTraits<Nektar::NekMatrix<DataType, lhsForm, BlockType, space>, 
//                               Nektar::NekVector<DataType, dim, space> >
//    {
//        public:
//            typedef Nektar::NekMatrix<DataType, lhsForm, BlockType, space> LhsType;
//            typedef Nektar::NekVector<DataType, dim, space> RhsType;
//            typedef Nektar::NekVectorMetadata MetadataType;
//            typedef Nektar::NekVector<DataType, dim, space> result_type;
//            static const bool HasOpEqual = false;
//            static const bool HasOpLeftEqual = false;
//            
//            static void Multiply(result_type& result, const LhsType& lhs, const RhsType& rhs)
//            {
//                ASSERTL0(result.GetRows() == lhs.GetColumns(), "Dimension error in matrix/vector multiply");
//#ifdef NEKTAR_USING_LAPACK
//                int m = lhs.GetRows();
//                int n = lhs.GetColumns();
//                Blas::Dgemv('T', m, n, 1.0, lhs.GetPtr().get(),  m, rhs.GetPtr().get(),
//                            1. 1.0, result.GetPtr().get(), 1);
//#else
//                for(unsigned int i = 0; i < result.GetRows(); ++i)
//                {
//                    result[i] = DataType(0);
//                    for(unsigned int j = 0; j < result.GetColumns(); ++j)
//                    {
//                        result[i] += lhs(i,j)*rhs(j);
//                    }
//                }
//#endif //NEKTAR_USING_LAPACK
//            }
//
//    };
////#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

}

#endif // NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

/**
    $Log: NekVector.hpp,v $
    Revision 1.16  2007/02/15 06:56:55  bnelson
    *** empty log message ***

    Revision 1.15  2007/02/04 00:15:41  bnelson
    *** empty log message ***

    Revision 1.14  2007/01/29 01:31:08  bnelson
    *** empty log message ***

    Revision 1.13  2007/01/23 03:12:50  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.12  2006/11/29 00:09:32  bnelson
    Removed the LPNorm stub.

    Revision 1.11  2006/11/20 03:38:44  bnelson
    Added /= and /

    Revision 1.10  2006/11/18 17:18:30  bnelson
    Added L1, L2, and Infinity norms

    Revision 1.9  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.8  2006/10/02 01:16:14  bnelson
    Started working on adding BLAS and LAPACK

    Revision 1.7  2006/09/30 15:18:37  bnelson
    no message

    Revision 1.6  2006/09/21 01:00:38  bnelson
    Fixed an expression template problem.

    Revision 1.5  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.4  2006/08/25 01:22:01  bnelson
    no message

    Revision 1.3  2006/08/14 02:29:49  bnelson
    Updated points, vectors, and matrix classes to work with ElVis.  Added a variety of methods to all of these classes.

    Revision 1.2  2006/06/01 13:44:29  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:12:42  kirby
    *** empty log message ***

    Revision 1.2  2006/05/25 03:02:40  bnelson
    Added Matrix/Vector multiplication.

    Revision 1.1  2006/05/04 18:57:44  kirby
    *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


**/

