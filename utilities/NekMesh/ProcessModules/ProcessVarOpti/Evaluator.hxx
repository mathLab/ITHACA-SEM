////////////////////////////////////////////////////////////////////////////////
//
//  File: Evaluator.hxx
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Inline header used to evaluate functional.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_NODEOPTI_EVALUATOR
#define UTILITIES_NEKMESH_NODEOPTI_EVALUATOR

namespace Nektar
{
namespace Utilities
{

using namespace std;

/**
 * @brief Calculate determinant of input matrix.
 *
 * Specialised versions of this function exist only for 1x1, 2x2 and 3x3 matrices.
 *
 * @param jac  Input matrix
 *
 * @return Jacobian of @p jac.
 */
template <int DIM> inline NekDouble Determinant(NekDouble jac[][DIM])
{
    NekDouble det = 0;

    for( unsigned int k=0; k<DIM; ++k )
    {
        int sign = (k%2) ? -1 : 1;

        NekDouble mat[DIM-1][DIM-1];

        for( unsigned int j=1, jj=0; j<DIM; ++j )
        {
            for( unsigned int i=0, ii=0; i<DIM; ++i )
            {
                if( i != k )
                {
                    mat[jj][ii] = jac[j][i];
                    ++ii;
                }
            }

            ++jj;
        }

        det += sign * jac[0][k] * Determinant<DIM-1>( mat );
    }

    return det;
}

template <> inline NekDouble Determinant<1>(NekDouble jac[][1])
{
    return jac[0][0];
}

template <> inline NekDouble Determinant<2>(NekDouble jac[][2])
{
    return jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
}

template <> inline NekDouble Determinant<3>(NekDouble jac[][3])
{
    return jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) -
           jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) +
           jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
}

/**
 * @brief Calculate inverse transpose of input matrix.
 *
 * Specialised versions of this function exist only for 2x2 and 3x3 matrices.
 *
 * @param in   Input matrix \f$ A \f$
 * @param out  Output matrix \f$ A^{-\top} \f$
 */
template <int DIM>
inline void InvTrans(NekDouble in[][DIM], NekDouble out[][DIM])
{
    boost::ignore_unused(in,out);
}

template <> inline void InvTrans<2>(NekDouble in[][2], NekDouble out[][2])
{
    NekDouble invDet = 1.0 / Determinant<2>(in);

    out[0][0] =  in[1][1] * invDet;
    out[1][0] = -in[0][1] * invDet;
    out[0][1] = -in[1][0] * invDet;
    out[1][1] =  in[0][0] * invDet;
}

template <> inline void InvTrans<3>(NekDouble in[][3], NekDouble out[][3])
{
    NekDouble invdet = 1.0 / Determinant<3>(in);

    out[0][0] =  (in[1][1] * in[2][2] - in[2][1] * in[1][2]) * invdet;
    out[1][0] = -(in[0][1] * in[2][2] - in[0][2] * in[2][1]) * invdet;
    out[2][0] =  (in[0][1] * in[1][2] - in[0][2] * in[1][1]) * invdet;
    out[0][1] = -(in[1][0] * in[2][2] - in[1][2] * in[2][0]) * invdet;
    out[1][1] =  (in[0][0] * in[2][2] - in[0][2] * in[2][0]) * invdet;
    out[2][1] = -(in[0][0] * in[1][2] - in[1][0] * in[0][2]) * invdet;
    out[0][2] =  (in[1][0] * in[2][1] - in[2][0] * in[1][1]) * invdet;
    out[1][2] = -(in[0][0] * in[2][1] - in[2][0] * in[0][1]) * invdet;
    out[2][2] =  (in[0][0] * in[1][1] - in[1][0] * in[0][1]) * invdet;
}

/**
 * @brief Calculate Scalar product of input vectors.
 */
/*template<int DIM>
inline NekDouble ScalarProd(NekDouble in1[DIM], NekDouble in2[DIM])
{
    return 0.0;
}

template <> inline NekDouble ScalarProd<2>(NekDouble in1[2], NekDouble in2[2])
{
    return    in1[0] * in2[0]
            + in1[1] * in2[1];
}

template <> inline NekDouble ScalarProd<3>(NekDouble in1[3], NekDouble in2[3])
{
    return    in1[0] * in2[0]
            + in1[1] * in2[1]
            + in1[2] * in2[2];
}*/
template<int DIM>
inline NekDouble ScalarProd(NekDouble (&in1)[DIM], NekDouble (&in2)[DIM])
{
    boost::ignore_unused(in1,in2);
    return 0.0;
}

template<>
inline NekDouble ScalarProd<2>(NekDouble (&in1)[2], NekDouble (&in2)[2])
{
    return    in1[0] * in2[0]
            + in1[1] * in2[1];
}
template<>
inline NekDouble ScalarProd<3>(NekDouble (&in1)[3], NekDouble (&in2)[3])
{
    return    in1[0] * in2[0]
            + in1[1] * in2[1]
            + in1[2] * in2[2];
}


/**
 * @brief Calculate \f$ E = F^\top F - I \f$ tensor used in derivation of linear
 * elasticity gradients.
 *
 * Specialised versions of this function exist only for 2x2 and 3x3 matrices.
 *
 * @param in   Input matrix \f$ F \f$
 * @param out  Output matrix \f$ F^\top F - I \f$
 */
template <int DIM>
inline void EMatrix(NekDouble in[][DIM], NekDouble out[][DIM])
{
    boost::ignore_unused(in,out);
}

template <> inline void EMatrix<2>(NekDouble in[][2], NekDouble out[][2])
{
    out[0][0] = 0.5 * (in[0][0] * in[0][0] + in[1][0] * in[1][0] - 1.0);
    out[1][0] = 0.5 * (in[0][0] * in[0][1] + in[1][0] * in[1][1]);
    out[0][1] = 0.5 * (in[0][0] * in[0][1] + in[1][0] * in[1][1]);
    out[1][1] = 0.5 * (in[1][1] * in[1][1] + in[0][1] * in[0][1] - 1.0);
}

template <> inline void EMatrix<3>(NekDouble in[][3], NekDouble out[][3])
{
    out[0][0] = 0.5 * (in[0][0] * in[0][0] + in[1][0] * in[1][0] +
                       in[2][0] * in[2][0] - 1.0);
    out[1][0] = 0.5 * (in[0][0] * in[1][0] + in[1][0] * in[1][1] +
                       in[2][0] * in[2][1]);
    out[0][1] = out[1][0];
    out[2][0] = 0.5 * (in[0][0] * in[0][2] + in[1][0] * in[1][2] +
                       in[2][0] * in[2][2]);
    out[0][2] = out[2][0];
    out[1][1] = 0.5 * (in[0][1] * in[0][1] + in[1][1] * in[1][1] +
                       in[2][1] * in[2][1] - 1.0);
    out[1][2] = 0.5 * (in[0][1] * in[0][2] + in[1][1] * in[1][2] +
                       in[2][1] * in[2][2]);
    out[2][1] = out[1][2];
    out[2][2] = 0.5 * (in[0][2] * in[0][2] + in[1][2] * in[1][2] +
                       in[2][2] * in[2][2] - 1.0);
}



/**
 * @brief Calculate Frobenius inner product of input matrices.
 */
template<int DIM>
inline NekDouble FrobProd(NekDouble in1[][DIM],
                          NekDouble in2[][DIM])
{
    boost::ignore_unused(in1,in2);
    return 0.0;
}

template<>
inline NekDouble FrobProd<2>(NekDouble in1[][2], NekDouble in2[][2])
{
    return    in1[0][0] * in2[0][0]
            + in1[0][1] * in2[0][1]
            + in1[1][0] * in2[1][0]
            + in1[1][1] * in2[1][1] ;
}

template<>
inline NekDouble FrobProd<3>(NekDouble in1[][3], NekDouble in2[][3])
{
    return    in1[0][0] * in2[0][0]
            + in1[0][1] * in2[0][1]
            + in1[0][2] * in2[0][2]
            + in1[1][0] * in2[1][0]
            + in1[1][1] * in2[1][1]
            + in1[1][2] * in2[1][2]
            + in1[2][0] * in2[2][0]
            + in1[2][1] * in2[2][1]
            + in1[2][2] * in2[2][2] ;
}

/**
 * @brief Calculate Frobenius norm \f$ \| A \|_f ^2 \f$ of a matrix \f$ A \f$.
 *
 * @param inarray   Input matrix \f$ A \f$
 */
template<int DIM>
inline NekDouble FrobeniusNorm(NekDouble inarray[][DIM])
{
    boost::ignore_unused(inarray);
    return 0.0;
}

template<>
inline NekDouble FrobeniusNorm<2>(NekDouble inarray[][2])
{
    return    inarray[0][0] * inarray[0][0]
            + inarray[0][1] * inarray[0][1]
            + inarray[1][0] * inarray[1][0]
            + inarray[1][1] * inarray[1][1] ;
}

template<>
inline NekDouble FrobeniusNorm<3>(NekDouble inarray[][3])
{
    return    inarray[0][0] * inarray[0][0]
            + inarray[0][1] * inarray[0][1]
            + inarray[0][2] * inarray[0][2]
            + inarray[1][0] * inarray[1][0]
            + inarray[1][1] * inarray[1][1]
            + inarray[1][2] * inarray[1][2]
            + inarray[2][0] * inarray[2][0]
            + inarray[2][1] * inarray[2][1]
            + inarray[2][2] * inarray[2][2] ;
}

}
}

#endif
