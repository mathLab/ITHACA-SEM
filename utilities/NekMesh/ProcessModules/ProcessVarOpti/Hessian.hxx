////////////////////////////////////////////////////////////////////////////////
//
//  File: Hessian.hxx
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
//  Description: Utility functions for Hessian matrices
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_NODEOPTI_HESSIAN
#define UTILITIES_NEKMESH_NODEOPTI_HESSIAN


namespace Nektar
{
namespace Utilities
{

/**
 * @brief Returns 1 if Hessian matrix is indefinite and 0 otherwise.
 *
 * Specialised versions of this function exist only for 2x2 and 3x3 matrices.
 */
template <int DIM> int NodeOpti::IsIndefinite()
{
    ASSERTL0(false, "DIM error");
    return 0;
}

template <> int NodeOpti::IsIndefinite<2>()
{
    vector<NekDouble> eigR(2);
    vector<NekDouble> eigI(2);
    NekMatrix<NekDouble> H(2, 2);
    H(0, 0) = m_grad[2];
    H(1, 0) = m_grad[3];
    H(0, 1) = H(1, 0);
    H(1, 1) = m_grad[4];

    // cout << H << endl << endl;

    int nVel   = 2;
    char jobvl = 'N', jobvr = 'N';
    int worklen = 8 * nVel, info;
    NekDouble dum;

    DNekMat eval(nVel, nVel, 0.0, eDIAGONAL);
    vector<NekDouble> vl(nVel * nVel);
    vector<NekDouble> work(worklen);
    vector<NekDouble> wi(nVel);

    Lapack::Dgeev(jobvl, jobvr, nVel, H.GetRawPtr(), nVel, &(eval.GetPtr())[0],
                  &wi[0], &vl[0], nVel, &dum, nVel, &work[0], worklen, info);

    ASSERTL0(!info, "dgeev failed");

    if (eval(0, 0) < 0.0 || eval(1, 1) < 0.0)
    {
        if (eval(0, 0) < 0.0 && eval(1, 1) < 0.0)
        {
            return 2;
        }
        else
        {
            return 1;
        }
    }

    return 0;
}

template <> int NodeOpti::IsIndefinite<3>()
{
    vector<NekDouble> eigR(3);
    vector<NekDouble> eigI(3);
    NekMatrix<NekDouble> H(3, 3);
    H(0, 0) = m_grad[3];
    H(1, 0) = m_grad[4];
    H(0, 1) = H(1, 0);
    H(2, 0) = m_grad[5];
    H(0, 2) = H(2, 0);
    H(1, 1) = m_grad[6];
    H(2, 1) = m_grad[7];
    H(1, 2) = H(2, 1);
    H(2, 2) = m_grad[8];

    int nVel   = 3;
    char jobvl = 'N', jobvr = 'N';
    int worklen = 8 * nVel, info;
    NekDouble dum;

    DNekMat eval(nVel, nVel, 0.0, eDIAGONAL);
    vector<NekDouble> vl(nVel * nVel);
    vector<NekDouble> work(worklen);
    vector<NekDouble> wi(nVel);

    Lapack::Dgeev(jobvl, jobvr, nVel, H.GetRawPtr(), nVel, &(eval.GetPtr())[0],
                  &wi[0], &vl[0], nVel, &dum, nVel, &work[0], worklen, info);

    ASSERTL0(!info, "dgeev failed");

    if (eval(0, 0) < 0.0 || eval(1, 1) < 0.0 || eval(2, 2) < 0.0)
    {
        if (eval(0, 0) < 0.0 && eval(1, 1) < 0.0 && eval(2, 2))
        {
            return 2;
        }
        else
        {
            return 1;
        }
    }

    return 0;
}

/**
 * @brief Calculates minimum eigenvalue of Hessian matrix.
 *
 * Specialised versions of this function exist only for 2x2 and 3x3 matrices.
 */
template <int DIM> void NodeOpti::MinEigen(NekDouble &val)
{
    boost::ignore_unused(val);
    ASSERTL0(false, "DIM error");
}

template <> void NodeOpti::MinEigen<2>(NekDouble &val)
{
    NekDouble H[2][2];
    H[0][0] = m_grad[2];
    H[1][0] = m_grad[3];
    //H[0][1] = H[1][0];
    H[1][1] = m_grad[4];

    NekDouble D = (H[0][0] - H[1][1]) * (H[0][0] - H[1][1]) + 4.0 * H[1][0] * H[1][0];
    NekDouble Dsqrt = sqrt(D);

    //eval[0] = (H[0][0] + H[1][1] + Dsqrt ) / 2.0;
    val = (H[0][0] + H[1][1] - Dsqrt ) / 2.0; // the minimum Eigenvalue
}

template <> void NodeOpti::MinEigen<3>(NekDouble &val)
{
    NekDouble H[3][3];
    H[0][0] = m_grad[3];
    H[1][0] = m_grad[4];
    H[0][1] = H[1][0];
    H[2][0] = m_grad[5];
    H[0][2] = H[2][0];
    H[1][1] = m_grad[6];
    H[2][1] = m_grad[7];
    H[1][2] = H[2][1];
    H[2][2] = m_grad[8];

    //double eval[3]; // the eigenvalues

    NekDouble p1 = H[0][1] * H[0][1] + H[0][2] * H[0][2] + H[1][2] * H[1][2];
    if (p1 == 0.0) // H is diagonal
    {
        // find the minimum Eigenvalue
        if(H[0][0] < H[1][1])
        {
            if(H[0][0] < H[2][2])
            {
                val = H[0][0];
            }
            else
            {
                val = H[2][2];
            }
        }
        else
        {
            if(H[1][1] < H[2][2])
            {
                val = H[1][1];
            }
            else
            {
                val = H[2][2];
            }
        }
    }
    else
    {
        NekDouble q  = (H[0][0] + H[1][1] + H[2][2]) / 3.0;
        NekDouble p2 =    (H[0][0] - q)*(H[0][0] - q)
                     + (H[1][1] - q)*(H[1][1] - q)
                     + (H[2][2] - q)*(H[2][2] - q)
                     + 2.0 * p1;
        NekDouble p = sqrt(p2 / 6.0);

        NekDouble B[3][3];   // B = (1.0 / p) * (H - q * I)   with I being the identity matrix
        NekDouble pinv = 1.0 / p;
        B[0][0] = pinv * (H[0][0] - q);
        B[1][1] = pinv * (H[1][1] - q);
        B[2][2] = pinv * (H[2][2] - q);
        B[0][1] = pinv * H[0][1];
        B[1][0] = B[0][1];
        B[0][2] = pinv * H[0][2];
        B[2][0] = B[0][2];
        B[1][2] = pinv * H[1][2];
        B[2][1] = B[1][2];

        NekDouble r = Determinant<3>(B) / 2.0;

        // In exact arithmetic for h symmetric matrix  -1 <= r <= 1
        // but computation error can leave it slightly outside this range.
        NekDouble phi;
        if (r <= -1)
        {
            phi = M_PI / 3.0;
        }
        else if (r >= 1)
        {
            phi = 0.0;
        }
        else
        {
            phi = acos(r) / 3.0;
        }

        // the eigenvalues satisfy eval[2] <= eval[1] <= eval[0]
        //eval[0] = q + 2.0 * p * cos(phi);
        val = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
        //eval[1] = 3.0 * q - eval[0] - eval[2];     // since trace(H) = eval[0] + eval[1] + eval[2]
    }
}

}
}

#endif
