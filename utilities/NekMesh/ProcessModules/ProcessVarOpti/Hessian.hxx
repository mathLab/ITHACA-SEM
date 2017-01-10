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
//  License for the specific language governing rights and limitations under
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
    Array<OneD, NekDouble> eigR(2);
    Array<OneD, NekDouble> eigI(2);
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
    Array<OneD, NekDouble> vl(nVel * nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi(nVel);

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
    Array<OneD, NekDouble> eigR(3);
    Array<OneD, NekDouble> eigI(3);
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
    Array<OneD, NekDouble> vl(nVel * nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi(nVel);

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
    ASSERTL0(false, "DIM error");
}

template <> void NodeOpti::MinEigen<2>(NekDouble &val)
{
    Array<OneD, NekDouble> eigR(2);
    Array<OneD, NekDouble> eigI(2);
    NekMatrix<NekDouble> H(2, 2);
    H(0, 0) = m_grad[2];
    H(1, 0) = m_grad[3];
    H(0, 1) = H(1, 0);
    H(1, 1) = m_grad[4];

    int nVel   = 2;
    char jobvl = 'N', jobvr = 'V';
    int worklen = 8 * nVel, info;

    DNekMat eval(nVel, nVel, 0.0, eDIAGONAL);
    DNekMat evec(nVel, nVel, 0.0, eFULL);
    Array<OneD, NekDouble> vl(nVel * nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi(nVel);

    Lapack::Dgeev(jobvl, jobvr, nVel, H.GetRawPtr(), nVel, &(eval.GetPtr())[0],
                  &wi[0], &vl[0], nVel, &(evec.GetPtr())[0], nVel, &work[0],
                  worklen, info);

    ASSERTL0(!info, "dgeev failed");

    int minI;
    NekDouble tmp = std::numeric_limits<double>::max();
    for (int i = 0; i < 2; i++)
    {
        if (eval(i, i) < tmp)
        {
            minI = i;
            tmp  = eval(i, i);
        }
    }

    val = eval(minI, minI);
}

template <> void NodeOpti::MinEigen<3>(NekDouble &val)
{
    Array<OneD, NekDouble> eigR(3);
    Array<OneD, NekDouble> eigI(3);
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
    char jobvl = 'N', jobvr = 'V';
    int worklen = 8 * nVel, info;

    DNekMat eval(nVel, nVel, 0.0, eDIAGONAL);
    DNekMat evec(nVel, nVel, 0.0, eFULL);
    Array<OneD, NekDouble> vl(nVel * nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi(nVel);

    Lapack::Dgeev(jobvl, jobvr, nVel, H.GetRawPtr(), nVel, &(eval.GetPtr())[0],
                  &wi[0], &vl[0], nVel, &(evec.GetPtr())[0], nVel, &work[0],
                  worklen, info);

    ASSERTL0(!info, "dgeev failed");

    int minI;
    NekDouble tmp = std::numeric_limits<double>::max();
    for (int i = 0; i < 3; i++)
    {
        if (eval(i, i) < tmp)
        {
            minI = i;
            tmp  = eval(i, i);
        }
    }

    val = eval(minI, minI);
}

}
}

#endif
