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
 * Specialised versions of this function exist only for 2x2 and 3x3 matrices.
 *
 * @param jac  Input matrix
 *
 * @return Jacobian of @p jac.
 */
template <int DIM> inline NekDouble Determinant(NekDouble jac[DIM][DIM])
{
    return 0.0;
}

template <> inline NekDouble Determinant<2>(NekDouble jac[2][2])
{
    return jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
}

template <> inline NekDouble Determinant<3>(NekDouble jac[3][3])
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
inline void InvTrans(NekDouble in[DIM][DIM], NekDouble out[DIM][DIM])
{
}

template <> inline void InvTrans<2>(NekDouble in[2][2], NekDouble out[2][2])
{
    NekDouble invDet = 1.0 / Determinant(in);

    out[0][0] =  in[1][1] * invDet;
    out[1][0] = -in[0][1] * invDet;
    out[0][1] = -in[1][0] * invDet;
    out[1][1] =  in[0][0] * invDet;
}

template <> inline void InvTrans<3>(NekDouble in[3][3], NekDouble out[3][3])
{
    NekDouble invdet = 1.0 / Determinant(in);

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
 * @brief Calculate \f$ E = F^\top F - I \f$ tensor used in derivation of linear
 * elasticity gradients.
 *
 * Specialised versions of this function exist only for 2x2 and 3x3 matrices.
 *
 * @param in   Input matrix \f$ F \f$
 * @param out  Output matrix \f$ F^\top F - I \f$
 */
template <int DIM>
inline void EMatrix(NekDouble in[DIM][DIM], NekDouble out[DIM][DIM])
{
}

template <> inline void EMatrix<2>(NekDouble in[2][2], NekDouble out[2][2])
{
    out[0][0] = 0.5 * (in[0][0] * in[0][0] + in[1][0] * in[1][0] - 1.0);
    out[1][0] = 0.5 * (in[0][0] * in[0][1] + in[1][0] * in[1][1]);
    out[0][1] = 0.5 * (in[0][0] * in[0][1] + in[1][0] * in[1][1]);
    out[1][1] = 0.5 * (in[1][1] * in[1][1] + in[0][1] * in[0][1] - 1.0);
}

template <> inline void EMatrix<3>(NekDouble in[3][3], NekDouble out[3][3])
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
 * @brief Auxiliary function used in the calculation of linear elasticity
 * gradients.
 */
template <int DIM>
inline void LEM2(NekDouble jacIdeal[DIM][DIM],
                 NekDouble jacDerivPhi[DIM][DIM][DIM],
                 NekDouble ret[DIM][DIM][DIM])
{
    for (int k = 0; k < DIM; k++)
    {
        NekDouble part1[DIM][DIM], part2[DIM][DIM];
        for (int m = 0; m < DIM; ++m)
        {
            for (int n = 0; n < DIM; ++n)
            {
                part1[m][n] = 0.0;
                part2[m][n] = 0.0;
                for (int l = 0; l < DIM; ++l)
                {
                    part1[m][n] += jacDerivPhi[k][l][m] * jacIdeal[l][n];
                    part2[m][n] += jacIdeal[l][m] * jacDerivPhi[k][l][n];
                }
            }
        }

        for (int m = 0; m < DIM; ++m)
        {
            for (int n = 0; n < DIM; ++n)
            {
                ret[k][m][n] = 0.5 * (part1[m][n] + part2[m][n]);
            }
        }
    }
}

/**
 * @brief Auxiliary function used in the calculation of linear elasticity
 * gradients.
 */
template <int DIM>
inline void LEM3(NekDouble jacDerivPhi[DIM][DIM][DIM],
                 NekDouble ret[DIM][DIM][DIM][DIM])
{
    for (int j = 0; j < DIM; j++)
    {
        for (int k = 0; k < DIM; k++)
        {
            NekDouble part1[DIM][DIM], part2[DIM][DIM];
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    part1[m][n] = 0.0;
                    part2[m][n] = 0.0;
                    for (int l = 0; l < DIM; ++l)
                    {
                        part1[m][n] +=
                            jacDerivPhi[j][l][m] * jacDerivPhi[k][l][n];
                        part2[m][n] +=
                            jacDerivPhi[k][l][m] * jacDerivPhi[j][l][n];
                    }
                }
            }

            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    ret[j][k][m][n] = 0.5 * (part1[m][n] + part2[m][n]);
                }
            }
        }
    }
}

/**
 * @brief Calculate Frobenius inner product of input matrices.
 */
template <int DIM>
inline NekDouble FrobProd(NekDouble in1[DIM][DIM], NekDouble in2[DIM][DIM])
{
    NekDouble ret = 0;
    for (int n = 0; n < DIM; ++n)
    {
        for (int l = 0; l < DIM; ++l)
        {
            ret += in1[n][l] * in2[n][l];
        }
    }
    return ret;
}

// Typedef for derivative storage, we use boost::multi_array so we can pass this
// to functions easily
typedef boost::multi_array<NekDouble, 4> DerivArray;

/**
 * @brief Calculate Jacobian matrix \f$ \nabla\phi =
 * \nabla\phi_M\nabla\phi_I^{-1} \f$ for each evaluation point of each element.
 *
 * @param elmt       Element to process
 * @param point      Index of evaluation point
 * @param deriv      Derivative array containing \f$ \nabla\phi_M \f$
 * @param data       Data array containing \f$ \nabla\phi_I^{-1} \f$
 * @param jacIdeal   Output Jacobian matrix
 */
template <int DIM>
inline NekDouble CalcIdealJac(int elmt, int point, DerivArray &deriv,
                              std::vector<ElUtilSharedPtr> &data,
                              NekDouble jacIdeal[DIM][DIM])
{
    for (int m = 0; m < DIM; ++m)
    {
        for (int n = 0; n < DIM; ++n)
        {
            jacIdeal[n][m] = 0.0;
            for (int l = 0; l < DIM; ++l)
            {
                jacIdeal[n][m] += deriv[l][elmt][n][point] *
                                  data[elmt]->maps[point][m * 3 + l];
            }
        }
    }

    return Determinant(jacIdeal);
}

/**
 * @brief Calculate Frobenius norm \f$ \| A \|_f ^2 \f$ of a matrix \f$ A \f$.
 *
 * @param inarray   Input matrix \f$ A \f$
 */
template <int DIM> inline NekDouble FrobeniusNorm(NekDouble inarray[DIM][DIM])
{
    NekDouble ret = 0.0, *start = &inarray[0][0];
    for (int i = 0; i < DIM * DIM; ++i, ++start)
    {
        ret += (*start) * (*start);
    }
    return ret;
}

/**
 * @brief Evaluate functional for elements connected to a node.
 *
 * @param minJacNew   Stores current minimum Jacobian for the element group
 * @param gradient    If true, calculate gradient.
 */
template <int DIM>
NekDouble NodeOpti::GetFunctional(NekDouble &minJacNew, bool gradient)
{
    map<LibUtilities::ShapeType, vector<ElUtilSharedPtr> >::iterator typeIt;
    map<LibUtilities::ShapeType, DerivArray> derivs;

    for (typeIt = m_data.begin(); typeIt != m_data.end(); typeIt++)
    {
        const int nElmt  = typeIt->second.size();
        const int totpts = m_derivUtils[typeIt->first]->ptsStd * nElmt;
        NekDouble X[DIM * totpts];

        // Store x/y components of each element sequentially in memory
        for (int i = 0, cnt = 0; i < nElmt; ++i)
        {
            for (int j = 0; j < m_derivUtils[typeIt->first]->ptsStd; ++j)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    X[cnt + d * m_derivUtils[typeIt->first]->ptsStd + j] =
                        *(typeIt->second[i]->nodes[j][d]);
                }
            }
            cnt += DIM * m_derivUtils[typeIt->first]->ptsStd;
        }

        // Storage for derivatives, ordered by:
        //   - standard coordinate direction
        //   - number of elements
        //   - cartesian coordinate direction
        //   - quadrature points
        derivs.insert(std::make_pair(
            typeIt->first,
            DerivArray(boost::extents[DIM][nElmt][DIM]
                                     [m_derivUtils[typeIt->first]->pts])));

        // Calculate x- and y-gradients
        for (int d = 0; d < DIM; ++d)
        {
            Blas::Dgemm('N', 'N', m_derivUtils[typeIt->first]->pts, DIM * nElmt,
                        m_derivUtils[typeIt->first]->ptsStd, 1.0,
                        m_derivUtils[typeIt->first]->VdmD[d].GetRawPtr(),
                        m_derivUtils[typeIt->first]->pts, X,
                        m_derivUtils[typeIt->first]->ptsStd, 0.0,
                        &derivs[typeIt->first][d][0][0][0],
                        m_derivUtils[typeIt->first]->pts);
        }
    }

    minJacNew          = std::numeric_limits<double>::max();
    NekDouble integral = 0.0;
    NekDouble ep =
        m_minJac < 0.0 ? sqrt(1e-8 + 0.04 * m_minJac * m_minJac) : 1e-4;
    NekDouble jacIdeal[DIM][DIM], jacDet;
    m_grad = Array<OneD, NekDouble>(DIM == 2 ? 5 : 9, 0.0);

    switch (m_opti)
    {
        case eLinEl:
        {
            const NekDouble nu = 0.45;
            const NekDouble mu = 1.0 / 2.0 / (1.0 + nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            for (typeIt = m_data.begin(); typeIt != m_data.end(); typeIt++)
            {
                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt->first]->quadW;
                for (int i = 0; i < typeIt->second.size(); ++i)
                {
                    for (int k = 0; k < m_derivUtils[typeIt->first]->pts; ++k)
                    {
                        jacDet = CalcIdealJac(i, k, derivs[typeIt->first],
                                              typeIt->second, jacIdeal);
                        minJacNew = min(minJacNew, jacDet);

                        NekDouble Emat[DIM][DIM];
                        EMatrix<DIM>(jacIdeal, Emat);

                        NekDouble trEtE = FrobProd<DIM>(Emat, Emat);
                        NekDouble sigma =
                            0.5 *
                            (jacDet + sqrt(jacDet * jacDet + 4.0 * ep * ep));

                        if (sigma < numeric_limits<double>::min() && !gradient)
                        {
                            return numeric_limits<double>::max();
                        }
                        ASSERTL0(sigma > numeric_limits<double>::min(),
                                 std::string("dividing by zero ") +
                                     boost::lexical_cast<string>(sigma) + " " +
                                     boost::lexical_cast<string>(jacDet) + " " +
                                     boost::lexical_cast<string>(ep));

                        NekDouble lsigma = log(sigma);
                        integral += quadW[k] *
                                    fabs(typeIt->second[i]->maps[k][9]) *
                                    (K * 0.5 * lsigma * lsigma + mu * trEtE);

                        if (gradient)
                        {
                            NekDouble jacInvTrans[DIM][DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble phiM[DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    phiM[n][m] =
                                        derivs[typeIt->first][m][i][n][k];
                                }
                            }

                            InvTrans<DIM>(phiM, jacInvTrans);
                            NekDouble derivDet = Determinant<DIM>(phiM);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt->first]->VdmD[m])(
                                    k, typeIt->second[i]->NodeId(m_node->m_id));
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] +=
                                        jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *=
                                    derivDet /
                                    fabs(typeIt->second[i]->maps[k][9]);
                            }

                            NekDouble jacDeriv[DIM][DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    NekDouble delta = m == n ? 1.0 : 0.0;
                                    for (int l = 0; l < DIM; ++l)
                                    {
                                        jacDeriv[m][n][l] =
                                            delta * basisDeriv[l];
                                    }
                                }
                            }

                            NekDouble jacDerivPhi[DIM][DIM][DIM];
                            for (int p = 0; p < DIM; ++p)
                            {
                                for (int m = 0; m < DIM; ++m)
                                {
                                    for (int n = 0; n < DIM; ++n)
                                    {
                                        jacDerivPhi[p][m][n] = 0.0;
                                        for (int l = 0; l < DIM; ++l)
                                        {
                                            // want phi_I^{-1} (l,n)
                                            jacDerivPhi[p][m][n] +=
                                                jacDeriv[p][m][l] *
                                                typeIt->second[i]
                                                    ->maps[k][l + 3 * n];
                                        }
                                    }
                                }
                            }

                            NekDouble M2[DIM][DIM][DIM];
                            LEM2<DIM>(jacIdeal, jacDerivPhi, M2);

                            NekDouble M3[DIM][DIM][DIM][DIM];
                            LEM3<DIM>(jacDerivPhi, M3);

                            NekDouble frobProdA[DIM];
                            NekDouble frobProdB[DIM][DIM];
                            NekDouble frobProdC[DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                frobProdA[m] = FrobProd<DIM>(M2[m], Emat);
                                for (int l = 0; l < DIM; l++)
                                {
                                    frobProdB[m][l] =
                                        FrobProd<DIM>(M3[m][l], Emat);
                                    frobProdC[m][l] =
                                        FrobProd<DIM>(M2[m], M2[l]);
                                }
                            }

                            for (int j = 0; j < DIM; ++j)
                            {
                                m_grad[j] +=
                                    quadW[k] *
                                    fabs(typeIt->second[i]->maps[k][9]) *
                                    (2.0 * mu * frobProdA[j] +
                                     K * lsigma * jacDetDeriv[j] /
                                         (2.0 * sigma - jacDet));
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    m_grad[ct + DIM] +=
                                        quadW[k] *
                                        fabs(typeIt->second[i]->maps[k][9]) *
                                        (2.0 * mu * frobProdB[m][l] +
                                         2.0 * mu * frobProdC[m][l] +
                                         jacDetDeriv[m] * jacDetDeriv[l] * K /
                                             (2.0 * sigma - jacDet) /
                                             (2.0 * sigma - jacDet) *
                                             (1.0 -
                                              jacDet * lsigma /
                                                  (2.0 * sigma - jacDet)));
                                }
                            }
                        }
                    }
                }
            }
            break;
        }

        case eHypEl:
        {
            const NekDouble nu = 0.45;
            const NekDouble mu = 1.0 / 2.0 / (1.0 + nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            for (typeIt = m_data.begin(); typeIt != m_data.end(); typeIt++)
            {
                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt->first]->quadW;
                for (int i = 0; i < typeIt->second.size(); ++i)
                {
                    for (int k = 0; k < m_derivUtils[typeIt->first]->pts; ++k)
                    {
                        jacDet = CalcIdealJac(i, k, derivs[typeIt->first],
                                              typeIt->second, jacIdeal);
                        minJacNew    = min(minJacNew, jacDet);
                        NekDouble I1 = FrobeniusNorm(jacIdeal);

                        NekDouble sigma =
                            0.5 *
                            (jacDet + sqrt(jacDet * jacDet + 4.0 * ep * ep));

                        if (sigma < numeric_limits<double>::min() && !gradient)
                        {
                            return numeric_limits<double>::max();
                        }

                        ASSERTL0(sigma > numeric_limits<double>::min(),
                                 std::string("dividing by zero ") +
                                     boost::lexical_cast<string>(sigma) + " " +
                                     boost::lexical_cast<string>(jacDet) + " " +
                                     boost::lexical_cast<string>(ep));

                        NekDouble lsigma = log(sigma);
                        integral += quadW[k] *
                                    fabs(typeIt->second[i]->maps[k][9]) *
                                    (0.5 * mu * (I1 - 3.0 - 2.0 * lsigma) +
                                     0.5 * K * lsigma * lsigma);

                        // Derivative of basis function in each direction
                        if (gradient)
                        {
                            NekDouble jacInvTrans[DIM][DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble phiM[DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    phiM[n][m] =
                                        derivs[typeIt->first][m][i][n][k];
                                }
                            }

                            InvTrans<DIM>(phiM, jacInvTrans);
                            NekDouble derivDet = Determinant<DIM>(phiM);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt->first]->VdmD[m])(
                                    k, typeIt->second[i]->NodeId(m_node->m_id));
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] +=
                                        jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *=
                                    derivDet /
                                    fabs(typeIt->second[i]->maps[k][9]);
                            }

                            NekDouble jacDeriv[DIM][DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    NekDouble delta = m == n ? 1.0 : 0.0;
                                    for (int l = 0; l < DIM; ++l)
                                    {
                                        jacDeriv[m][n][l] =
                                            delta * basisDeriv[l];
                                    }
                                }
                            }

                            NekDouble jacDerivPhi[DIM][DIM][DIM];
                            for (int p = 0; p < DIM; ++p)
                            {
                                for (int m = 0; m < DIM; ++m)
                                {
                                    for (int n = 0; n < DIM; ++n)
                                    {
                                        jacDerivPhi[p][m][n] = 0.0;
                                        for (int l = 0; l < DIM; ++l)
                                        {
                                            // want phi_I^{-1} (l,n)
                                            jacDerivPhi[p][m][n] +=
                                                jacDeriv[p][m][l] *
                                                typeIt->second[i]
                                                    ->maps[k][l + 3 * n];
                                        }
                                    }
                                }
                            }

                            NekDouble frobProd[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                frobProd[m] =
                                    FrobProd<DIM>(jacIdeal, jacDerivPhi[m]);
                            }

                            for (int j = 0; j < DIM; ++j)
                            {
                                m_grad[j] +=
                                    quadW[k] *
                                    fabs(typeIt->second[i]->maps[k][9]) *
                                    (mu * frobProd[j] +
                                     (jacDetDeriv[j] / (2.0 * sigma - jacDet) *
                                      (K * lsigma - mu)));
                            }

                            NekDouble frobProdHes[DIM][DIM]; // holder for the
                                                             // hessian
                                                             // frobprods
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l)
                                {
                                    frobProdHes[m][l] = FrobProd<DIM>(
                                        jacDerivPhi[m], jacDerivPhi[l]);
                                }
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    m_grad[ct + DIM] +=
                                        quadW[k] *
                                        fabs(typeIt->second[i]->maps[k][9]) *
                                        (mu * frobProdHes[m][l] +
                                         jacDetDeriv[m] * jacDetDeriv[l] /
                                             (2.0 * sigma - jacDet) /
                                             (2.0 * sigma - jacDet) *
                                             (K -
                                              jacDet * (K * lsigma - mu) /
                                                  (2.0 * sigma - jacDet)));
                                }
                            }
                        }
                    }
                }
            }
            break;
        }

        case eRoca:
        {
            for (typeIt = m_data.begin(); typeIt != m_data.end(); typeIt++)
            {
                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt->first]->quadW;
                for (int i = 0; i < typeIt->second.size(); ++i)
                {
                    for (int k = 0; k < m_derivUtils[typeIt->first]->pts; ++k)
                    {
                        jacDet = CalcIdealJac(i, k, derivs[typeIt->first],
                                              typeIt->second, jacIdeal);
                        minJacNew      = min(minJacNew, jacDet);
                        NekDouble frob = FrobeniusNorm(jacIdeal);
                        NekDouble sigma =
                            0.5 *
                            (jacDet + sqrt(jacDet * jacDet + 4.0 * ep * ep));

                        if (sigma < numeric_limits<double>::min() && !gradient)
                        {
                            return numeric_limits<double>::max();
                        }

                        ASSERTL0(sigma > numeric_limits<double>::min(),
                                 std::string("dividing by zero ") +
                                     boost::lexical_cast<string>(sigma) + " " +
                                     boost::lexical_cast<string>(jacDet) + " " +
                                     boost::lexical_cast<string>(ep));

                        NekDouble W = frob / DIM / pow(fabs(sigma), 2.0 / DIM);
                        integral +=
                            quadW[k] * fabs(typeIt->second[i]->maps[k][9]) * W;

                        // Derivative of basis function in each direction
                        if (gradient)
                        {
                            NekDouble jacInvTrans[DIM][DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble phiM[DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    phiM[n][m] =
                                        derivs[typeIt->first][m][i][n][k];
                                }
                            }

                            InvTrans<DIM>(phiM, jacInvTrans);
                            NekDouble derivDet = Determinant<DIM>(phiM);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt->first]->VdmD[m])(
                                    k, typeIt->second[i]->NodeId(m_node->m_id));
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] +=
                                        jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *=
                                    derivDet /
                                    fabs(typeIt->second[i]->maps[k][9]);
                            }

                            NekDouble jacDeriv[DIM][DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    NekDouble delta = m == n ? 1.0 : 0.0;
                                    for (int l = 0; l < DIM; ++l)
                                    {
                                        jacDeriv[m][n][l] =
                                            delta * basisDeriv[l];
                                    }
                                }
                            }

                            NekDouble jacDerivPhi[DIM][DIM][DIM];
                            for (int p = 0; p < DIM; ++p)
                            {
                                for (int m = 0; m < DIM; ++m)
                                {
                                    for (int n = 0; n < DIM; ++n)
                                    {
                                        jacDerivPhi[p][m][n] = 0.0;
                                        for (int l = 0; l < DIM; ++l)
                                        {
                                            // want phi_I^{-1} (l,n)
                                            jacDerivPhi[p][m][n] +=
                                                jacDeriv[p][m][l] *
                                                typeIt->second[i]
                                                    ->maps[k][l + 3 * n];
                                        }
                                    }
                                }
                            }

                            NekDouble frobProd[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                frobProd[m] =
                                    FrobProd<DIM>(jacIdeal, jacDerivPhi[m]);
                            }

                            for (int j = 0; j < DIM; ++j)
                            {
                                m_grad[j] +=
                                    quadW[k] *
                                    fabs(typeIt->second[i]->maps[k][9]) *
                                    (2.0 * W * (frobProd[j] / frob -
                                                jacDetDeriv[j] / DIM /
                                                    (2.0 * sigma - jacDet)));
                            }

                            NekDouble frobProdHes[DIM][DIM]; // holder for the
                                                             // hessian
                                                             // frobprods
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l)
                                {
                                    frobProdHes[m][l] = FrobProd<DIM>(
                                        jacDerivPhi[m], jacDerivPhi[l]);
                                }
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    m_grad[ct + DIM] +=
                                        quadW[k] *
                                        fabs(typeIt->second[i]->maps[k][9]) *
                                        (m_grad[m] * m_grad[l] / W +
                                         2.0 * W *
                                             (frobProdHes[m][l] / frob -
                                              2.0 * frobProd[m] * frobProd[l] /
                                                  frob / frob +
                                              jacDetDeriv[m] * jacDetDeriv[l] *
                                                  jacDet /
                                                  (2.0 * sigma - jacDet) /
                                                  (2.0 * sigma - jacDet) /
                                                  (2.0 * sigma - jacDet) /
                                                  DIM));
                                }
                            }
                        }
                    }
                }
            }
            break;
        }

        case eWins:
        {
            for (typeIt = m_data.begin(); typeIt != m_data.end(); typeIt++)
            {
                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt->first]->quadW;
                for (int i = 0; i < typeIt->second.size(); ++i)
                {
                    for (int k = 0; k < m_derivUtils[typeIt->first]->pts; ++k)
                    {
                        jacDet = CalcIdealJac(i, k, derivs[typeIt->first],
                                              typeIt->second, jacIdeal);
                        minJacNew      = min(minJacNew, jacDet);
                        NekDouble frob = FrobeniusNorm(jacIdeal);
                        NekDouble sigma =
                            0.5 *
                            (jacDet + sqrt(jacDet * jacDet + 4.0 * ep * ep));

                        if (sigma < numeric_limits<double>::min() && !gradient)
                        {
                            return numeric_limits<double>::max();
                        }

                        ASSERTL0(sigma > numeric_limits<double>::min(),
                                 std::string("dividing by zero ") +
                                     boost::lexical_cast<string>(sigma) + " " +
                                     boost::lexical_cast<string>(jacDet) + " " +
                                     boost::lexical_cast<string>(ep));

                        NekDouble W = frob / sigma;
                        integral +=
                            quadW[k] * fabs(typeIt->second[i]->maps[k][9]) * W;

                        // Derivative of basis function in each direction
                        if (gradient)
                        {
                            NekDouble jacInvTrans[DIM][DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble phiM[DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    phiM[n][m] =
                                        derivs[typeIt->first][m][i][n][k];
                                }
                            }

                            InvTrans<DIM>(phiM, jacInvTrans);
                            NekDouble derivDet = Determinant<DIM>(phiM);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt->first]->VdmD[m])(
                                    k, typeIt->second[i]->NodeId(m_node->m_id));
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] +=
                                        jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *=
                                    derivDet /
                                    fabs(typeIt->second[i]->maps[k][9]);
                            }

                            NekDouble jacDeriv[DIM][DIM][DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    NekDouble delta = m == n ? 1.0 : 0.0;
                                    for (int l = 0; l < DIM; ++l)
                                    {
                                        jacDeriv[m][n][l] =
                                            delta * basisDeriv[l];
                                    }
                                }
                            }

                            NekDouble jacDerivPhi[DIM][DIM][DIM];
                            for (int p = 0; p < DIM; ++p)
                            {
                                for (int m = 0; m < DIM; ++m)
                                {
                                    for (int n = 0; n < DIM; ++n)
                                    {
                                        jacDerivPhi[p][m][n] = 0.0;
                                        for (int l = 0; l < DIM; ++l)
                                        {
                                            // want phi_I^{-1} (l,n)
                                            jacDerivPhi[p][m][n] +=
                                                jacDeriv[p][m][l] *
                                                typeIt->second[i]
                                                    ->maps[k][l + 3 * n];
                                        }
                                    }
                                }
                            }

                            NekDouble frobProd[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                frobProd[m] =
                                    FrobProd<DIM>(jacIdeal, jacDerivPhi[m]);
                            }

                            for (int j = 0; j < DIM; ++j)
                            {
                                m_grad[j] +=
                                    quadW[k] *
                                    fabs(typeIt->second[i]->maps[k][9]) *
                                    (W *
                                     (2.0 * frobProd[j] / frob -
                                      jacDetDeriv[j] / (2.0 * sigma - jacDet)));
                            }

                            NekDouble frobProdHes[DIM][DIM]; // holder for the
                                                             // hessian
                                                             // frobprods
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l)
                                {
                                    frobProdHes[m][l] = FrobProd<DIM>(
                                        jacDerivPhi[m], jacDerivPhi[l]);
                                }
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    m_grad[ct + DIM] +=
                                        quadW[k] *
                                        fabs(typeIt->second[i]->maps[k][9]) *
                                        (m_grad[m] * m_grad[l] / W +
                                         2.0 * W *
                                             (frobProdHes[m][l] / frob -
                                              2.0 * frobProd[m] * frobProd[l] /
                                                  frob / frob +
                                              0.5 * jacDetDeriv[m] *
                                                  jacDetDeriv[l] * jacDet /
                                                  (2.0 * sigma - jacDet) /
                                                  (2.0 * sigma - jacDet) /
                                                  (2.0 * sigma - jacDet)));
                                }
                            }
                        }
                    }
                }
            }
            break;
        }
    }

    // ASSERTL0(std::isfinite(integral),"inf in integral");

    return integral;
    // return sqrt(m_grad[0]*m_grad[0] + m_grad[1]*m_grad[1]);
}
}
}

#endif
