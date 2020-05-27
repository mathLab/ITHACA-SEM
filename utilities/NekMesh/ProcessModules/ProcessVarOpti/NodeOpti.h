////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_NODEOPTI
#define UTILITIES_NEKMESH_NODEOPTI

#include <ostream>
#include <mutex>

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LibUtilities/BasicUtils/Thread.h>

#include "ProcessVarOpti.h"

#include <utilities/NekMesh/ProcessModules/ProcessVarOpti/Evaluator.hxx>

namespace Nektar
{
namespace Utilities
{

class NodeOptiJob;

class NodeOpti
{
    // Typedef for derivative storage, we use boost::multi_array so we can pass
    // this to functions easily
    typedef boost::multi_array<NekDouble, 4> DerivArray;

public:
    NodeOpti(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
             ResidualSharedPtr r,
             std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
             optiType o, int dim)
        : m_node(n), m_res(r), m_derivUtils(d), m_opti(o)
    {
        // filter element types within d vector
        for (int i = 0; i < e.size(); i++)
        {
            m_data[e[i]->GetEl()->GetShapeType()].push_back(e[i]);
        }

        // Set up storage for GetFunctional to avoid reallocation on each call.
        size_t storageCount = 0;

        // Count total storage needed.
        for (auto &typeIt : m_data)
        {
            const int pts    = m_derivUtils[typeIt.first]->pts;
            const int nElmt  = typeIt.second.size();

            storageCount = std::max(storageCount,
                                    dim * m_derivUtils[typeIt.first]->ptsStd *
                                    typeIt.second.size());

            m_derivs.insert(std::make_pair(
                                typeIt.first,
                                DerivArray(boost::extents[dim][nElmt][dim][pts])));
        }

        m_tmpStore.resize(storageCount);
    }

    virtual ~NodeOpti(){};

    void CalcMinJac();

    virtual void Optimise() = 0;
    NodeOptiJob *GetJob();

    template <int DIM>
    NekDouble GetFunctional(NekDouble &minJacNew, bool gradient = true);

    template <int DIM> void MinEigen(NekDouble &val);

protected:
    NodeSharedPtr m_node;
    std::mutex mtx;
    std::map<LibUtilities::ShapeType, std::vector<ElUtilSharedPtr> > m_data;
    std::vector<NekDouble> m_grad;
    std::vector<NekDouble> m_tmpStore;
    std::unordered_map<LibUtilities::ShapeType, DerivArray, EnumHash> m_derivs;


    template <int DIM> int IsIndefinite();

    NekDouble m_minJac;
    ResidualSharedPtr m_res;
    std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> m_derivUtils;
    optiType m_opti;

    static NekDouble c1()
    {
        return 1e-3;
    }
    static NekDouble gradTol()
    {
        return 1e-8;
    }
    static NekDouble alphaTol()
    {
        return 1e-8;
    }
};

typedef std::shared_ptr<NodeOpti> NodeOptiSharedPtr;
typedef LibUtilities::NekFactory<
    int, NodeOpti, NodeSharedPtr, std::vector<ElUtilSharedPtr>,
    ResidualSharedPtr, std::map<LibUtilities::ShapeType, DerivUtilSharedPtr>,
    optiType>
    NodeOptiFactory;

NodeOptiFactory &GetNodeOptiFactory();

class NodeOpti3D3D : public NodeOpti // 1D optimsation in 3D space
{
public:
    NodeOpti3D3D(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
                 ResidualSharedPtr r,
                 std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
                 optiType o)
        : NodeOpti(n, e, r, d, o, 3)
    {
    }

    ~NodeOpti3D3D(){};

    void Optimise();

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n, std::vector<ElUtilSharedPtr> e, ResidualSharedPtr r,
        std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d, optiType o)
    {
        return NodeOptiSharedPtr(new NodeOpti3D3D(n, e, r, d, o));
    }

private:
};

class NodeOpti2D2D : public NodeOpti // 1D optimsation in 3D space
{
public:
    NodeOpti2D2D(NodeSharedPtr n, std::vector<ElUtilSharedPtr> e,
                 ResidualSharedPtr r,
                 std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d,
                 optiType o)
        : NodeOpti(n, e, r, d, o, 2)
    {
    }

    ~NodeOpti2D2D(){};

    void Optimise();

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n, std::vector<ElUtilSharedPtr> e, ResidualSharedPtr r,
        std::map<LibUtilities::ShapeType, DerivUtilSharedPtr> d, optiType o)
    {
        return NodeOptiSharedPtr(new NodeOpti2D2D(n, e, r, d, o));
    }

private:
};

class NodeOptiJob : public Thread::ThreadJob
{
public:
    NodeOptiJob(NodeOpti *no) : node(no)
    {
    }

    void Run()
    {
        node->Optimise();
    }

private:
    NodeOpti *node;
};
/**
 * @brief Evaluate functional for elements connected to a node.
 *
 * @param minJacNew   Stores current minimum Jacobian for the element group
 * @param gradient    If true, calculate gradient.
 */
template <int DIM>
NekDouble NodeOpti::GetFunctional(NekDouble &minJacNew, bool gradient)
{
    for (auto &typeIt : m_data)
    {
        const int ptsStd = m_derivUtils[typeIt.first]->ptsStd;
        const int pts    = m_derivUtils[typeIt.first]->pts;
        const int nElmt  = typeIt.second.size();

        NekDouble* X = &m_tmpStore[0];

        // Store x/y components of each element sequentially in memory
        for (int i = 0, cnt = 0; i < nElmt; ++i)
        {
            for (int j = 0; j < ptsStd; ++j)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    X[cnt + d * ptsStd + j] =
                        *(typeIt.second[i]->nodes[j][d]);
                }
            }
            cnt += DIM * ptsStd;
        }

        // Storage for derivatives, ordered by:
        //   - standard coordinate direction
        //   - number of elements
        //   - cartesian coordinate direction
        //   - quadrature points

        // Calculate x- and y-gradients
        for (int d = 0; d < DIM; ++d)
        {
            Blas::Dgemm('N', 'N', pts, DIM * nElmt, ptsStd, 1.0,
                        m_derivUtils[typeIt.first]->VdmD[d].GetRawPtr(),
                        pts, X, ptsStd, 0.0,
                        &m_derivs[typeIt.first][d][0][0][0], pts);
        }
    }

    minJacNew          = std::numeric_limits<double>::max();
    NekDouble integral = 0.0;
    NekDouble ep =
        m_minJac < 0.0 ? sqrt(1e-8 + 0.04 * m_minJac * m_minJac) : 1e-4;
    NekDouble jacIdeal[DIM][DIM], jacDet;
    m_grad = vector<NekDouble>(DIM == 2 ? 5 : 9, 0.0);

    switch (m_opti)
    {
        case eLinEl:
        {
            const NekDouble nu = 0.4;
            const NekDouble mu = 1.0 / 2.0 / (1.0 + nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            for (auto &typeIt : m_data)
            {
                const int nElmt = typeIt.second.size();
                const int pts = m_derivUtils[typeIt.first]->pts;

                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt.first]->quadW;

                for (int i = 0; i < nElmt; ++i)
                {
                    for (int k = 0; k < pts; ++k)
                    {
                        NekDouble phiM[DIM][DIM];
                        for (int l = 0; l < DIM; ++l)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                phiM[n][l] =
                                    m_derivs[typeIt.first][l][i][n][k];
                            }
                        }

                        // begin CalcIdealJac
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacIdeal[n][m] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacIdeal[n][m] += phiM[n][l] *
                                        typeIt.second[i]->maps[k][m * 3 + l];
                                }
                            }
                        }
                        jacDet = Determinant(jacIdeal);
                        // end CalcIdealJac

                        NekDouble absIdealMapDet = fabs(typeIt.second[i]->maps[k][9]);
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
                                    absIdealMapDet *
                                    (K * 0.5 * lsigma * lsigma + mu * trEtE);

                        if (gradient)
                        {
                            NekDouble jacDerivPhi[DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble derivDet = Determinant<DIM>(phiM);
                            NekDouble jacInvTrans[DIM][DIM];
                            InvTrans<DIM>(phiM, jacInvTrans);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt.first]->VdmD[m])(
                                    k, typeIt.second[i]->NodeId(m_node->m_id));
                            }
                            // jacDeriv is actually a tensor,
                            // but can be stored as a vector, as 18 out of 27 entries are zero
                            // and the other 9 entries are three triplets
                            // this is due to the delta function in jacDeriv
                            NekDouble jacDeriv[DIM];
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacDeriv[l] = basisDeriv[l];
                            }

                            // jacDerivPhi is actually a tensor,
                            // but can be stored as a vector due to the simple form of jacDeriv
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacDerivPhi[n] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacDerivPhi[n] += jacDeriv[l] *
                                            typeIt.second[i]->maps[k][l + 3 * n];
                                }
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *= derivDet / absIdealMapDet;
                            }
                            // end of common part to all four versionsNekDouble


                            NekDouble M2[DIM][DIM][DIM];
                            // use the delta function in jacDeriv and do some tensor calculus
                            // to come up with this simplified expression for:
                            // LEM2<DIM>(jacIdeal, jacDerivPhi, M2);
                            for(int d = 0; d < DIM; d++)
                            {
                                for (int m = 0; m < DIM; ++m)
                                {
                                    for (int n = 0; n < DIM; ++n)
                                    {
                                        M2[d][m][n] = 0.5*(jacDerivPhi[m] * jacIdeal[d][n]
                                                            + jacIdeal[d][m] * jacDerivPhi[n]);
                                    }
                                }
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                NekDouble frobProdA = FrobProd<DIM>(M2[m], Emat);

                                m_grad[m] +=
                                    quadW[k] * absIdealMapDet *
                                    (2.0 * mu * frobProdA +
                                     K * lsigma * jacDetDeriv[m] /
                                         (2.0 * sigma - jacDet));
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    NekDouble frobProdBC = FrobProd<DIM>(M2[m], M2[l]);
                                    NekDouble M3[DIM][DIM];
                                    // use the delta function in jacDeriv and do some tensor calculus
                                    // to come up with this simplified expression for:
                                    // LEM3<DIM>(jacDerivPhi, M3);
                                    if (m == l)
                                    {
                                        for (int p = 0; p < DIM; ++p)
                                        {
                                            for (int q = 0; q < DIM; ++q)
                                            {
                                                M3[p][q] = jacDerivPhi[p] * jacDerivPhi[q];
                                            }
                                        }
                                        frobProdBC += FrobProd<DIM>(M3,Emat);
                                    }

                                    m_grad[ct + DIM] +=
                                        quadW[k] * absIdealMapDet *
                                        (2.0 * mu * frobProdBC +
                                         jacDetDeriv[m] * jacDetDeriv[l] * K /
                                             (2.0 * sigma - jacDet) /
                                             (2.0 * sigma - jacDet) *
                                             (1.0 - jacDet * lsigma /
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
            const NekDouble nu = 0.4;
            const NekDouble mu = 1.0 / 2.0 / (1.0 + nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            for (auto &typeIt : m_data)
            {
                const int nElmt = typeIt.second.size();
                const int pts   = m_derivUtils[typeIt.first]->pts;

                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt.first]->quadW;

                for (int i = 0; i < nElmt; ++i)
                {
                    for (int k = 0; k < pts; ++k)
                    {
                        NekDouble phiM[DIM][DIM];
                        for (int l = 0; l < DIM; ++l)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                phiM[n][l] =
                                    m_derivs[typeIt.first][l][i][n][k];
                            }
                        }
                        // begin CalcIdealJac
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacIdeal[n][m] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacIdeal[n][m] += phiM[n][l] *
                                            typeIt.second[i]->maps[k][m * 3 + l];
                                }
                            }
                        }
                        jacDet = Determinant(jacIdeal);
                        // end CalcIdealJac

                        minJacNew    = min(minJacNew, jacDet);

                        NekDouble absIdealMapDet = fabs(typeIt.second[i]->maps[k][9]);

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
                        integral += quadW[k] * absIdealMapDet *
                                    (0.5 * mu * (I1 - 3.0 - 2.0 * lsigma) +
                                     0.5 * K * lsigma * lsigma);

                        // Derivative of basis function in each direction
                        if (gradient)
                        {
                            NekDouble jacDerivPhi[DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble derivDet = Determinant<DIM>(phiM);
                            NekDouble jacInvTrans[DIM][DIM];
                            InvTrans<DIM>(phiM, jacInvTrans);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] =
                                    *(m_derivUtils[typeIt.first]->VdmD[m].GetRawPtr() +
                                      typeIt.second[i]->NodeId(m_node->m_id) * pts + k);
                            }

                            // jacDeriv is actually a tensor,
                            // but can be stored as a vector, as 18 out of 27 entries are zero
                            // and the other 9 entries are three triplets
                            // this is due to the delta function in jacDeriv
                            NekDouble jacDeriv[DIM];
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacDeriv[l] = basisDeriv[l];
                            }

                            // jacDerivPhi is actually a tensor,
                            // but can be stored as a vector due to the simple form of jacDeriv
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacDerivPhi[n] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacDerivPhi[n] += jacDeriv[l] *
                                            typeIt.second[i]->maps[k][l + 3 * n];
                                }
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *= derivDet / absIdealMapDet;
                            }
                            // end of common part to all four versionsNekDouble

                            for (int m = 0; m < DIM; ++m)
                            {
                                // because of the zero entries of the tensor jacDerivPhi,
                                // the Frobenius-product becomes a scalar product
                                NekDouble frobProd =
                                    ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);

                                m_grad[m] +=
                                    quadW[k] * absIdealMapDet *
                                    (mu * frobProd +
                                     (jacDetDeriv[m] / (2.0 * sigma - jacDet) *
                                      (K * lsigma - mu)));
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    NekDouble frobProdHes = 0.0;
                                    // because of the zero entries of the tensor jacDerivPhi,
                                    // the matrix frobProdHes has only diagonal entries
                                    if (m == l)
                                    {
                                        // because of the zero entries of the tensor jacDerivPhi,
                                        // the Frobenius-product becomes a scalar product
                                        frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);
                                    }

                                    m_grad[ct + DIM] +=
                                        quadW[k] * absIdealMapDet *
                                        (mu * frobProdHes +
                                         jacDetDeriv[m] * jacDetDeriv[l] /
                                             (2.0 * sigma - jacDet) /
                                             (2.0 * sigma - jacDet) *
                                             (K - jacDet * (K * lsigma - mu) /
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
            for (auto &typeIt : m_data)
            {
                const int nElmt = typeIt.second.size();
                const int pts = m_derivUtils[typeIt.first]->pts;

                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt.first]->quadW;

                for (int i = 0; i < nElmt; ++i)
                {
                    for (int k = 0; k < pts; ++k)
                    {
                        NekDouble phiM[DIM][DIM];
                        for (int l = 0; l < DIM; ++l)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                phiM[n][l] =
                                    m_derivs[typeIt.first][l][i][n][k];
                            }
                        }
                        // begin CalcIdealJac
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacIdeal[n][m] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacIdeal[n][m] += phiM[n][l] *
                                            typeIt.second[i]->maps[k][m * 3 + l];
                                }
                            }
                        }
                        jacDet = Determinant(jacIdeal);
                        // end CalcIdealJac

                        NekDouble absIdealMapDet = fabs(typeIt.second[i]->maps[k][9]);
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
                            quadW[k] * absIdealMapDet * W;

                        // Derivative of basis function in each direction
                        if (gradient)
                        {
                            NekDouble jacDerivPhi[DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble derivDet = Determinant<DIM>(phiM);
                            NekDouble jacInvTrans[DIM][DIM];
                            InvTrans<DIM>(phiM, jacInvTrans);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt.first]->VdmD[m])(
                                    k, typeIt.second[i]->NodeId(m_node->m_id));
                            }
                            // jacDeriv is actually a tensor,
                            // but can be stored as a vector, as 18 out of 27 entries are zero
                            // and the other 9 entries are three triplets
                            // this is due to the delta function in jacDeriv
                            NekDouble jacDeriv[DIM];
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacDeriv[l] = basisDeriv[l];
                            }

                            // jacDerivPhi is actually a tensor,
                            // but can be stored as a vector due to the simple form of jacDeriv
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacDerivPhi[n] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacDerivPhi[n] += jacDeriv[l] *
                                            typeIt.second[i]->maps[k][l + 3 * n];
                                }
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *= derivDet / absIdealMapDet;
                            }
                            // end of common part to all four versionsNekDouble

                            NekDouble frobProd[DIM];
                            NekDouble inc[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                // because of the zero entries of the tensor jacDerivPhi,
                                // the Frobenius-product becomes a scalar product
                                frobProd[m] = ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);

                                inc[m] =
                                    quadW[k] * absIdealMapDet *
                                    (2.0 * W * (frobProd[m] / frob -
                                                jacDetDeriv[m] / DIM /
                                                    (2.0 * sigma - jacDet)));
                                m_grad[m] += inc[m];
                            }



                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    NekDouble frobProdHes = 0.0;
                                    // because of the zero entries of the tensor jacDerivPhi,
                                    // the matrix frobProdHes has only diagonal entries
                                    if (m == l)
                                    {
                                        // because of the zero entries of the tensor jacDerivPhi,
                                        // the Frobenius-product becomes a scalar product
                                        frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);
                                    }

                                    m_grad[ct + DIM] +=
                                        quadW[k] * absIdealMapDet *
                                        (inc[m] * inc[l] / W +
                                         2.0 * W *
                                             (frobProdHes / frob -
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
            for (auto &typeIt : m_data)
            {
                const int nElmt = typeIt.second.size();
                const int pts = m_derivUtils[typeIt.first]->pts;

                NekVector<NekDouble> &quadW =
                    m_derivUtils[typeIt.first]->quadW;

                for (int i = 0; i < nElmt; ++i)
                {
                    for (int k = 0; k < pts; ++k)
                    {
                        NekDouble phiM[DIM][DIM];
                        for (int l = 0; l < DIM; ++l)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                phiM[n][l] =
                                    m_derivs[typeIt.first][l][i][n][k];
                            }
                        }
                        // begin CalcIdealJac
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacIdeal[n][m] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacIdeal[n][m] += phiM[n][l] *
                                            typeIt.second[i]->maps[k][m * 3 + l];
                                }
                            }
                        }
                        jacDet = Determinant(jacIdeal);
                        // end CalcIdealJac

                        NekDouble absIdealMapDet = fabs(typeIt.second[i]->maps[k][9]);
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
                            quadW[k] * absIdealMapDet * W;

                        // Derivative of basis function in each direction
                        if (gradient)
                        {
                            NekDouble jacDerivPhi[DIM];
                            NekDouble jacDetDeriv[DIM];

                            NekDouble derivDet = Determinant<DIM>(phiM);
                            NekDouble jacInvTrans[DIM][DIM];
                            InvTrans<DIM>(phiM, jacInvTrans);

                            NekDouble basisDeriv[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                basisDeriv[m] = *(
                                    m_derivUtils[typeIt.first]->VdmD[m])(
                                    k, typeIt.second[i]->NodeId(m_node->m_id));
                            }
                            // jacDeriv is actually a tensor,
                            // but can be stored as a vector, as 18 out of 27 entries are zero
                            // and the other 9 entries are three triplets
                            // this is due to the delta function in jacDeriv
                            NekDouble jacDeriv[DIM];
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacDeriv[l] = basisDeriv[l];
                            }

                            // jacDerivPhi is actually a tensor,
                            // but can be stored as a vector due to the simple form of jacDeriv
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacDerivPhi[n] = 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacDerivPhi[n] += jacDeriv[l] *
                                            typeIt.second[i]->maps[k][l + 3 * n];
                                }
                            }

                            for (int m = 0; m < DIM; ++m)
                            {
                                jacDetDeriv[m] = 0.0;
                                for (int n = 0; n < DIM; ++n)
                                {
                                    jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                                }
                                jacDetDeriv[m] *= derivDet / absIdealMapDet;
                            }
                            // end of common part to all four versionsNekDouble

                            NekDouble frobProd[DIM];
                            NekDouble inc[DIM];
                            for (int m = 0; m < DIM; ++m)
                            {
                                // because of the zero entries of the tensor jacDerivPhi,
                                // the Frobenius-product becomes a scalar product
                                frobProd[m] = ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);

                                inc[m] =
                                    quadW[k] *
                                    absIdealMapDet *
                                    (W *
                                     (2.0 * frobProd[m] / frob -
                                      jacDetDeriv[m] / (2.0 * sigma - jacDet)));
                                m_grad[m] += inc[m];
                            }

                            int ct = 0;
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int l = m; l < DIM; ++l, ct++)
                                {
                                    NekDouble frobProdHes = 0.0;
                                    // because of the zero entries of the tensor jacDerivPhi,
                                    // the matrix frobProdHes has only diagonal entries
                                    if (m == l)
                                    {
                                        // because of the zero entries of the tensor jacDerivPhi,
                                        // the Frobenius-product becomes a scalar product
                                        frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);
                                    }

                                    m_grad[ct + DIM] +=
                                        quadW[k] * absIdealMapDet *
                                        (inc[m] * inc[l] / W +
                                         2.0 * W *
                                             (frobProdHes / frob -
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
