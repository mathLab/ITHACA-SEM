///////////////////////////////////////////////////////////////////////////////
//
// File NodalUtil.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: 2D Nodal Triangle Fekete Utilities --
//              Basis function, Interpolation, Integral, Derivation, etc
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <limits>

#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/NodalUtil.h>

#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

namespace Nektar
{
namespace LibUtilities
{

/**
 * @brief Obtain the integration weights for the given nodal distribution.
 *
 * This routine constructs the integration weights for the given nodal
 * distribution inside NodalUtil::m_xi. To do this we solve the linear system
 * \f$ \mathbf{V}^\top \mathbf{w} = \mathbf{g} \f$ where \f$ \mathbf{g}_i =
 * \int_\Omega \psi_i(x)\, dx \f$, and we use the fact that under the definition
 * of the orthogonal basis for each element type, \f$ \mathbf{g}_i = 0 \f$ for
 * \f$ i > 0 \f$. We use NodalUtil::v_ModeZeroIntegral to return the analytic
 * value of \f$ \mathbf{g}_0 \f$.
 *
 * @return Vector of integration weights for the integration points.
 */
NekVector<NekDouble> NodalUtil::GetWeights()
{
    // Get number of modes in orthogonal basis
    int numModes = v_NumModes();

    // If we have the same number of nodes as modes, then we can solve the
    // linear system V^T w = (1,0,...)
    if (numModes == m_numPoints)
    {
        NekVector<NekDouble> g(m_numPoints, 0.0);
        g(0) = v_ModeZeroIntegral();

        SharedMatrix V = GetVandermonde();

        // Solve the system V^T w = g to obtain weights.
        LinearSystem matL(V);
        return matL.SolveTranspose(g);
    }
    else
    {
        // System is either over- or under-determined. Need to do least squares
        // here using SVD.
        return NekVector<NekDouble>();
    }
}

/**
 * @brief Return the Vandermonde matrix for the nodal distribution.
 *
 * This routine constructs and returns the Vandermonde matrix \f$\mathbf{V}\f$,
 * with each entry as \f$\mathbf{V}_{ij} = (\psi_i(\xi_j))\f$ where \f$ \psi_i
 * \f$ is the orthogonal basis obtained through the abstract function
 * NodalUtil::v_OrthoBasis.
 *
 * @return The Vandermonde matrix.
 */
SharedMatrix NodalUtil::GetVandermonde()
{
    int rows = m_numPoints, cols = v_NumModes();
    SharedMatrix matV = MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(
        rows, cols, 0.0);

    for (int j = 0; j < cols; ++j)
    {
        NekVector<NekDouble> col = v_OrthoBasis(j);
        for (int i = 0; i < rows; ++i)
        {
            (*matV)(i, j) = col[i];
        }
    }

    return matV;
}
/**
 * @brief Return the Vandermonde matrix of the derivative of the basis functions
 * for the nodal distribution.
 *
 * This routine constructs and returns the Vandermonde matrix for the derivative
 * of the basis functions \f$\mathbf{V}_d\f$ for coordinate directions \f$ 0
 * \leq d \leq 2 \f$, with each entry as \f$\mathbf{V}_{ij} = (\partial_d
 * \psi_i(\xi_j))\f$ where \f$ \partial_d\psi_i \f$ is the derivative of the
 * orthogonal basis obtained through the abstract function
 * NodalUtil::v_OrthoBasisDeriv.
 *
 * @param dir  Direction of derivative in the standard element.
 *
 * @return Vandermonde matrix corresponding with derivative of the basis
 *         functions in direction @p dir.
 */
SharedMatrix NodalUtil::GetVandermondeForDeriv(int dir)
{
    int rows = m_numPoints, cols = v_NumModes();
    SharedMatrix matV = MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(
        rows, cols, 0.0);

    for (int j = 0; j < cols; ++j)
    {
        NekVector<NekDouble> col = v_OrthoBasisDeriv(dir, j);
        for (int i = 0; i < rows; ++i)
        {
            (*matV)(i, j) = col[i];
        }
    }

    return matV;
}

/**
 * @brief Return the derivative matrix for the nodal distribution.
 *
 * This routine constructs and returns the derivative matrices
 * \f$\mathbf{D}_d\f$ for coordinate directions \f$ 0 \leq d \leq 2 \f$, which
 * can be used to evaluate the derivative of a nodal expansion at the points
 * defined by NodalUtil::m_xi. These are calculated as \f$ \mathbf{D}_d =
 * \mathbf{V}_d \mathbf{V}^{-1} \f$, where \f$ \mathbf{V}_d \f$ is the
 * derivative Vandermonde matrix and \f$ \mathbf{V} \f$ is the Vandermonde
 * matrix.
 *
 * @param dir  Coordinate direction in which to evaluate the derivative.
 *
 * @return The derivative matrix for direction @p dir.
 */
SharedMatrix NodalUtil::GetDerivMatrix(int dir)
{
    SharedMatrix V  = GetVandermonde();
    SharedMatrix Vd = GetVandermondeForDeriv(dir);
    SharedMatrix D = MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(
        V->GetRows(), V->GetColumns(), 0.0);

    V->Invert();

    *D = (*Vd) * (*V);

    return D;
}

/**
 * @brief Construct the interpolation matrix used to evaluate the basis at the
 * points @p xi inside the element.
 *
 * This routine returns a matrix \f$ \mathbf{I}(\mathbf{a}) \f$ that can be used
 * to evaluate the nodal basis at the points defined by the parameter @p xi,
 * which is denoted by \f$ \mathbf{a} = (a_1, \dots, a_N) \f$ and \f$ N \f$ is
 * the number of points in @p xi.
 *
 * In particular, if the array \f$ \mathbf{u} \f$ with components \f$
 * \mathbf{u}_i = u^\delta(\xi_i) \f$ represents the polynomial approximation of
 * a function \f$ u \f$ evaluated at the nodal points NodalUtil::m_xi, then the
 * evaluation of \f$ u^\delta \f$ evaluated at the input points \f$ \mathbf{a}
 * \f$ is given by \f$ \mathbf{I}(\mathbf{a})\mathbf{u} \f$.
 *
 * @param xi  An array of first size number of spatial dimensions \f$ d \f$ and
 *            secondary size the number of points to interpolate.
 *
 * @return The interpolation matrix for the points @p xi.
 */
SharedMatrix NodalUtil::GetInterpolationMatrix(
    Array<OneD, Array<OneD, NekDouble> > &xi)
{
    std::shared_ptr<NodalUtil> subUtil = v_CreateUtil(xi);
    SharedMatrix matS = GetVandermonde();
    SharedMatrix matT = subUtil->GetVandermonde();
    SharedMatrix D = MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(
        matT->GetRows(), matS->GetColumns(), 0.0);

    matS->Invert();

    *D = (*matT) * (*matS);
    return D;
}

/**
 * @brief Construct the nodal utility class for a triangle.
 *
 * The constructor of this class sets up two member variables used in the
 * evaluation of the orthogonal basis:
 *
 * - NodalUtilTriangle::m_eta is used to construct the collapsed coordinate
 *   locations of the nodal points \f$ (\eta_1, \eta_2) \f$ inside the square
 *   \f$[-1,1]^2\f$ on which the orthogonal basis functions are defined.
 * - NodalUtilTriangle::m_ordering constructs a mapping from the index set \f$ I
 *   = \{ (i,j)\ |\ 0\leq i,j \leq P, i+j \leq P \}\f$ to an ordering \f$ 0 \leq
 *   m(ij) \leq (P+1)(P+2)/2 \f$ that defines the monomials \f$ \xi_1^i \xi_2^j
 *   \f$ that span the triangular space. This is then used to calculate which
 *   \f$ (i,j) \f$ pair corresponding to a column of the Vandermonde matrix when
 *   calculating the orthogonal polynomials.
 *
 * @param degree  Polynomial order of this nodal triangle.
 * @param r       \f$ \xi_1 \f$-coordinates of nodal points in the standard
 *                element.
 * @param s       \f$ \xi_2 \f$-coordinates of nodal points in the standard
 *                element.
 */
NodalUtilTriangle::NodalUtilTriangle(int                    degree,
                                     Array<OneD, NekDouble> r,
                                     Array<OneD, NekDouble> s)
    : NodalUtil(degree, 2), m_eta(2)
{
    // Set up parent variables.
    m_numPoints = r.size();
    m_xi[0] = r;
    m_xi[1] = s;

    // Construct a mapping (i,j) -> m from the triangular tensor product space
    // (i,j) to a single ordering m.
    for (int i = 0; i <= m_degree; ++i)
    {
        for (int j = 0; j <= m_degree - i; ++j)
        {
            m_ordering.push_back(std::make_pair(i,j));
        }
    }

    // Calculate collapsed coordinates from r/s values
    m_eta[0] = Array<OneD, NekDouble>(m_numPoints);
    m_eta[1] = Array<OneD, NekDouble>(m_numPoints);

    for (int i = 0; i < m_numPoints; ++i)
    {
        if (fabs(m_xi[1][i]-1.0) < NekConstants::kNekZeroTol)
        {
            m_eta[0][i] = -1.0;
            m_eta[1][i] =  1.0;
        }
        else
        {
            m_eta[0][i] = 2*(1+m_xi[0][i])/(1-m_xi[1][i])-1.0;
            m_eta[1][i] = m_xi[1][i];
        }
    }
}

/**
 * @brief Return the value of the modal functions for the triangular element at
 * the nodal points #m_xi for a given mode.
 *
 * In a triangle, we use the orthogonal basis
 *
 * \f[
 * \psi_{m(ij)} = \sqrt{2} P^{(0,0)}_i(\xi_1) P_j^{(2i+1,0)}(\xi_2) (1-\xi_2)^i
 * \f]
 *
 * where \f$ m(ij) \f$ is the mapping defined in NodalUtilTriangle::m_ordering
 * and \f$ J_n^{(\alpha,\beta)}(z) \f$ denotes the standard Jacobi polynomial.
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 *
 * @return Vector containing orthogonal basis evaluated at the points #m_xi.
 */
NekVector<NekDouble> NodalUtilTriangle::v_OrthoBasis(const int mode)
{
    std::vector<NekDouble> jacobi_i(m_numPoints), jacobi_j(m_numPoints);
    std::pair<int, int> modes = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_eta[0][0], &jacobi_i[0], NULL, modes.first, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[1][0], &jacobi_j[0], NULL, modes.second,
        2.0 * modes.first + 1.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);
    NekDouble sqrt2 = sqrt(2.0);

    for (int i = 0; i < m_numPoints; ++i)
    {
        ret[i] = sqrt2 * jacobi_i[i] * jacobi_j[i] *
            pow(1.0 - m_eta[1][i], modes.first);
    }

    return ret;
}

/**
 * @brief Return the value of the derivative of the modal functions for the
 * triangular element at the nodal points #m_xi for a given mode.
 *
 * Note that this routine must use the chain rule combined with the collapsed
 * coordinate derivatives as described in Sherwin & Karniadakis (2nd edition),
 * pg 150.
 *
 * @param dir   Coordinate direction in which to evaluate the derivative.
 * @param mode  The mode of the orthogonal basis to evaluate.
 *
 * @return Vector containing the derivative of the orthogonal basis evaluated at
 *         the points #m_xi.
 */
NekVector<NekDouble> NodalUtilTriangle::v_OrthoBasisDeriv(
    const int dir, const int mode)
{
    std::vector<NekDouble> jacobi_i(m_numPoints), jacobi_j(m_numPoints);
    std::vector<NekDouble> jacobi_di(m_numPoints), jacobi_dj(m_numPoints);
    std::pair<int, int> modes = m_ordering[mode];

    // Calculate Jacobi polynomials and their derivatives. Note that we use both
    // jacobfd and jacobd since jacobfd is only valid for derivatives in the
    // open interval (-1,1).
    Polylib::jacobfd(
        m_numPoints, &m_eta[0][0], &jacobi_i[0], NULL, modes.first, 0.0,
        0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[1][0], &jacobi_j[0], NULL, modes.second,
        2.0*modes.first + 1.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[0][0], &jacobi_di[0], modes.first, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[1][0], &jacobi_dj[0], modes.second,
        2.0*modes.first + 1.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);
    NekDouble sqrt2 = sqrt(2.0);

    if (dir == 0)
    {
        // d/d(\xi_1) = 2/(1-\eta_2) d/d(\eta_1)
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = 2.0 * sqrt2 * jacobi_di[i] * jacobi_j[i];
            if (modes.first > 0)
            {
                ret[i] *= pow(1.0 - m_eta[1][i], modes.first - 1.0);
            }
        }
    }
    else
    {
        // d/d(\xi_2) = 2(1+\eta_1)/(1-\eta_2) d/d(\eta_1) + d/d(eta_2)
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = (1 + m_eta[0][i]) * sqrt2 * jacobi_di[i] * jacobi_j[i];
            if (modes.first > 0)
            {
                ret[i] *= pow(1.0 - m_eta[1][i], modes.first - 1.0);
            }

            NekDouble tmp = jacobi_dj[i] * pow(1.0 - m_eta[1][i], modes.first);
            if (modes.first > 0)
            {
                tmp -= modes.first * jacobi_j[i] *
                    pow(1.0 - m_eta[1][i], modes.first-1);
            }

            ret[i] += sqrt2 * jacobi_i[i] * tmp;
        }
    }

    return ret;
}

/**
 * @brief Construct the nodal utility class for a tetrahedron.
 *
 * The constructor of this class sets up two member variables used in the
 * evaluation of the orthogonal basis:
 *
 * - NodalUtilTetrahedron::m_eta is used to construct the collapsed coordinate
 *   locations of the nodal points \f$ (\eta_1, \eta_2, \eta_3) \f$ inside the
 *   cube \f$[-1,1]^3\f$ on which the orthogonal basis functions are defined.
 * - NodalUtilTetrahedron::m_ordering constructs a mapping from the index set
 *   \f$ I = \{ (i,j,k)\ |\ 0\leq i,j,k \leq P, i+j \leq P, i+j+k \leq P \}\f$
 *   to an ordering \f$ 0 \leq m(ijk) \leq (P+1)(P+2)(P+3)/6 \f$ that defines
 *   the monomials \f$ \xi_1^i \xi_2^j \xi_3^k \f$ that span the tetrahedral
 *   space. This is then used to calculate which \f$ (i,j,k) \f$ triple
 *   (represented as a tuple) corresponding to a column of the Vandermonde
 *   matrix when calculating the orthogonal polynomials.
 *
 * @param degree  Polynomial order of this nodal tetrahedron
 * @param r       \f$ \xi_1 \f$-coordinates of nodal points in the standard
 *                element.
 * @param s       \f$ \xi_2 \f$-coordinates of nodal points in the standard
 *                element.
 * @param t       \f$ \xi_3 \f$-coordinates of nodal points in the standard
 *                element.
 */
NodalUtilTetrahedron::NodalUtilTetrahedron(int                    degree,
                                           Array<OneD, NekDouble> r,
                                           Array<OneD, NekDouble> s,
                                           Array<OneD, NekDouble> t)
    : NodalUtil(degree, 3), m_eta(3)
{
    m_numPoints = r.size();
    m_xi[0] = r;
    m_xi[1] = s;
    m_xi[2] = t;

    for (int i = 0; i <= m_degree; ++i)
    {
        for (int j = 0; j <= m_degree - i; ++j)
        {
            for (int k = 0; k <= m_degree - i - j; ++k)
            {
                m_ordering.push_back(Mode(i, j, k));
            }
        }
    }

    // Calculate collapsed coordinates from r/s values
    m_eta[0] = Array<OneD, NekDouble>(m_numPoints);
    m_eta[1] = Array<OneD, NekDouble>(m_numPoints);
    m_eta[2] = Array<OneD, NekDouble>(m_numPoints);

    for (int i = 0; i < m_numPoints; ++i)
    {
        if (fabs(m_xi[2][i] - 1.0) < NekConstants::kNekZeroTol)
        {
            // Very top point of the tetrahedron
            m_eta[0][i] = -1.0;
            m_eta[1][i] = -1.0;
            m_eta[2][i] = m_xi[2][i];
        }
        else
        {
            if (fabs(m_xi[1][i] - 1.0) <  NekConstants::kNekZeroTol)
            {
                // Distant diagonal edge shared by all eta_x coordinate planes:
                // the xi_y == -xi_z line
                m_eta[0][i] = -1.0;
            }
            else if (fabs(m_xi[1][i] + m_xi[2][i]) < NekConstants::kNekZeroTol)
            {
                m_eta[0][i] = -1.0;
            }
            else
            {
                m_eta[0][i] = 2.0 * (1.0 + m_xi[0][i]) /
                    (-m_xi[1][i] - m_xi[2][i]) - 1.0;
            }
            m_eta[1][i] = 2.0 * (1.0 + m_xi[1][i]) / (1.0 - m_xi[2][i]) - 1.0;
            m_eta[2][i] = m_xi[2][i];
        }
    }
}

/**
 * @brief Return the value of the modal functions for the tetrahedral element at
 * the nodal points #m_xi for a given mode.
 *
 * In a tetrahedron, we use the orthogonal basis
 *
 * \f[ \psi_{m(ijk)} = \sqrt{8} P^{(0,0)}_i(\xi_1) P_j^{(2i+1,0)}(\xi_2)
 * P_k^{(2i+2j+2,0)}(\xi_3) (1-\xi_2)^i (1-\xi_3)^{i+j} \f]
 *
 * where \f$ m(ijk) \f$ is the mapping defined in #m_ordering and \f$
 * J_n^{(\alpha,\beta)}(z) \f$ denotes the standard Jacobi polynomial.
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 *
 * @return Vector containing orthogonal basis evaluated at the points #m_xi.
 */
NekVector<NekDouble> NodalUtilTetrahedron::v_OrthoBasis(const int mode)
{
    std::vector<NekDouble> jacA(m_numPoints), jacB(m_numPoints);
    std::vector<NekDouble> jacC(m_numPoints);

    int I, J, K;
    std::tie(I, J, K) = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_eta[0][0], &jacA[0], NULL, I, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[1][0], &jacB[0], NULL, J, 2.0 * I + 1.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[2][0], &jacC[0], NULL, K, 2.0 * (I+J) + 2.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);
    NekDouble sqrt8 = sqrt(8.0);

    for (int i = 0; i < m_numPoints; ++i)
    {
        ret[i] = sqrt8 * jacA[i] * jacB[i] * jacC[i] *
            pow(1.0 - m_eta[1][i], I) * pow(1.0 - m_eta[2][i], I + J);
    }

    return ret;
}

/**
 * @brief Return the value of the derivative of the modal functions for the
 * tetrahedral element at the nodal points #m_xi for a given mode.
 *
 * Note that this routine must use the chain rule combined with the collapsed
 * coordinate derivatives as described in Sherwin & Karniadakis (2nd edition),
 * pg 152.
 *
 * @param dir   Coordinate direction in which to evaluate the derivative.
 * @param mode  The mode of the orthogonal basis to evaluate.
 *
 * @return Vector containing the derivative of the orthogonal basis evaluated at
 *         the points #m_xi.
 */
NekVector<NekDouble> NodalUtilTetrahedron::v_OrthoBasisDeriv(
    const int dir, const int mode)
{
    std::vector<NekDouble> jacA(m_numPoints), jacB(m_numPoints);
    std::vector<NekDouble> jacC(m_numPoints);
    std::vector<NekDouble> jacDerivA(m_numPoints), jacDerivB(m_numPoints);
    std::vector<NekDouble> jacDerivC(m_numPoints);

    int I, J, K;
    std::tie(I, J, K) = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_eta[0][0], &jacA[0], NULL, I, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[1][0], &jacB[0], NULL, J, 2.0 * I + 1.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[2][0], &jacC[0], NULL, K, 2.0 * (I+J) + 2.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[0][0], &jacDerivA[0], I, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[1][0], &jacDerivB[0], J, 2.0 * I + 1.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[2][0], &jacDerivC[0], K, 2.0 * (I+J) + 2.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);
    NekDouble sqrt8 = sqrt(8.0);

    // Always compute x-derivative since this term appears in the latter two
    // terms.
    for (int i = 0; i < m_numPoints; ++i)
    {
        ret[i] = 4.0 * sqrt8 * jacDerivA[i] * jacB[i] * jacC[i];

        if (I > 0)
        {
            ret[i] *= pow(1 - m_eta[1][i], I - 1);
        }

        if (I + J > 0)
        {
            ret[i] *= pow(1 - m_eta[2][i], I + J - 1);
        }
    }

    if (dir >= 1)
    {
        // Multiply by (1+a)/2
        NekVector<NekDouble> tmp(m_numPoints);

        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] *= 0.5 * (m_eta[0][i] + 1.0);

            tmp[i] = 2.0 * sqrt8 * jacA[i] * jacC[i];
            if (I + J > 0)
            {
                tmp[i] *= pow(1.0 - m_eta[2][i], I + J - 1);
            }

            NekDouble tmp2 = jacDerivB[i] * pow(1.0 - m_eta[1][i], I);
            if (I > 0)
            {
                tmp2 -= I * jacB[i] * pow(1.0 - m_eta[1][i], I - 1);
            }

            tmp[i] *= tmp2;
        }

        if (dir == 1)
        {
            ret += tmp;
            return ret;
        }

        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] += 0.5 * (1.0 + m_eta[1][i]) * tmp[i];

            NekDouble tmp2 = jacDerivC[i] * pow(1.0 - m_eta[2][i], I + J);
            if (I + J > 0)
            {
                tmp2 -= jacC[i] * (I + J) * pow(1.0 - m_eta[2][i], I + J - 1);
            }

            ret[i] += sqrt8 * jacA[i] * jacB[i] * pow(1.0 - m_eta[1][i], I) *
                tmp2;
        }
    }

    return ret;
}

/**
 * @brief Construct the nodal utility class for a prism.
 *
 * The constructor of this class sets up two member variables used in the
 * evaluation of the orthogonal basis:
 *
 * - NodalUtilPrism::m_eta is used to construct the collapsed coordinate
 *   locations of the nodal points \f$ (\eta_1, \eta_2, \eta_3) \f$ inside the
 *   cube \f$[-1,1]^3\f$ on which the orthogonal basis functions are defined.
 * - NodalUtilPrism::m_ordering constructs a mapping from the index set
 *   \f$ I = \{ (i,j,k)\ |\ 0\leq i,j,k \leq P, i+k \leq P \}\f$ to an ordering
 *   \f$ 0 \leq m(ijk) \leq (P+1)(P+1)(P+2)/2 \f$ that defines the monomials \f$
 *   \xi_1^i \xi_2^j \xi_3^k \f$ that span the prismatic space. This is then
 *   used to calculate which \f$ (i,j,k) \f$ triple (represented as a tuple)
 *   corresponding to a column of the Vandermonde matrix when calculating the
 *   orthogonal polynomials.
 *
 * @param degree  Polynomial order of this nodal tetrahedron
 * @param r       \f$ \xi_1 \f$-coordinates of nodal points in the standard
 *                element.
 * @param s       \f$ \xi_2 \f$-coordinates of nodal points in the standard
 *                element.
 * @param t       \f$ \xi_3 \f$-coordinates of nodal points in the standard
 *                element.
 */
NodalUtilPrism::NodalUtilPrism(int                    degree,
                               Array<OneD, NekDouble> r,
                               Array<OneD, NekDouble> s,
                               Array<OneD, NekDouble> t)
    : NodalUtil(degree, 3), m_eta(3)
{
    m_numPoints = r.size();
    m_xi[0] = r;
    m_xi[1] = s;
    m_xi[2] = t;

    for (int i = 0; i <= m_degree; ++i)
    {
        for (int j = 0; j <= m_degree; ++j)
        {
            for (int k = 0; k <= m_degree - i; ++k)
            {
                m_ordering.push_back(Mode(i, j, k));
            }
        }
    }

    // Calculate collapsed coordinates from r/s values
    m_eta[0] = Array<OneD, NekDouble>(m_numPoints);
    m_eta[1] = Array<OneD, NekDouble>(m_numPoints);
    m_eta[2] = Array<OneD, NekDouble>(m_numPoints);

    for (int i = 0; i < m_numPoints; ++i)
    {
        if (fabs(m_xi[2][i] - 1.0) < NekConstants::kNekZeroTol)
        {
            // Very top point of the prism
            m_eta[0][i] = -1.0;
            m_eta[1][i] = m_xi[1][i];
            m_eta[2][i] = 1.0;
        }
        else
        {
            // Third basis function collapsed to "pr" direction instead of "qr"
            // direction
            m_eta[0][i] = 2.0*(1.0 + m_xi[0][i])/(1.0 - m_xi[2][i]) - 1.0;
            m_eta[1][i] = m_xi[1][i];
            m_eta[2][i] = m_xi[2][i];
        }
    }
}

/**
 * @brief Return the value of the modal functions for the prismatic element at
 * the nodal points #m_xi for a given mode.
 *
 * In a prism, we use the orthogonal basis
 *
 * \f[ \psi_{m(ijk)} = \sqrt{2} P^{(0,0)}_i(\xi_1) P_j^{(0,0)}(\xi_2)
 * P_k^{(2i+1,0)}(\xi_3) (1-\xi_3)^i \f]
 *
 * where \f$ m(ijk) \f$ is the mapping defined in #m_ordering and \f$
 * J_n^{(\alpha,\beta)}(z) \f$ denotes the standard Jacobi polynomial.
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 *
 * @return Vector containing orthogonal basis evaluated at the points #m_xi.
 */
NekVector<NekDouble> NodalUtilPrism::v_OrthoBasis(const int mode)
{
    std::vector<NekDouble> jacA(m_numPoints), jacB(m_numPoints);
    std::vector<NekDouble> jacC(m_numPoints);

    int I, J, K;
    std::tie(I, J, K) = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_eta[0][0], &jacA[0], NULL, I, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[1][0], &jacB[0], NULL, J, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[2][0], &jacC[0], NULL, K, 2.0 * I + 1.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);
    NekDouble sqrt2 = sqrt(2.0);

    for (int i = 0; i < m_numPoints; ++i)
    {
        ret[i] = sqrt2 * jacA[i] * jacB[i] * jacC[i] *
            pow(1.0 - m_eta[2][i], I);
    }

    return ret;
}

/**
 * @brief Return the value of the derivative of the modal functions for the
 * prismatic element at the nodal points #m_xi for a given mode.
 *
 * Note that this routine must use the chain rule combined with the collapsed
 * coordinate derivatives as described in Sherwin & Karniadakis (2nd edition),
 * pg 152.
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 * @param dir   Coordinate direction in which to evaluate the derivative.
 *
 * @return Vector containing the derivative of the orthogonal basis evaluated at
 *         the points #m_xi.
 */
NekVector<NekDouble> NodalUtilPrism::v_OrthoBasisDeriv(
    const int dir, const int mode)
{
    std::vector<NekDouble> jacA(m_numPoints), jacB(m_numPoints);
    std::vector<NekDouble> jacC(m_numPoints);
    std::vector<NekDouble> jacDerivA(m_numPoints), jacDerivB(m_numPoints);
    std::vector<NekDouble> jacDerivC(m_numPoints);

    int I, J, K;
    std::tie(I, J, K) = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_eta[0][0], &jacA[0], NULL, I, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[1][0], &jacB[0], NULL, J, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_eta[2][0], &jacC[0], NULL, K, 2.0 * I + 1.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[0][0], &jacDerivA[0], I, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[1][0], &jacDerivB[0], J, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_eta[2][0], &jacDerivC[0], K, 2.0 * I + 1.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);
    NekDouble sqrt2 = sqrt(2.0);

    if (dir == 1)
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = sqrt2 * jacA[i] * jacDerivB[i] * jacC[i] *
                pow(1.0 - m_eta[2][i], I);
        }
    }
    else
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = 2.0 * sqrt2 * jacDerivA[i] * jacB[i] * jacC[i];

            if (I > 0)
            {
                ret[i] *= pow(1.0 - m_eta[2][i], I - 1);
            }
        }

        if (dir == 0)
        {
            return ret;
        }

        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] *= 0.5 * (1.0 + m_eta[0][i]);

            NekDouble tmp = jacDerivC[i] * pow(1.0 - m_eta[2][i], I);

            if (I > 0)
            {
                tmp -= jacC[i] * I * pow(1.0 - m_eta[2][i], I - 1);
            }

            ret[i] += sqrt2 * jacA[i] * jacB[i] * tmp;
        }
    }

    return ret;
}

/**
 * @brief Construct the nodal utility class for a quadrilateral.
 *
 * The constructor of this class sets up the #m_ordering member variable used in
 * the evaluation of the orthogonal basis.
 *
 * @param degree  Polynomial order of this nodal quad.
 * @param r       \f$ \xi_1 \f$-coordinates of nodal points in the standard
 *                element.
 * @param s       \f$ \xi_2 \f$-coordinates of nodal points in the standard
 *                element.
 */
NodalUtilQuad::NodalUtilQuad(int                    degree,
                             Array<OneD, NekDouble> r,
                             Array<OneD, NekDouble> s)
    : NodalUtil(degree, 2)
{
    // Set up parent variables.
    m_numPoints = r.size();
    m_xi[0] = r;
    m_xi[1] = s;

    // Construct a mapping (i,j) -> m from the tensor product space (i,j) to a
    // single ordering m.
    for (int j = 0; j <= m_degree; ++j)
    {
        for (int i = 0; i <= m_degree; ++i)
        {
            m_ordering.push_back(std::make_pair(i,j));
        }
    }
}

/**
 * @brief Return the value of the modal functions for the quad element at
 * the nodal points #m_xi for a given mode.
 *
 * In a quad, we use the orthogonal basis
 *
 * \f[
 * \psi_{m(ij)} = P^{(0,0)}_i(\xi_1) P_j^{(0,0)}(\xi_2)
 * \f]
 *
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 */
NekVector<NekDouble> NodalUtilQuad::v_OrthoBasis(const int mode)
{
    std::vector<NekDouble> jacobi_i(m_numPoints), jacobi_j(m_numPoints);
    std::pair<int, int> modes = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_xi[0][0], &jacobi_i[0], NULL, modes.first, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_xi[1][0], &jacobi_j[0], NULL, modes.second, 0.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);

    for (int i = 0; i < m_numPoints; ++i)
    {
        ret[i] = jacobi_i[i] * jacobi_j[i];
    }

    return ret;
}

/**
 * @brief Return the value of the derivative of the modal functions for the
 * quadrilateral element at the nodal points #m_xi for a given mode.
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 */
NekVector<NekDouble> NodalUtilQuad::v_OrthoBasisDeriv(
    const int dir, const int mode)
{
    std::vector<NekDouble> jacobi_i(m_numPoints), jacobi_j(m_numPoints);
    std::vector<NekDouble> jacobi_di(m_numPoints), jacobi_dj(m_numPoints);
    std::pair<int, int> modes = m_ordering[mode];

    // Calculate Jacobi polynomials and their derivatives. Note that we use both
    // jacobfd and jacobd since jacobfd is only valid for derivatives in the
    // open interval (-1,1).
    Polylib::jacobfd(
        m_numPoints, &m_xi[0][0], &jacobi_i[0], NULL, modes.first, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_xi[1][0], &jacobi_j[0], NULL, modes.second, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_xi[0][0], &jacobi_di[0], modes.first, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_xi[1][0], &jacobi_dj[0], modes.second, 0.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);

    if (dir == 0)
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = jacobi_di[i] * jacobi_j[i];
        }
    }
    else
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = jacobi_i[i] * jacobi_dj[i];
        }
    }

    return ret;
}


/**
 * @brief Construct the nodal utility class for a hexahedron.
 *
 * The constructor of this class sets up the #m_ordering member variable used in
 * the evaluation of the orthogonal basis.
 *
 * @param degree  Polynomial order of this nodal hexahedron.
 * @param r       \f$ \xi_1 \f$-coordinates of nodal points in the standard
 *                element.
 * @param s       \f$ \xi_2 \f$-coordinates of nodal points in the standard
 *                element.
 */
NodalUtilHex::NodalUtilHex(int degree,
                           Array<OneD, NekDouble> r,
                           Array<OneD, NekDouble> s,
                           Array<OneD, NekDouble> t)
    : NodalUtil(degree, 3)
{
    // Set up parent variables.
    m_numPoints = r.size();
    m_xi[0] = r;
    m_xi[1] = s;
    m_xi[2] = t;

    // Construct a mapping (i,j,k) -> m from the tensor product space (i,j,k) to
    // a single ordering m.
    for (int k = 0; k <= m_degree; ++k)
    {
        for (int j = 0; j <= m_degree; ++j)
        {
            for (int i = 0; i <= m_degree; ++i)
            {
                m_ordering.push_back(Mode(i, j, k));
            }
        }
    }
}

/**
 * @brief Return the value of the modal functions for the hex element at
 * the nodal points #m_xi for a given mode.
 *
 * In a quad, we use the orthogonal basis
 *
 * \f[
 * \psi_{m(ijk)} = P^{(0,0)}_i(\xi_1) P_j^{(0,0)}(\xi_2) P_k^{(0,0)}(\xi_3)
 * \f]
 *
 * @param mode  The mode of the orthogonal basis to evaluate.
 */
NekVector<NekDouble> NodalUtilHex::v_OrthoBasis(const int mode)
{
    std::vector<NekDouble> jacobi_i(m_numPoints), jacobi_j(m_numPoints);
    std::vector<NekDouble> jacobi_k(m_numPoints);

    int I, J, K;
    std::tie(I, J, K) = m_ordering[mode];

    // Calculate Jacobi polynomials
    Polylib::jacobfd(
        m_numPoints, &m_xi[0][0], &jacobi_i[0], NULL, I, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_xi[1][0], &jacobi_j[0], NULL, J, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_xi[2][0], &jacobi_k[0], NULL, K, 0.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);

    for (int i = 0; i < m_numPoints; ++i)
    {
        ret[i] = jacobi_i[i] * jacobi_j[i] * jacobi_k[i];
    }

    return ret;
}

NekVector<NekDouble> NodalUtilHex::v_OrthoBasisDeriv(
    const int dir, const int mode)
{
    std::vector<NekDouble> jacobi_i(m_numPoints), jacobi_j(m_numPoints);
    std::vector<NekDouble> jacobi_k(m_numPoints);
    std::vector<NekDouble> jacobi_di(m_numPoints), jacobi_dj(m_numPoints);
    std::vector<NekDouble> jacobi_dk(m_numPoints);

    int I, J, K;
    std::tie(I, J, K) = m_ordering[mode];

    // Calculate Jacobi polynomials and their derivatives. Note that we use both
    // jacobfd and jacobd since jacobfd is only valid for derivatives in the
    // open interval (-1,1).
    Polylib::jacobfd(
        m_numPoints, &m_xi[0][0], &jacobi_i[0], NULL, I, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_xi[1][0], &jacobi_j[0], NULL, J, 0.0, 0.0);
    Polylib::jacobfd(
        m_numPoints, &m_xi[2][0], &jacobi_k[0], NULL, K, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_xi[0][0], &jacobi_di[0], I, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_xi[1][0], &jacobi_dj[0], J, 0.0, 0.0);
    Polylib::jacobd(
        m_numPoints, &m_xi[2][0], &jacobi_dk[0], K, 0.0, 0.0);

    NekVector<NekDouble> ret(m_numPoints);

    if (dir == 0)
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = jacobi_di[i] * jacobi_j[i] * jacobi_k[i];
        }
    }
    else if (dir == 1)
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = jacobi_dj[i] * jacobi_i[i] * jacobi_k[i];
        }
    }
    else
    {
        for (int i = 0; i < m_numPoints; ++i)
        {
            ret[i] = jacobi_i[i] * jacobi_j[i] * jacobi_dk[i];
        }
    }

    return ret;
}

}
}

