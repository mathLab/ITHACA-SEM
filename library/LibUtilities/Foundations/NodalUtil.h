///////////////////////////////////////////////////////////////////////////////
//
// File NodalUtil.h
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
// Description: 2D and 3D Nodal Triangle and Tetrahedron Utilities header file
//              Basis function, Interpolation, Integral, Derivation, etc.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALUTIL_H
#define NODALUTIL_H

#include <tuple>

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/Foundations/Points.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{
namespace LibUtilities
{

typedef std::shared_ptr<NekMatrix<NekDouble> > SharedMatrix;

/**
 * @brief A class to assist in the construction of nodal simplex and hybrid
 * elements in two and three dimensions.
 *
 * The NodalUtil class and its subclasses are designed to take care of some
 * common issues that arise when considering triangles, tetrahedra and prismatic
 * elements that are equipped with a nodal Lagrangian basis, defined using a set
 * of nodal points \f$ \xi_i \f$ that we store in the array
 * NodalUtil::m_xi. Since one cannot write this basis analytically, we instead
 * construct the Vandermonde matrix
 *
 * \f[ \mathbf{V}_{ij} = \psi_j(\xi_i) \f]
 *
 * where \f$ \psi_j \f$ is a basis that spans the polynomial space of the
 * element. The Vandermonde matrix can then be used to construct the integration
 * weights, derivative and interpolation matrices. Although this can be any
 * basis, such as the monomial basis \f$ x^i y^j z^k \f$, in practice this is
 * numerically unstable at high polynomial orders. Elements are therefore
 * expected to use the 'traditional' modal orthogonal basis. See Sherwin &
 * Karniadakis or Hesthaven & Warburton for further details of this basis and
 * the surrounding numerical issues.
 *
 * This class therefore contains the generic logic needed to construct various
 * matrices, and subclasses override virtual functions that define the
 * orthogonal basis and its derivatives for a particular element type.
 */
class NodalUtil
{
public:
    LIB_UTILITIES_EXPORT virtual ~NodalUtil() = default;
    LIB_UTILITIES_EXPORT NekVector<NekDouble> GetWeights();
    LIB_UTILITIES_EXPORT SharedMatrix GetVandermonde();
    LIB_UTILITIES_EXPORT SharedMatrix GetVandermondeForDeriv(int dir);
    LIB_UTILITIES_EXPORT SharedMatrix GetDerivMatrix(int dir);
    LIB_UTILITIES_EXPORT SharedMatrix GetInterpolationMatrix(
        Array<OneD, Array<OneD, NekDouble> > &xi);

protected:
    /**
     * @brief Set up the NodalUtil object.
     *
     * @param dim     Dimension of the element.
     * @param degree  Polynomial degree of the element.
     */
    NodalUtil(int degree, int dim) : m_dim(dim), m_degree(degree), m_xi(dim)
    {
    }

    /// Dimension of the nodal element
    int m_dim;
    /// Degree of the nodal element
    int m_degree;
    /// Total number of nodal points
    int m_numPoints;
    /// Coordinates of the nodal points defining the basis
    Array<OneD, Array<OneD, NekDouble> > m_xi;

    /**
     * @brief Return the values of the orthogonal basis at the nodal points for
     * a given mode.
     *
     * @param mode  Mode number, which is between 0 and NodalUtil::v_NumModes()
     *              - 1.
     *
     * @return Orthogonal mode @p mode evaluated at the nodal points.
     */
    virtual NekVector<NekDouble> v_OrthoBasis(const int mode) = 0;

    /**
     * @brief Return the values of the derivative of the orthogonal basis at the
     * nodal points for a given mode.
     *
     * @param dir   Coordinate direction of derivative.
     * @param mode  Mode number, which is between 0 and NodalUtil::v_NumModes()
     *              - 1.
     */
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode) = 0;

    /**
     * @brief Construct a NodalUtil object of the appropriate element type for a
     * given set of points.
     *
     * This function is used inside NodalUtil::GetInterpolationMatrix so that
     * the (potentially non-square) Vandermonde matrix can be constructed to
     * create the interpolation matrix at an arbitrary set of points in the
     * domain.
     *
     * @param xi  Distribution of nodal points to create utility with.
     */
    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi) = 0;

    /**
     * @brief Return the value of the integral of the zero-th mode for this
     * element.
     *
     * Note that for the orthogonal basis under consideration, all modes
     * integrate to zero asides from the zero-th mode. This function is used in
     * NodalUtil::GetWeights to determine integration weights.
     */
    virtual NekDouble v_ModeZeroIntegral() = 0;

    /**
     * @brief Calculate the number of degrees of freedom for this element.
     */
    virtual int v_NumModes() = 0;
};

/**
 * @brief Specialisation of the NodalUtil class to support nodal triangular
 * elements.
 */
class NodalUtilTriangle : public NodalUtil
{
public:
    LIB_UTILITIES_EXPORT NodalUtilTriangle(int degree,
                                           Array<OneD, NekDouble> r,
                                           Array<OneD, NekDouble> s);

    LIB_UTILITIES_EXPORT virtual ~NodalUtilTriangle()
    {
    }

protected:
    /// Mapping from the \f$ (i,j) \f$ indexing of the basis to a continuous
    /// ordering.
    std::vector<std::pair<int, int> > m_ordering;

    /// Collapsed coordinates \f$ (\eta_1, \eta_2) \f$ of the nodal points.
    Array<OneD, Array<OneD, NekDouble> > m_eta;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi)
    {
        return MemoryManager<NodalUtilTriangle>::AllocateSharedPtr(
            m_degree, xi[0], xi[1]);
    }

    virtual NekDouble v_ModeZeroIntegral()
    {
        return 2.0 * sqrt(2.0);
    }

    virtual int v_NumModes()
    {
        return (m_degree + 1) * (m_degree + 2) / 2;
    }
};

/**
 * @brief Specialisation of the NodalUtil class to support nodal tetrahedral
 * elements.
 */
class NodalUtilTetrahedron : public NodalUtil
{
    typedef std::tuple<int, int, int> Mode;

public:
    LIB_UTILITIES_EXPORT NodalUtilTetrahedron(int degree,
                                              Array<OneD, NekDouble> r,
                                              Array<OneD, NekDouble> s,
                                              Array<OneD, NekDouble> t);

    LIB_UTILITIES_EXPORT virtual ~NodalUtilTetrahedron()
    {
    }

protected:
    /// Mapping from the \f$ (i,j,k) \f$ indexing of the basis to a continuous
    /// ordering.
    std::vector<Mode> m_ordering;

    /// Collapsed coordinates \f$ (\eta_1, \eta_2, \eta_3) \f$ of the nodal
    /// points.
    Array<OneD, Array<OneD, NekDouble> > m_eta;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi)
    {
        return MemoryManager<NodalUtilTetrahedron>::AllocateSharedPtr(
            m_degree, xi[0], xi[1], xi[2]);
    }

    virtual NekDouble v_ModeZeroIntegral()
    {
        return 8.0 * sqrt(2.0) / 3.0;
    }

    virtual int v_NumModes()
    {
        return (m_degree + 1) * (m_degree + 2) * (m_degree + 3) / 6;
    }
};

/**
 * @brief Specialisation of the NodalUtil class to support nodal prismatic
 * elements.
 */
class NodalUtilPrism : public NodalUtil
{
    typedef std::tuple<int, int, int> Mode;

public:
    LIB_UTILITIES_EXPORT NodalUtilPrism(int degree,
                                        Array<OneD, NekDouble> r,
                                        Array<OneD, NekDouble> s,
                                        Array<OneD, NekDouble> t);

    LIB_UTILITIES_EXPORT virtual ~NodalUtilPrism()
    {
    }

protected:
    /// Mapping from the \f$ (i,j) \f$ indexing of the basis to a continuous
    /// ordering.
    std::vector<Mode> m_ordering;

    /// Collapsed coordinates \f$ (\eta_1, \eta_2, \eta_3) \f$ of the nodal
    /// points.
    Array<OneD, Array<OneD, NekDouble> > m_eta;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi)
    {
        return MemoryManager<NodalUtilPrism>::AllocateSharedPtr(
            m_degree, xi[0], xi[1], xi[2]);
    }

    virtual NekDouble v_ModeZeroIntegral()
    {
        return 4.0 * sqrt(2.0);
    }

    virtual int v_NumModes()
    {
        return (m_degree + 1) * (m_degree + 1) * (m_degree + 2) / 2;
    }
};

/**
 * @brief Specialisation of the NodalUtil class to support nodal quad elements.
 */
class NodalUtilQuad : public NodalUtil
{
public:
    LIB_UTILITIES_EXPORT NodalUtilQuad(int degree,
                                       Array<OneD, NekDouble> r,
                                       Array<OneD, NekDouble> s);

    LIB_UTILITIES_EXPORT virtual ~NodalUtilQuad()
    {
    }

protected:
    /// Mapping from the \f$ (i,j) \f$ indexing of the basis to a continuous
    /// ordering.
    std::vector<std::pair<int, int> > m_ordering;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi)
    {
        return MemoryManager<NodalUtilQuad>::AllocateSharedPtr(
            m_degree, xi[0], xi[1]);
    }

    virtual NekDouble v_ModeZeroIntegral()
    {
        return 4.0;
    }

    virtual int v_NumModes()
    {
        return (m_degree + 1) * (m_degree + 1);
    }
};

/**
 * @brief Specialisation of the NodalUtil class to support nodal hex elements.
 */
class NodalUtilHex : public NodalUtil
{
    typedef std::tuple<int, int, int> Mode;

public:
    LIB_UTILITIES_EXPORT NodalUtilHex(int degree,
                                      Array<OneD, NekDouble> r,
                                      Array<OneD, NekDouble> s,
                                      Array<OneD, NekDouble> t);

    LIB_UTILITIES_EXPORT virtual ~NodalUtilHex()
    {
    }

protected:
    /// Mapping from the \f$ (i,j,k) \f$ indexing of the basis to a continuous
    /// ordering.
    std::vector<Mode> m_ordering;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi)
    {
        return MemoryManager<NodalUtilHex>::AllocateSharedPtr(
            m_degree, xi[0], xi[1], xi[2]);
    }

    virtual NekDouble v_ModeZeroIntegral()
    {
        return 8.0;
    }

    virtual int v_NumModes()
    {
        return (m_degree + 1) * (m_degree + 1) * (m_degree + 1);
    }
};


}
}

#endif //NODALUTIL_H
