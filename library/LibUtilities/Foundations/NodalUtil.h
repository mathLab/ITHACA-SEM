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
// Description: 2D and 3D Nodal Triangle and Tetrahedron Utilities header file
//              Basis function, Interpolation, Integral, Derivation, etc.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALUTIL_H
#define NODALUTIL_H

#include <boost/tuple/tuple.hpp>

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/Foundations/Points.h>


//#include <LibUtilities/BasicUtils/BasicUtilsFwd.hpp>  // for NekManager
#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{
namespace LibUtilities
{

typedef boost::shared_ptr<NekMatrix<NekDouble> > SharedMatrix;

/**
 * @brief A class to assist in the construction of nodal simplex and hybrid
 * elements in two and three dimensions.
 *
 * The NodalUtil class and its subclasses are designed to take care of some
 * common issues that arise when considering triangles, tetrahedra and prismatic
 * elements that are equipped with a nodal Lagrangian basis, defined using a set
 * of nodal points that we store in the array NodalUtil::m_xi. Since one cannot
 * write this basis analytically, we instead construct the Vandermonde matrix
 *
 * \f[ \mathbf{V}_{ij} \f]
 */
class NodalUtil
{
public:
    NodalUtil(int degree, int dim) : m_dim(dim), m_degree(degree), m_xi(dim)
    {
    }

    LIB_UTILITIES_EXPORT NekVector<NekDouble> GetWeights();
    LIB_UTILITIES_EXPORT SharedMatrix GetVandermonde();
    LIB_UTILITIES_EXPORT SharedMatrix GetVandermondeForDeriv(int dir);
    LIB_UTILITIES_EXPORT SharedMatrix GetDerivMatrix(int dir);
    LIB_UTILITIES_EXPORT SharedMatrix GetInterpolationMatrix(
        Array<OneD, Array<OneD, NekDouble> > &xi);

protected:
    int m_dim;
    int m_degree;
    int m_numPoints;
    Array<OneD, Array<OneD, NekDouble> > m_xi;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode) = 0;
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode) = 0;
    virtual boost::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi) = 0;
    virtual NekDouble v_ModeZeroIntegral() = 0;
    virtual int v_NumModes() = 0;
};

class NodalUtilTriangle : public NodalUtil
{
public:
    NodalUtilTriangle(int degree,
                      Array<OneD, NekDouble> r,
                      Array<OneD, NekDouble> s);

protected:
    std::vector<std::pair<int, int> > m_ordering;
    Array<OneD, Array<OneD, NekDouble> > m_eta;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual boost::shared_ptr<NodalUtil> v_CreateUtil(
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

class NodalUtilTetrahedron : public NodalUtil
{
    typedef boost::tuple<int, int, int> Mode;

public:
    NodalUtilTetrahedron(int degree,
                         Array<OneD, NekDouble> r,
                         Array<OneD, NekDouble> s,
                         Array<OneD, NekDouble> t);

protected:
    std::vector<Mode> m_ordering;
    Array<OneD, Array<OneD, NekDouble> > m_eta;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual boost::shared_ptr<NodalUtil> v_CreateUtil(
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

class NodalUtilPrism : public NodalUtil
{
    typedef boost::tuple<int, int, int> Mode;

public:
    NodalUtilPrism(int degree,
                   Array<OneD, NekDouble> r,
                   Array<OneD, NekDouble> s,
                   Array<OneD, NekDouble> t);

protected:
    std::vector<Mode> m_ordering;
    Array<OneD, Array<OneD, NekDouble> > m_eta;

    virtual NekVector<NekDouble> v_OrthoBasis(const int mode);
    virtual NekVector<NekDouble> v_OrthoBasisDeriv(
        const int dir, const int mode);

    virtual boost::shared_ptr<NodalUtil> v_CreateUtil(
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

} // end of LibUtilities namespace
} // end of Nektar namespace

#endif //NODALUTIL_H
