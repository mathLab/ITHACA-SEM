///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors2D.cpp
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
// Description: Implementation of 2D geometric factors.
//
///////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactors2D.h>
#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * @class GeomFactors2D
         *
         * Computes and stores geometric factors and pointwise geometry
         * information for 2D expansions.
         *
         * In the case that a three-dimensional coordinate system is
         * indicated, the Jacobian is not square so cannot be directly inverted.
         *
         *
         * @see GeomFactors
         */

        /**
         * Initialises the class and populates the geometric factors. The
         * physical coordinates of the mapping are computed, along with the
         * derivative of the mapping. This is interpolated onto the target
         * expansion basis and passed to SetUpJacGmat2D for the remaining terms
         * to be computed.
         *
         * @param   gtype       Type of geometry (e.g. regular or deformed).
         * @param   coordim     Dimension of coordinate system.
         * @param   Coords      Coordinate information.
         * @param   tbasis      Target basis for geometric information.
         */
        GeomFactors2D::GeomFactors2D(
                        const GeomType gtype,
                        const int coordim,
                        const Array<OneD, const StdRegions
                                            ::StdExpansion2DSharedPtr> &Coords,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                             &tbasis) :
            GeomFactors(gtype, 2, coordim)
        {
            // Sanity checks.
            ASSERTL1(coordim == 2 || coordim == 3,
                     "The coordinate dimension should be equal to two or three"
                     "for two-dimensional elements");
            ASSERTL1(tbasis.num_elements() == 2,
                     "tbasis should be an array of size two");

            // Copy shared pointers to coordinates.
            for (int i = 0; i < m_coordDim; ++i)
            {
                m_coords[i] = Coords[i];
            }

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            m_pointsKey[0] = tbasis[0]->GetPointsKey();
            m_pointsKey[1] = tbasis[1]->GetPointsKey();

            // Calculate the Jacobian and metric terms.
            SetUpJacGmat2D();

            // Test if the element is valid.
            CheckIfValid();
        }


        /**
         *
         */
        GeomFactors2D::GeomFactors2D(const GeomFactors2D& S) :
            GeomFactors(S)
        {
        }


        /**
         *
         */
        GeomFactors2D::~GeomFactors2D()
        {
        }


        /**
         *
         */
        void GeomFactors2D::SetUpJacGmat2D()
        {
            DerivStorage deriv = GetDeriv(m_pointsKey);

            // Check the number of derivative dimensions matches the coordinate
            // space.
            ASSERTL1(deriv[0].num_elements()==m_coordDim,
                     "The dimension of array d1 does not match the coordinate "
                     "dimension");
            ASSERTL1(deriv[1].num_elements()==m_coordDim,
                     "The dimension of array d2 does not match the coordinate "
                     "dimension");

            // Compute total number of points in target basis
            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();
            int pts = 1;
            int i, j, k, l;

            // Check each derivative combination has the correct sized storage
            // for all quadrature points.
            for (i = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_coordDim; ++j)
                {
                    ASSERTL1(deriv[i][j].num_elements() == nqtot,
                             "Number of quadrature points do not match");
                }
            }

            // Only compute single values if the element is regular.
            pts = (m_type == eRegular || m_type == eMovingRegular) ? 1 : nqtot;

            // Allocate storage for Jacobian, metric terms and derivative
            // factors.
            m_jac  = Array<OneD, NekDouble>(pts,                      0.0);
            m_gmat = Array<TwoD, NekDouble>(m_expDim*m_expDim,   pts, 0.0);
            m_derivFactors =
                     Array<TwoD, NekDouble>(m_expDim*m_coordDim, pts, 0.0);

            // Allocate temporary array for storing the g_{ij} terms
            Array<TwoD, NekDouble> tmp(m_expDim*m_expDim, pts, 0.0);

            // Compute g_{ij} as t_i \cdot t_j and store in tmp.
            for (i = 0, l = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j, ++l)
                {
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        Vmath::Vvtvp(pts, &deriv[i][k][0], 1,
                                          &deriv[j][k][0], 1,
                                          &tmp[l][0],        1,
                                          &tmp[l][0],        1);
                    }
                }
            }

            // Compute g = det(g_{ij}) (= Jacobian squared) and store
            // temporarily in m_jac.
            Vmath::Vvtvvtm(pts, &tmp[0][0], 1, &tmp[3][0], 1,
                                &tmp[1][0], 1, &tmp[2][0], 1, &m_jac[0], 1);

            // Compute inverse metric terms g^{ij}, which is easily calculated
            // for the 2x2 matrix g_ij.
            Vmath::Vdiv(pts, &tmp[3][0], 1, &m_jac[0],     1, &m_gmat[0][0], 1);
            Vmath::Vdiv(pts, &tmp[1][0], 1, &m_jac[0],     1, &m_gmat[1][0], 1);
            Vmath::Vdiv(pts, &tmp[2][0], 1, &m_jac[0],     1, &m_gmat[2][0], 1);
            Vmath::Vdiv(pts, &tmp[0][0], 1, &m_jac[0],     1, &m_gmat[3][0], 1);
            Vmath::Smul(pts, -1.0,          &m_gmat[1][0], 1, &m_gmat[1][0], 1);
            Vmath::Smul(pts, -1.0,          &m_gmat[2][0], 1, &m_gmat[2][0], 1);

            // Compute the Jacobian = sqrt(g)
            Vmath::Vsqrt(pts, &m_jac[0], 1, &m_jac[0], 1);

            // Compute the derivative factors
            for (k = 0, l = 0; k < m_coordDim; ++k)
            {
                for (j = 0; j < m_expDim; ++j, ++l)
                {
                    for (i = 0; i < m_expDim; ++i)
                    {
                        Vmath::Vvtvp(pts, &deriv[i][k][0],        1,
                                          &m_gmat[m_expDim*i+j][0], 1,
                                          &m_derivFactors[l][0],    1,
                                          &m_derivFactors[l][0],    1);
                    }
                }
            }
        }


        /**
         * Check if the element is valid. This check only applies to elements
         * in a 2D coordinate space. The Jacobian determinant is computed
         * directly, rather than the square-root of the metric tensor
         * determinant.
         */
        void GeomFactors2D::CheckIfValid()
        {
            // Jacobian test only makes sense in 2D coordinates
            if (GetCoordim() != 2)
            {
                return;
            }

            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();
            int pts = (m_type == eRegular || m_type == eMovingRegular)
                            ? 1 : nqtot;

            DerivStorage deriv = GetDeriv(m_pointsKey);

            Array<OneD, NekDouble> jac(pts, 0.0);
            Vmath::Vvtvvtm(pts, &deriv[0][0][0], 1, &deriv[1][1][0], 1,
                                &deriv[1][0][0], 1, &deriv[0][1][0], 1,
                                &jac[0],           1);

            if(Vmath::Vmin(pts, &jac[0], 1) < 0)
            {
                m_valid = false;
            }
        }


        /**
         * Computes the normal at each point in a 2D element, given two
         * tangential vectors, evaluated at each point.
         * @param   tbasis1     Tangential vector evaluated at each point in
         *                      the element.
         * @param   tbasis2     Second tangential vector evaluated at each
         *                      point in the element.
         */
        void GeomFactors2D::v_ComputeSurfaceNormals()
        {
            DerivStorage deriv = GetDeriv(m_pointsKey);

            // Number of dimensions.
            int coordim = deriv[0].num_elements();
            // Total number of points in the element.
            int nqtot = deriv[0][0].num_elements();

            // Allocate temporary storage.
            Array<OneD, NekDouble> temp(nqtot,0.0);

            // Initialization of storage for normals
            m_normal = Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
            for (int i = 0; i < m_coordDim; ++i)
            {
                m_normal[i] = Array<OneD,NekDouble>(nqtot);
            }

            // Derive Geometric Normal vectors by cross product of two
            // tangential basis vectors.
            Vmath::Vmul (nqtot, deriv[0][2],     1,
                                deriv[1][1],     1,
                                temp,              1);
            Vmath::Vvtvm(nqtot, deriv[0][1],     1,
                                deriv[1][2],     1,
                                temp,              1,
                                m_normal[0],       1);

            Vmath::Vmul (nqtot, deriv[0][0],     1,
                                deriv[1][2],     1,
                                temp,              1);
            Vmath::Vvtvm(nqtot, deriv[0][2],     1,
                                deriv[1][0],     1,
                                temp,              1,
                                m_normal[1],       1);

            Vmath::Vmul (nqtot, deriv[0][1],     1,
                                deriv[1][0],     1,
                                temp,              1);
            Vmath::Vvtvm(nqtot, deriv[0][0],     1,
                                deriv[1][1],     1,
                                temp,              1,
                                m_normal[2],       1);

            // Reset temp array
            temp = Array<OneD, NekDouble>(nqtot,0.0);

            // Normalization of Surface Normal
            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vvtvp(nqtot, m_normal[i],   1,
                                    m_normal[i],   1,
                                    temp,          1,
                                    temp,          1);
            }
            Vmath::Vsqrt(nqtot, temp, 1, temp, 1);

            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vdiv(nqtot,  m_normal[i],   1,
                                    temp,          1,
                                    m_normal[i],   1);
            }
        }


        void GeomFactors2D::v_Interp(
                    const PointsKeyArray &map_points,
                    const Array<OneD, const NekDouble> &src,
                    const PointsKeyArray &tpoints,
                    Array<OneD, NekDouble> &tgt) const
        {
            LibUtilities::Interp2D(map_points[0], map_points[1], src,
                                   tpoints[0], tpoints[1], tgt);
        }

        void GeomFactors2D::v_Adjoint(
                    const Array<TwoD, const NekDouble>& src,
                    Array<TwoD, NekDouble>& tgt) const
        {
            int n = src[0].num_elements();

            Vmath::Vcopy(n, &src[3][0], 1, &tgt[0][0], 1);
            Vmath::Smul (n, -1.0, &src[1][0], 1, &tgt[1][0], 1);
            Vmath::Smul (n, -1.0, &src[2][0], 1, &tgt[2][0], 1);
            Vmath::Vcopy(n, &src[0][0], 1, &tgt[3][0], 1);
        }


    }
}
