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

            PointsKeyArray p(m_expDim);
            p[0] = m_coords[0]->GetBasis(0)->GetPointsKey();
            p[1] = m_coords[0]->GetBasis(1)->GetPointsKey();
            int nqtot = p[0].GetNumPoints() *
                        p[1].GetNumPoints();
            int pts = (m_type == eRegular || m_type == eMovingRegular)
                            ? 1 : nqtot;

            DerivStorage deriv = GetDeriv(p);

            Array<OneD, NekDouble> jac(pts, 0.0);
            Vmath::Vvtvvtm(pts, &deriv[0][0][0], 1, &deriv[1][1][0], 1,
                                &deriv[1][0][0], 1, &deriv[0][1][0], 1,
                                &jac[0],           1);

            if(Vmath::Vmin(pts, &jac[0], 1) < 0)
            {
                m_valid = false;
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
