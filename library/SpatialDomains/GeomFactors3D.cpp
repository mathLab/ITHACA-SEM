///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors3D.cpp
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
// Description: Implementation of 3D geometric factors.
//
///////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactors3D.h>
#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * @class GeomFactors3D
         */

        /**
         *  The argument 'tbasis' contains the information about the quadrature
         *  points
         *              on which the weighted metric terms should be specified
         * @param   gtype       Type of geometry.
         * @param   coordim     Dimension of coordinate system.
         * @param   Coords      ?
         * @param   tbasis      Basis for tangential vectors.
         */
        GeomFactors3D::GeomFactors3D(const GeomType gtype,
                          const int coordim,
                          const Array<OneD, const StdRegions
                                            ::StdExpansion3DSharedPtr> &Coords,
                          const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis) :
            GeomFactors(gtype, 3, coordim)
        {
            ASSERTL1(coordim == 3,
                     "The coordinate dimension should be to three"
                     "for three-dimensional elements");
            ASSERTL1(tbasis.num_elements() == 3,
                     "tbasis should be an array of size three");

            for (int i = 0; i < m_coordDim; ++i)
            {
                m_coords[i] = Coords[i];
            }

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            m_pointsKey[0] = tbasis[0]->GetPointsKey();
            m_pointsKey[1] = tbasis[1]->GetPointsKey();
            m_pointsKey[2] = tbasis[2]->GetPointsKey();

            // Test if the element is valid.
            CheckIfValid();
        }



        /**
         * Computes the Jacobian of the 3D map directly from the derivatives of
         * the map to determine if it is negative and thus element is invalid.
         */
        void GeomFactors3D::CheckIfValid()
        {
            PointsKeyArray p(m_expDim);
            p[0] = m_coords[0]->GetBasis(0)->GetPointsKey();
            p[1] = m_coords[0]->GetBasis(1)->GetPointsKey();
            p[2] = m_coords[0]->GetBasis(2)->GetPointsKey();
            int nqtot = p[0].GetNumPoints() *
                        p[1].GetNumPoints() *
                        p[2].GetNumPoints();
            int pts = (m_type == eRegular || m_type == eMovingRegular)
                            ? 1 : nqtot;

            Array<OneD, NekDouble> jac(pts, 0.0);
            Array<OneD, NekDouble> tmp(pts, 0.0);
            DerivStorage deriv = GetDeriv(p);

            // J3D - Spencers book page 158
            Vmath::Vvtvvtm(pts, &deriv[1][1][0], 1, &deriv[2][2][0], 1,
                                &deriv[2][1][0], 1, &deriv[1][2][0], 1,
                                &tmp[0],           1);
            Vmath::Vvtvp  (pts, &deriv[0][0][0], 1, &tmp[0],           1,
                                &jac[0],           1, &jac[0],           1);

            Vmath::Vvtvvtm(pts, &deriv[2][1][0], 1, &deriv[0][2][0], 1,
                                &deriv[0][1][0], 1, &deriv[2][2][0], 1,
                                &tmp[0],           1);
            Vmath::Vvtvp  (pts, &deriv[1][0][0], 1, &tmp[0],           1,
                                &jac[0],           1, &jac[0],           1);

            Vmath::Vvtvvtm(pts, &deriv[0][1][0], 1, &deriv[1][2][0], 1,
                                &deriv[1][1][0], 1, &deriv[0][2][0], 1,
                                &tmp[0],           1);
            Vmath::Vvtvp  (pts, &deriv[2][0][0], 1, &tmp[0],           1,
                                &jac[0],           1, &jac[0],           1);

            if (Vmath::Vmin(pts, &jac[0], 1) < 0)
            {
                m_valid = false;
            }
        }

        void GeomFactors3D::v_Interp(
                    const PointsKeyArray &map_points,
                    const Array<OneD, const NekDouble> &src,
                    const PointsKeyArray &tpoints,
                    Array<OneD, NekDouble> &tgt) const
        {
            LibUtilities::Interp3D(map_points[0], map_points[1], map_points[2],
                                   src,
                                   tpoints[0], tpoints[1], tpoints[2], tgt);
        }

        void GeomFactors3D::v_Adjoint(
                    const Array<TwoD, const NekDouble>& src,
                    Array<TwoD, NekDouble>& tgt) const
        {
            int a, b, c, d, e, i, j;
            int n = src[0].num_elements();
            // Compute g^{ij} by computing Cofactors(g_ij)^T
            for (i = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j)
                {
                    a = ((i+1)%m_expDim)*m_expDim + ((j+1)%m_expDim);
                    b = ((i+1)%m_expDim)*m_expDim + ((j+2)%m_expDim);
                    c = ((i+2)%m_expDim)*m_expDim + ((j+1)%m_expDim);
                    d = ((i+2)%m_expDim)*m_expDim + ((j+2)%m_expDim);
                    e = j*m_expDim + i;
                    Vmath::Vvtvvtm(n, &src[a][0], 1, &src[d][0], 1,
                                        &src[b][0], 1, &src[c][0], 1,
                                        &tgt[e][0], 1);
                }
            }
        }
    }
}
