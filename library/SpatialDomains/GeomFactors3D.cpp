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

            // The quadrature points of the mapping
            // (as specified in Coords)
            LibUtilities::PointsKey pkey0_map(Coords[0]->GetBasis(0)->GetPointsKey());
            LibUtilities::PointsKey pkey1_map(Coords[0]->GetBasis(1)->GetPointsKey());
            LibUtilities::PointsKey pkey2_map(Coords[0]->GetBasis(2)->GetPointsKey());
            int nquad0_map = pkey0_map.GetNumPoints();
            int nquad1_map = pkey1_map.GetNumPoints();
            int nquad2_map = pkey2_map.GetNumPoints();
            int nqtot_map  = nquad0_map*nquad1_map*nquad2_map;

            // The quadrature points at the points at which we
            // want to know the metrics (as specified in tbasis)
            LibUtilities::PointsKey pkey0_tbasis(tbasis[0]->GetPointsKey());
            LibUtilities::PointsKey pkey1_tbasis(tbasis[1]->GetPointsKey());
            LibUtilities::PointsKey pkey2_tbasis(tbasis[2]->GetPointsKey());
            int nquad0_tbasis = pkey0_tbasis.GetNumPoints();
            int nquad1_tbasis = pkey1_tbasis.GetNumPoints();
            int nquad2_tbasis = pkey2_tbasis.GetNumPoints();
            int nqtot_tbasis  = nquad0_tbasis*nquad1_tbasis*nquad2_tbasis;

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            m_pointsKey[0] = pkey0_tbasis;
            m_pointsKey[1] = pkey1_tbasis;
            m_pointsKey[2] = pkey2_tbasis;

            // setup temp storage
            Array<OneD, Array<OneD,NekDouble> > d1_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d2_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d3_map   (coordim);

            m_deriv = Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(m_expDim);
            m_deriv[0] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
            m_deriv[1] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
            m_deriv[2] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);

            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                d1_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d2_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d3_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                m_deriv[0][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                m_deriv[1][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                m_deriv[2][i] = Array<OneD,NekDouble>(nqtot_tbasis);

                // Transform from coefficient space to physical space
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),Coords[i]->UpdatePhys());
                // Take the derivative (calculated at the points as specified   in 'Coords')
                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1_map[i],d2_map[i],d3_map[i]);
              
                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics      ('tbasis')
                if( (pkey0_map == pkey0_tbasis) &&
                    (pkey1_map == pkey1_tbasis) &&
                    (pkey2_map == pkey2_tbasis) )
                {
                    m_deriv[0][i] = d1_map[i];
                    m_deriv[1][i] = d2_map[i];
                    m_deriv[2][i] = d3_map[i];
                }
                else
                {
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,             pkey2_map,    d1_map[i],
                                           pkey0_tbasis, pkey1_tbasis,          pkey2_tbasis, m_deriv[0][i]);
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,             pkey2_map,    d2_map[i],
                                           pkey0_tbasis, pkey1_tbasis,          pkey2_tbasis, m_deriv[1][i]);
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,             pkey2_map,    d3_map[i],
                                           pkey0_tbasis, pkey1_tbasis,          pkey2_tbasis, m_deriv[2][i]);
                }
            }

            // Calculate the Jacobian and metric terms
            SetUpJacGmat3D(m_deriv[0],m_deriv[1],m_deriv[2]);

            // Test if the element is valid.
            CheckIfValid();
        }


        /**
         *
         */
        void GeomFactors3D::SetUpJacGmat3D(
                        const Array<OneD, Array<OneD, NekDouble> > d1,
                        const Array<OneD, Array<OneD, NekDouble> > d2,
                        const Array<OneD, Array<OneD, NekDouble> > d3)
        {
            ASSERTL1(d1.num_elements()==m_coordDim,
                     "The dimension of array d1 does not match the coordinate "
                     "dimension");
            ASSERTL1(d2.num_elements()==m_coordDim,
                     "The dimension of array d2 does not match the coordinate "
                     "dimension");
            ASSERTL1(d3.num_elements()==m_coordDim,
                     "The dimension of array d3 does not match the coordinate "
                     "dimension");

            // Compute total number of points in target basis
            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints() *
                        m_pointsKey[2].GetNumPoints();

            ASSERTL1(d1[0].num_elements() == nqtot,
                     "Number of quadrature points do not match");
            ASSERTL1(d2[0].num_elements() == nqtot,
                     "Number of quadrature points do not match");
            ASSERTL1(d3[0].num_elements() == nqtot,
                     "Number of quadrature points do not match");

            int pts = 1;
            int i, j, k, l;

            // Check each derivative combination has the correct sized storage
            // for all quadrature points.
            for (i = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_coordDim; ++j)
                {
                    ASSERTL1(m_deriv[i][j].num_elements() == nqtot,
                             "Number of quadrature points do not match");
                }
            }

            DerivStorage deriv;
            FillDeriv(deriv, m_pointsKey);

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

            // Compute g_{ij} as t_i \cdot t_j and store in tmp
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

            // Compute g^{ij} by computing Cofactors(g_ij)^T
            for (i = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j)
                {
                    int a = ((i+1)%m_expDim)*m_expDim + ((j+1)%m_expDim);
                    int b = ((i+1)%m_expDim)*m_expDim + ((j+2)%m_expDim);
                    int c = ((i+2)%m_expDim)*m_expDim + ((j+1)%m_expDim);
                    int d = ((i+2)%m_expDim)*m_expDim + ((j+2)%m_expDim);
                    int e = j*m_expDim + i;
                    Vmath::Vvtvvtm(pts, &tmp[a][0], 1, &tmp[d][0], 1,
                                        &tmp[b][0], 1, &tmp[c][0], 1,
                                        &m_gmat[e][0], 1);
                }
            }

            // Compute g = det(g_{ij}) (= Jacobian squared) and store
            // temporarily in m_jac.
            Vmath::Vvtvvtp(pts, &tmp[0][0], 1, &m_gmat[0][0],          1,
                                &tmp[1][0], 1, &m_gmat[m_expDim][0],   1,
                                &m_jac[0],  1);
            Vmath::Vvtvp  (pts, &tmp[2][0], 1, &m_gmat[2*m_expDim][0], 1,
                                &m_jac[0],  1, &m_jac[0],              1);

            for (i = 0; i < m_expDim*m_expDim; ++i)
            {
                Vmath::Vdiv(pts, &m_gmat[i][0], 1, &m_jac[0], 1, &m_gmat[i][0], 1);
            }

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
         * Computes the Jacobian of the 3D map directly from the derivatives of
         * the map to determine if it is negative and thus element is invalid.
         */
        void GeomFactors3D::CheckIfValid()
        {
            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints() *
                        m_pointsKey[2].GetNumPoints();
            int pts = (m_type == eRegular || m_type == eMovingRegular)
                            ? 1 : nqtot;

            Array<OneD, NekDouble> jac(pts, 0.0);
            Array<OneD, NekDouble> tmp(pts, 0.0);

            DerivStorage deriv;
            FillDeriv(deriv, m_pointsKey);

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
                    const Array<OneD, const LibUtilities::PointsKey> &map_points,
                    const DerivStorage &src,
                    const Array<OneD, const LibUtilities::PointsKey> &tpoints,
                    DerivStorage &tgt) const
        {
            for (int i = 0; i < m_coordDim; ++i)
            {
                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics
                //   ('tbasis')
                if( (map_points[0] == tpoints[0]) &&
                    (map_points[1] == tpoints[1]) &&
                    (map_points[2] == tpoints[2]))
                {
                    tgt[0][i] = src[0][i];
                    tgt[1][i] = src[1][i];
                    tgt[2][i] = src[2][i];
                }
                else
                {
                    LibUtilities::Interp3D(map_points[0],    map_points[1],             map_points[2],    src[0][i],
                                           tpoints[0], tpoints[1],          tpoints[2], tgt[0][i]);
                    LibUtilities::Interp3D(map_points[0],    map_points[1],             map_points[2],    src[1][i],
                                           tpoints[0], tpoints[1],          tpoints[2], tgt[1][i]);
                    LibUtilities::Interp3D(map_points[0],    map_points[1],             map_points[2],    src[2][i],
                            tpoints[0], tpoints[1],          tpoints[2], tgt[2][i]);
                }
            }
        }
    }
}
