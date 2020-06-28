////////////////////////////////////////////////////////////////////////////////
//
//  File: GeomFactors.cpp
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
//  Description: Geometric factors base class.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactors.h>
#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * @class GeomFactors
         *
         * This class stores the various geometric factors associated with a
         * specific element, necessary for fundamental integration and
         * differentiation operations as well as edge and surface normals.
         *
         * Initially, these algorithms are provided with a mapping from the
         * reference region element to the physical element. Practically, this
         * is represented using a corresponding reference region element for
         * each coordinate component. Note that for straight-sided elements,
         * these elements will be of linear order. Curved elements are
         * represented using higher-order coordinate mappings. This geometric
         * order is in contrast to the order of the spectral/hp expansion order
         * on the element.
         *
         * For application of the chain rule during differentiation we require
         * the partial derivatives \f[\frac{\partial \xi_i}{\partial \chi_j}\f]
         * evaluated at the physical points of the expansion basis. We also
         * construct the inverse metric tensor \f$g^{ij}\f$ which, in the case
         * of a domain embedded in a higher-dimensional space, supports the
         * conversion of covariant quantities to contravariant quantities.
         * When the expansion dimension is equal to the coordinate dimension the
         * Jacobian of the mapping \f$\chi_j\f$ is a square matrix and
         * consequently the required terms are the entries of the inverse of the
         * Jacobian. However, in general this is not the case, so we therefore
         * implement the construction of these terms following the derivation
         * in Cantwell, et. al. \cite CaYaKiPeSh13. Given the coordinate maps
         * \f$\chi_i\f$, this comprises of five steps
         * -# Compute the terms of the Jacobian
         *    \f$\frac{\partial \chi_i}{\partial \xi_j}\f$.
         * -# Compute the metric tensor
         *    \f$g_{ij}=\mathbf{t}_{(i)}\cdot\mathbf{t}_{(j)}\f$.
         * -# Compute the square of the Jacobian determinant
         *    \f$g=|\mathbf{g}|\f$.
         * -# Compute the inverse metric tensor \f$g^{ij}\f$.
         * -# Compute the terms \f$\frac{\partial \xi_i}{\partial \chi_j}\f$.
         */

        /**
         * @param   gtype       Specified whether the geometry is regular or
         *                      deformed.
         * @param   coordim     Specifies the dimension of the coordinate
         *                      system.
         * @param   Coords      Coordinate maps of the element.
         */
        GeomFactors::GeomFactors(
                const GeomType                              gtype,
                const int                                   coordim,
                const StdRegions::StdExpansionSharedPtr    &xmap,
                const Array<OneD, Array<OneD, NekDouble> > &coords) :
            m_type(gtype),
            m_expDim(xmap->GetShapeDimension()),
            m_coordDim(coordim),
            m_valid(true),
            m_xmap(xmap),
            m_coords(coords)
        {
            CheckIfValid();
        }


        /**
         * @param   S           An instance of a GeomFactors class from which
         *                      to construct a new instance.
         */
        GeomFactors::GeomFactors(const GeomFactors &S) :
            m_type(S.m_type),
            m_expDim(S.m_expDim),
            m_coordDim(S.m_coordDim),
            m_valid(S.m_valid),
            m_xmap(S.m_xmap),
            m_coords(S.m_coords)
        {
        }


        /**
         *
         */
        GeomFactors::~GeomFactors()
        {
        }


        /**
         * Member data equivalence is tested in the following order: shape type,
         * expansion dimension, coordinate dimension and coordinates.
         */
        bool operator==(const GeomFactors &lhs, const GeomFactors &rhs)
        {
            if(!(lhs.m_type == rhs.m_type))
            {
                return false;
            }

            if(!(lhs.m_expDim == rhs.m_expDim))
            {
                return false;
            }

            if(!(lhs.m_coordDim == rhs.m_coordDim))
            {
                return false;
            }

            const Array<OneD, const NekDouble> jac_lhs =
                            lhs.ComputeJac(lhs.m_xmap->GetPointsKeys());
            const Array<OneD, const NekDouble> jac_rhs =
                            rhs.ComputeJac(rhs.m_xmap->GetPointsKeys());
            if(!(jac_lhs == jac_rhs))
            {
                return false;
            }

            return true;
        }


        /**
         * Derivatives are computed at the geometry point distributions and
         * interpolated to the target point distributions.
         *
         * @param       tpoints     Target point distributions.
         * @returns                 Derivative of coordinate map evaluated at
         *                          target point distributions.
         */
        DerivStorage GeomFactors::ComputeDeriv(
                const LibUtilities::PointsKeyVector &keyTgt) const
        {
            ASSERTL1(keyTgt.size() == m_expDim,
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");

            int i = 0, j = 0;
            int nqtot_map      = 1;
            int nqtot_tbasis   = 1;
            DerivStorage deriv = DerivStorage(m_expDim);
            DerivStorage d_map = DerivStorage(m_expDim);
            LibUtilities::PointsKeyVector map_points(m_expDim);

            // Allocate storage and compute number of points
            for (i = 0; i < m_expDim; ++i)
            {
                map_points[i]  = m_xmap->GetBasis(i)->GetPointsKey();
                nqtot_map     *= map_points[i].GetNumPoints();
                nqtot_tbasis  *= keyTgt[i].GetNumPoints();
                deriv[i] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
                d_map[i] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
            }

            // Calculate local derivatives
            for(i = 0; i < m_coordDim; ++i)
            {
                Array<OneD, NekDouble> tmp(nqtot_map);
                // Transform from coefficient space to physical space
                m_xmap->BwdTrans(m_coords[i], tmp);
                
                // Allocate storage and take the derivative (calculated at the
                // points as specified in 'Coords')
                for (j = 0; j < m_expDim; ++j)
                {
                    d_map[j][i] = Array<OneD,NekDouble>(nqtot_map);
                    deriv[j][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                    m_xmap->StdPhysDeriv(j, tmp, d_map[j][i]);
                }
            }

            for (i = 0; i < m_coordDim; ++i)
            {
                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points at which we want to know the metrics
                //   ('tbasis')
                bool same = true;
                for (j = 0; j < m_expDim; ++j)
                {
                    same = same && (map_points[j] == keyTgt[j]);
                }
                if( same )
                {
                    for (j = 0; j < m_expDim; ++j)
                    {
                        deriv[j][i] = d_map[j][i];
                    }
                }
                else
                {
                    for (j = 0; j < m_expDim; ++j)
                    {
                        Interp(map_points, d_map[j][i], keyTgt, deriv[j][i]);
                    }
                }
            }

            return deriv;
        }


        /**
         * This routine returns an array of values specifying the Jacobian
         * of the mapping at quadrature points in the element. The array
         * is either of size 1 in the case of elements having #GeomType
         * #eRegular, or of size equal to the number of quadrature points for
         * #eDeformed elements.
         *
         * @returns             Array containing the Jacobian of the coordinate
         *                      mapping at the quadrature points of the element.
         * @see                 GeomType
         */
        Array<OneD, NekDouble> GeomFactors::ComputeJac(
                const LibUtilities::PointsKeyVector &keyTgt) const
        {
            ASSERTL1(keyTgt.size() == m_expDim,
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");

            int i = 0, j = 0, k = 0, l = 0;
            int ptsTgt   = 1;

            if (m_type == eDeformed)
            {
                // Allocate storage and compute number of points
                for (i = 0; i < m_expDim; ++i)
                {
                    ptsTgt   *= keyTgt[i].GetNumPoints();
                }
            }

            // Get derivative at geometry points
            DerivStorage deriv = ComputeDeriv(keyTgt);

            Array<TwoD, NekDouble> tmp (m_expDim*m_expDim, ptsTgt, 0.0);
            Array<TwoD, NekDouble> gmat(m_expDim*m_expDim, ptsTgt, 0.0);
            Array<OneD, NekDouble> jac (ptsTgt, 0.0);

            // Compute g_{ij} as t_i \cdot t_j and store in tmp
            for (i = 0, l = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j, ++l)
                {
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        Vmath::Vvtvp(ptsTgt, &deriv[i][k][0], 1,
                                             &deriv[j][k][0], 1,
                                             &tmp[l][0],      1,
                                             &tmp[l][0],      1);
                    }
                }
            }

            Adjoint(tmp, gmat);

            // Compute g = det(g_{ij}) (= Jacobian squared) and store
            // temporarily in m_jac.
            for (i = 0; i < m_expDim; ++i)
            {
                Vmath::Vvtvp(ptsTgt, &tmp[i][0], 1, &gmat[i*m_expDim][0], 1,
                                     &jac[0], 1, &jac[0], 1);
            }

            // Compute the Jacobian = sqrt(g)
            Vmath::Vsqrt(ptsTgt, &jac[0], 1, &jac[0], 1);

            return jac;
        }


        /**
         * This routine returns a two-dimensional array of values specifying
         * the inverse metric terms associated with the coordinate mapping of
         * the corresponding reference region to the physical element. These
         * terms correspond to the \f$g^{ij}\f$ terms in \cite CaYaKiPeSh13 and,
         * in the case of an embedded manifold, map covariant quantities to
         * contravariant quantities. The leading index of the array is the index
         * of the term in the tensor numbered as
         * \f[\left(\begin{array}{ccc}
         *    0 & 1 & 2 \\
         *    1 & 3 & 4 \\
         *    2 & 4 & 5
         * \end{array}\right)\f].
         * The second dimension is either of size 1 in the case of elements
         * having #GeomType #eRegular, or of size equal to the number of
         * quadrature points for #eDeformed elements.
         *
         * @see [Wikipedia "Covariance and Contravariance of Vectors"]
         *      (http://en.wikipedia.org/wiki/Covariance_and_contravariance_of_vectors)
         * @returns             Two-dimensional array containing the inverse
         *                      metric tensor of the coordinate mapping.
         */
        Array<TwoD, NekDouble> GeomFactors::ComputeGmat(
                const LibUtilities::PointsKeyVector &keyTgt) const
        {
            ASSERTL1(keyTgt.size() == m_expDim,
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");

            int i = 0, j = 0, k = 0, l = 0;
            int ptsTgt   = 1;

            if (m_type == eDeformed)
            {
                // Allocate storage and compute number of points
                for (i = 0; i < m_expDim; ++i)
                {
                    ptsTgt   *= keyTgt[i].GetNumPoints();
                }
            }

            // Get derivative at geometry points
            DerivStorage deriv = ComputeDeriv(keyTgt);

            Array<TwoD, NekDouble> tmp (m_expDim*m_expDim, ptsTgt, 0.0);
            Array<TwoD, NekDouble> gmat(m_expDim*m_expDim, ptsTgt, 0.0);
            Array<OneD, NekDouble> jac (ptsTgt, 0.0);

            // Compute g_{ij} as t_i \cdot t_j and store in tmp
            for (i = 0, l = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j, ++l)
                {
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        Vmath::Vvtvp(ptsTgt, &deriv[i][k][0], 1,
                                             &deriv[j][k][0], 1,
                                             &tmp[l][0],      1,
                                             &tmp[l][0],      1);
                    }
                }
            }

            Adjoint(tmp, gmat);

            // Compute g = det(g_{ij}) (= Jacobian squared) and store
            // temporarily in m_jac.
            for (i = 0; i < m_expDim; ++i)
            {
                Vmath::Vvtvp(ptsTgt, &tmp[i][0], 1, &gmat[i*m_expDim][0], 1,
                                     &jac[0], 1, &jac[0], 1);
            }

            for (i = 0; i < m_expDim*m_expDim; ++i)
            {
                Vmath::Vdiv(ptsTgt, &gmat[i][0], 1, &jac[0], 1, &gmat[i][0], 1);
            }

            return gmat;
        }


        /**
         * @param   keyTgt      Target point distributions.
         * @returns             Derivative factors evaluated at the target
         *                      point distributions.
         * A 1D example: /f$ Jac =(\partial x/ \partial \xi) /f$ ; 
         *               /f$ factor = 1/Jac = (\partial \xi/ \partial x) /f$ 
         */
        Array<TwoD, NekDouble> GeomFactors::ComputeDerivFactors(
                const LibUtilities::PointsKeyVector& keyTgt) const
        {
            ASSERTL1(keyTgt.size() == m_expDim,
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");

            int i = 0, j = 0, k = 0, l = 0;
            int ptsTgt   = 1;

            if (m_type == eDeformed)
            {
                // Allocate storage and compute number of points
                for (i = 0; i < m_expDim; ++i)
                {
                    ptsTgt   *= keyTgt[i].GetNumPoints();
                }
            }

            // Get derivative at geometry points
            DerivStorage deriv = ComputeDeriv(keyTgt);

            Array<TwoD, NekDouble> tmp (m_expDim*m_expDim, ptsTgt, 0.0);
            Array<TwoD, NekDouble> gmat(m_expDim*m_expDim, ptsTgt, 0.0);
            Array<OneD, NekDouble> jac (ptsTgt, 0.0);
            Array<TwoD, NekDouble> factors(m_expDim*m_coordDim, ptsTgt, 0.0);

            // Compute g_{ij} as t_i \cdot t_j and store in tmp
            for (i = 0, l = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j, ++l)
                {
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        Vmath::Vvtvp(ptsTgt, &deriv[i][k][0], 1,
                                             &deriv[j][k][0], 1,
                                             &tmp[l][0],      1,
                                             &tmp[l][0],      1);
                    }
                }
            }

            Adjoint(tmp, gmat);

            // Compute g = det(g_{ij}) (= Jacobian squared) and store
            // temporarily in m_jac.
            for (i = 0; i < m_expDim; ++i)
            {
                Vmath::Vvtvp(ptsTgt, &tmp[i][0], 1, &gmat[i*m_expDim][0], 1,
                                     &jac[0],    1, &jac[0],              1);
            }

            for (i = 0; i < m_expDim*m_expDim; ++i)
            {
                Vmath::Vdiv(ptsTgt, &gmat[i][0], 1, &jac[0], 1, &gmat[i][0], 1);
            }

            // Compute the Jacobian = sqrt(g)
            Vmath::Vsqrt(ptsTgt, &jac[0], 1, &jac[0], 1);

            // Compute the derivative factors
            for (k = 0, l = 0; k < m_coordDim; ++k)
            {
                for (j = 0; j < m_expDim; ++j, ++l)
                {
                    for (i = 0; i < m_expDim; ++i)
                    {
                        Vmath::Vvtvp(ptsTgt, &deriv[i][k][0],        1,
                                             &gmat[m_expDim*i+j][0], 1,
                                             &factors[l][0],         1,
                                             &factors[l][0],         1);
                    }
                }
            }

            return factors;
        }

        void GeomFactors::ComputeMovingFrames(
            const LibUtilities::PointsKeyVector& keyTgt,
            const SpatialDomains::GeomMMF MMFdir,
            const Array<OneD, const NekDouble> &factors,
                  Array<OneD, Array<OneD, NekDouble> > &movingframes)
        {
            ASSERTL1(keyTgt.size() == m_expDim,
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");

            int i = 0, k = 0;
            int ptsTgt   = 1;
            int nq = 1;

            for (i = 0; i < m_expDim; ++i)
            {
                nq   *= keyTgt[i].GetNumPoints();
            }

            if (m_type == eDeformed)
            {
                // Allocate storage and compute number of points
                for (i = 0; i < m_expDim; ++i)
                {
                    ptsTgt   *= keyTgt[i].GetNumPoints();
                }
            }

            // Get derivative at geometry points
            DerivStorage deriv = ComputeDeriv(keyTgt);

            // number of moving frames is requited to be 3, even for surfaces
            int MFdim = 3;

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > MFtmp(MFdim);

            // Compute g_{ij} as t_i \cdot t_j and store in tmp
            for (i = 0; i < MFdim; ++i)
            {
                MFtmp[i] = Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
                for (k = 0; k < m_coordDim; ++k)
                {
                    MFtmp[i][k] = Array<OneD, NekDouble>(nq);
                }
            }

            // Compute g_{ij} as t_i \cdot t_j and store in tmp
            for (i = 0; i < MFdim-1; ++i)
            {
                for (k = 0; k < m_coordDim; ++k)
                {
                    if (m_type == eDeformed)
                    {
                        Vmath::Vcopy(ptsTgt, &deriv[i][k][0], 1,
                                             &MFtmp[i][k][0], 1);
                    }
                    else
                    {
                        Vmath::Fill(nq, deriv[i][k][0], MFtmp[i][k], 1);
                    }
                }
            }

            // Direction of MF1 is preserved: MF2 is considered in the same
            // tangent plane as MF1. MF3 is computed by cross product of MF1
            // and MF2. MF2 is consequently computed as the cross product of
            // MF3 and MF1.
            Array<OneD, Array<OneD,NekDouble> > PrincipleDir(m_coordDim);
            for (int k=0; k<m_coordDim; k++)
            {
                PrincipleDir[k] = Array<OneD, NekDouble>(nq);
            }

            if( !(MMFdir == eLOCAL) )
            {
                ComputePrincipleDirection(keyTgt, MMFdir,
                                          factors, PrincipleDir);
            }

            // MF3 = MF1 \times MF2
            VectorCrossProd(MFtmp[0], MFtmp[1], MFtmp[2]);

            // Normalizing MF3
            VectorNormalise(MFtmp[2]);

            if ( !(MMFdir == eLOCAL) )
            {
                Array<OneD, NekDouble> temp(nq, 0.0);

                // Reorient MF1 along the PrincipleDir
                for (i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vvtvp(nq, MFtmp[2][i],     1,
                                     PrincipleDir[i], 1,
                                     temp,            1,
                                     temp,            1);
                }
                Vmath::Neg(nq, temp, 1);

                // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
                for (i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vvtvp(nq, temp,            1,
                                     MFtmp[2][i],     1,
                                     PrincipleDir[i], 1,
                                     MFtmp[0][i],     1);
                }
            }

            // Normalizing MF1
            VectorNormalise(MFtmp[0]);

            // MF2 = MF3 \times MF1
            VectorCrossProd(MFtmp[2], MFtmp[0], MFtmp[1]);

            // Normalizing MF2
            VectorNormalise(MFtmp[1]);

            for (i = 0; i < MFdim; ++i)
            {
                for (k = 0; k < m_coordDim; ++k)
                {
                    Vmath::Vcopy(nq, &MFtmp[i][k][0],                  1,
                                     &movingframes[i*m_coordDim+k][0], 1);
                }
            }
        }


        /**
         * Constructs the Jacobian as per Spencer's book p158 and tests if
         * negative.
         */
        void GeomFactors::CheckIfValid()
        {
            // Jacobian test only makes sense when expdim = coorddim
            // If one-dimensional then element is valid.
            if (m_coordDim != m_expDim || m_expDim == 1)
            {
                m_valid = true;
                return;
            }

            LibUtilities::PointsKeyVector p(m_expDim);
            int nqtot = 1;
            for (int i = 0; i < m_expDim; ++i)
            {
                p[i] = m_xmap->GetBasis(i)->GetPointsKey();
                nqtot *= p[i].GetNumPoints();
            }
            int pts = (m_type == eRegular || m_type == eMovingRegular)
                            ? 1 : nqtot;
            Array<OneD, NekDouble> jac(pts, 0.0);

            DerivStorage deriv = GetDeriv(p);

            switch (m_expDim)
            {
                case 2:
                {
                    Vmath::Vvtvvtm(pts, &deriv[0][0][0], 1, &deriv[1][1][0], 1,
                                        &deriv[1][0][0], 1, &deriv[0][1][0], 1,
                                        &jac[0],         1);
                    break;
                }
                case 3:
                {
                    Array<OneD, NekDouble> tmp(pts, 0.0);

                    Vmath::Vvtvvtm(pts, &deriv[1][1][0], 1, &deriv[2][2][0], 1,
                                        &deriv[2][1][0], 1, &deriv[1][2][0], 1,
                                        &tmp[0],         1);
                    Vmath::Vvtvp  (pts, &deriv[0][0][0], 1, &tmp[0],         1,
                                        &jac[0],         1, &jac[0],         1);

                    Vmath::Vvtvvtm(pts, &deriv[2][1][0], 1, &deriv[0][2][0], 1,
                                        &deriv[0][1][0], 1, &deriv[2][2][0], 1,
                                        &tmp[0],         1);
                    Vmath::Vvtvp  (pts, &deriv[1][0][0], 1, &tmp[0],         1,
                                        &jac[0],         1, &jac[0],         1);

                    Vmath::Vvtvvtm(pts, &deriv[0][1][0], 1, &deriv[1][2][0], 1,
                                        &deriv[1][1][0], 1, &deriv[0][2][0], 1,
                                        &tmp[0],         1);
                    Vmath::Vvtvp  (pts, &deriv[2][0][0], 1, &tmp[0],         1,
                                        &jac[0],         1, &jac[0],         1);

                    break;
                }
            }

            if (Vmath::Vmin(pts, &jac[0], 1) < 0)
            {
                m_valid = false;
            }
        }


        /**
         * @param   map_points  Source data point distribution.
         * @param   src         Source data to be interpolated.
         * @param   tpoints     Target data point distribution.
         * @param   tgt         Target data storage.
         */
        void GeomFactors::Interp(
                    const LibUtilities::PointsKeyVector &src_points,
                    const Array<OneD, const NekDouble> &src,
                    const LibUtilities::PointsKeyVector &tgt_points,
                    Array<OneD, NekDouble> &tgt) const
        {
            ASSERTL1(src_points.size() == tgt_points.size(),
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");

            switch (m_expDim)
            {
                case 1:
                    LibUtilities::Interp1D(src_points[0], src,
                                           tgt_points[0], tgt);
                    break;
                case 2:
                    LibUtilities::Interp2D(src_points[0], src_points[1], src,
                                           tgt_points[0], tgt_points[1], tgt);
                    break;
                case 3:
                    LibUtilities::Interp3D(src_points[0], src_points[1],
                                           src_points[2], src,
                                           tgt_points[0], tgt_points[1],
                                           tgt_points[2], tgt);
                    break;
            }
        }


        /**
         * Input and output arrays are of dimension
         * (m_expDim*m_expDim) x num_points. The first index of the input and
         * output arrays are ordered row-by-row.
         * @param   src         Input data array.
         * @param   tgt         Storage for adjoint matrix data.
         */
        void GeomFactors::Adjoint(
                    const Array<TwoD, const NekDouble>& src,
                    Array<TwoD, NekDouble>& tgt) const
        {
            ASSERTL1(src.size() == tgt.size(),
                     "Source matrix is of different size to destination"
                     "matrix for computing adjoint.");

            int n = src[0].size();
            switch (m_expDim)
            {
                case 1:
                    Vmath::Fill (n,  1.0, &tgt[0][0], 1);
                    break;
                case 2:
                    Vmath::Vcopy(n,       &src[3][0], 1, &tgt[0][0], 1);
                    Vmath::Smul (n, -1.0, &src[1][0], 1, &tgt[1][0], 1);
                    Vmath::Smul (n, -1.0, &src[2][0], 1, &tgt[2][0], 1);
                    Vmath::Vcopy(n,       &src[0][0], 1, &tgt[3][0], 1);
                    break;
                case 3:
                {
                    int a, b, c, d, e, i, j;

                    // Compute g^{ij} by computing Cofactors(g_ij)^T
                    for (i = 0; i < m_expDim; ++i)
                    {
                        for (j = 0; j < m_expDim; ++j)
                        {
                            a = ((i+1)%m_expDim) * m_expDim + ((j+1)%m_expDim);
                            b = ((i+1)%m_expDim) * m_expDim + ((j+2)%m_expDim);
                            c = ((i+2)%m_expDim) * m_expDim + ((j+1)%m_expDim);
                            d = ((i+2)%m_expDim) * m_expDim + ((j+2)%m_expDim);
                            e = j*m_expDim + i;
                            Vmath::Vvtvvtm(n, &src[a][0], 1, &src[d][0], 1,
                                              &src[b][0], 1, &src[c][0], 1,
                                              &tgt[e][0], 1);
                        }
                    }
                    break;
                }
            }
        }


        /**
         *
         */
        void GeomFactors::ComputePrincipleDirection(
            const LibUtilities::PointsKeyVector& keyTgt,
            const SpatialDomains::GeomMMF MMFdir,
            const Array<OneD, const NekDouble> &factors,
                  Array<OneD, Array<OneD,NekDouble> > &output)
        {
            int nq = output[0].size();

            output = Array<OneD,Array<OneD,NekDouble> >(m_coordDim);
            for (int i = 0; i < m_coordDim; ++i)
            {
                output[i] = Array<OneD, NekDouble> (nq, 0.0);
            }

            // Construction of Connection
            switch(MMFdir)
            {
                // projection to x-axis
                case eTangentX:
                {
                    Vmath::Fill(nq, 1.0, output[0], 1);
                    break;
                }
                case eTangentY:
                {
                    Vmath::Fill(nq, 1.0, output[1], 1);
                    break;
                }
                case eTangentXY:
                {
                    Vmath::Fill(nq, sqrt(2.0), output[0], 1);
                    Vmath::Fill(nq, sqrt(2.0), output[1], 1);
                    break;
                }
                case eTangentZ:
                {
                    Vmath::Fill(nq, 1.0, output[2], 1);
                    break;
                }
                case eTangentCircular:
                {
                    // Tangent direction depends on spatial location.
                    Array<OneD, Array<OneD, NekDouble> > x(m_coordDim);
                    for (int k = 0; k < m_coordDim; k++)
                    {
                        x[k] = Array<OneD, NekDouble>(nq);
                    }

                    // m_coords are StdExpansions which store the mapping
                    // between the std element and the local element. Bwd
                    // transforming the std element minimum basis gives a
                    // minimum physical basis for geometry. Need to then
                    // interpolate this up to the quadrature basis.
                    int nqtot_map = 1;
                    LibUtilities::PointsKeyVector map_points(m_expDim);
                    for (int i = 0; i < m_expDim; ++i)
                    {
                        map_points[i]  = m_xmap->GetBasis(i)->GetPointsKey();
                        nqtot_map     *= map_points[i].GetNumPoints();
                    }
                    Array<OneD, NekDouble> tmp(nqtot_map);
                    for (int k = 0; k < m_coordDim; k++)
                    {
                        m_xmap->BwdTrans(m_coords[k], tmp);
                        Interp(map_points, tmp, keyTgt, x[k]);
                    }

                    // circular around the center of the domain
                    NekDouble radius, xc=0.0, yc=0.0, xdis, ydis;
                    NekDouble la, lb;

                    ASSERTL1(factors.size() >= 4,
                             "factors is too short.");

                    la = factors[0];
                    lb = factors[1];
                    xc = factors[2];
                    yc = factors[3];

                    for (int i = 0; i < nq; i++)
                    {
                        xdis = x[0][i]-xc;
                        ydis = x[1][i]-yc;
                        radius = sqrt(xdis*xdis/la/la+ydis*ydis/lb/lb);
                        output[0][i] = ydis/radius;
                        output[1][i] = -1.0*xdis/radius;
                    }
                    break;
                }
                case eTangentIrregular:
                {
                    // Tangent direction depends on spatial location.
                    Array<OneD, Array<OneD, NekDouble> > x(m_coordDim);
                    for (int k = 0; k < m_coordDim; k++)
                    {
                        x[k] = Array<OneD, NekDouble>(nq);
                    }

                    int nqtot_map = 1;
                    LibUtilities::PointsKeyVector map_points(m_expDim);
                    for (int i = 0; i < m_expDim; ++i)
                    {
                        map_points[i]  = m_xmap->GetBasis(i)->GetPointsKey();
                        nqtot_map     *= map_points[i].GetNumPoints();
                    }
                    Array<OneD, NekDouble> tmp(nqtot_map);
                    for (int k = 0; k < m_coordDim; k++)
                    {
                        m_xmap->BwdTrans(m_coords[k], tmp);
                        Interp(map_points, tmp, keyTgt, x[k]);
                    }

                    // circular around the center of the domain
                    NekDouble xtan, ytan, mag;
                    for (int i = 0; i < nq; i++)
                    {
                        xtan = -1.0*(x[1][i]*x[1][i]*x[1][i] + x[1][i]);
                        ytan = 2.0*x[0][i];
                        mag  = sqrt(xtan*xtan + ytan*ytan);
                        output[0][i] = xtan/mag;
                        output[1][i] = ytan/mag;
                    }
                    break;
                }
                case eTangentNonconvex:
                {
                    // Tangent direction depends on spatial location.
                    Array<OneD, Array<OneD, NekDouble> > x(m_coordDim);
                    for (int k = 0; k < m_coordDim; k++)
                    {
                        x[k] = Array<OneD, NekDouble>(nq);
                    }

                    int nqtot_map = 1;
                    LibUtilities::PointsKeyVector map_points(m_expDim);
                    for (int i = 0; i < m_expDim; ++i)
                    {
                        map_points[i]  = m_xmap->GetBasis(i)->GetPointsKey();
                        nqtot_map     *= map_points[i].GetNumPoints();
                    }
                    Array<OneD, NekDouble> tmp(nqtot_map);
                    for (int k = 0; k < m_coordDim; k++)
                    {
                        m_xmap->BwdTrans(m_coords[k], tmp);
                        Interp(map_points, tmp, keyTgt, x[k]);
                    }

                    // circular around the center of the domain
                    NekDouble xtan, ytan, mag;
                    for (int i = 0; i < nq; i++)
                    {
                        xtan = -2.0*x[1][i]*x[1][i]*x[1][i] + x[1][i];
                        ytan = sqrt(3.0)*x[0][i];
                        mag  = sqrt(xtan*xtan + ytan*ytan);
                        output[0][i] = xtan/mag;
                        output[1][i] = ytan/mag;
                    }
                    break;
                }
                default:
                {
                    break;
                }
            }
        }


        /**
         *
         */
        void GeomFactors::VectorNormalise(
            Array<OneD, Array<OneD, NekDouble> > &array)
        {
            int ndim = array.size();
            ASSERTL0(ndim > 0, "Number of components must be > 0.");
            for (int i = 1; i < ndim; ++i)
            {
                ASSERTL0(array[i].size() == array[0].size(),
                         "Array size mismatch in coordinates.");
            }

            int nq        = array[0].size();
            Array<OneD, NekDouble> norm (nq, 0.0);

            // Compute the norm of each vector.
            for (int i = 0; i < ndim; ++i)
            {
                Vmath::Vvtvp(nq, array[i],  1,
                                 array[i],  1,
                                 norm,      1,
                                 norm,      1);
            }

            Vmath::Vsqrt(nq, norm, 1, norm, 1);

            // Normalise the vectors by the norm
            for (int i = 0; i < ndim; ++i)
            {
                Vmath::Vdiv(nq, array[i], 1, norm, 1, array[i], 1);
            }
        }


        /**
         * @brief Computes the vector cross-product in 3D of \a v1 and \a v2,
         * storing the result in \a v3.
         *
         * @param   v1          First input vector.
         * @param   v2          Second input vector.
         * @param   v3          Output vector computed to be orthogonal to
         *                      both \a v1 and \a v2.
         */
        void GeomFactors::VectorCrossProd(
            const Array<OneD, const Array<OneD, NekDouble> > &v1,
            const Array<OneD, const Array<OneD, NekDouble> > &v2,
                  Array<OneD,       Array<OneD, NekDouble> > &v3)
        {
            ASSERTL0(v1.size() == 3,
                     "Input 1 has dimension not equal to 3.");
            ASSERTL0(v2.size() == 3,
                     "Input 2 has dimension not equal to 3.");
            ASSERTL0(v3.size() == 3,
                     "Output vector has dimension not equal to 3.");

            int nq = v1[0].size();
            Array<OneD, NekDouble> temp(nq);

            Vmath::Vmul (nq, v1[2], 1, v2[1], 1, temp, 1);
            Vmath::Vvtvm(nq, v1[1], 1, v2[2], 1, temp, 1, v3[0], 1);

            Vmath::Vmul (nq, v1[0], 1, v2[2], 1, temp, 1);
            Vmath::Vvtvm(nq, v1[2], 1, v2[0], 1, temp, 1, v3[1], 1);

            Vmath::Vmul (nq, v1[1], 1, v2[0], 1, temp, 1);
            Vmath::Vvtvm(nq, v1[0], 1, v2[1], 1, temp, 1, v3[2], 1);
        }

    } //end of namespace
} //end of namespace
