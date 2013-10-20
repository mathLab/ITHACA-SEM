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
         */

        /**
         * The argument 'tbasis' contains the information about the quadrature
         * points on which the weighted metric terms should be specified.
         * The geometric factors are evaluated by considering the mapping to a
         * coordinate system based on the local tangent vectors (which are the
         * local derivatives of the global coordinates) and the normal
         * \f$ \bf g \f$ to these two tangent vectors. We therefore use the
         * 3 x 3 relationships but assume that
         * \f$ \partial x_1/\partial \xi_3 = { g_1},\,
         * \partial x_2/\partial \xi_3 = { g_2},\, \partial x_3/\partial
         * \xi_3 = { g_3} \f$ i.e.
         *
         * \f$ {\bf g }= \left [ \begin{array}{c} g_1 \\ g_2 \\ g_3 \end{array}
         * \right ] = \frac{\partial {\bf x}}{\partial \xi_1} \times
         * \frac{\partial {\bf x}}{\partial \xi_2} =   \left [ \begin{array}{c}
         * \frac{\partial x_2}{\partial\xi_1}\frac{\partial x_3}{\partial\xi_2}-
         * \frac{\partial x_3}{\partial\xi_1}\frac{\partial x_2}{\partial\xi_2}
         * \\
         * \frac{\partial x_3}{\partial\xi_1}\frac{\partial x_1}{\partial\xi_2}-
         * \frac{\partial x_1}{\partial\xi_1}\frac{\partial x_3}{\partial\xi_2}
         * \\
         * \frac{\partial x_1}{\partial\xi_1}\frac{\partial x_2}{\partial\xi_2}-
         * \frac{\partial x_2}{\partial\xi_1}\frac{\partial x_1}{\partial\xi_2}
         * \end{array} \right ] \f$
         *
         * The geometric factors are then given by:
         *
         * \f$ \begin{array}{cc}
         * \frac{\partial \xi_1}{\partial x_1} = \frac{1}{J_{3D}} \left (
         *      \frac{\partial x_2}{\partial \xi_2} {g_3}
         *      - \frac{\partial x_3}{\partial \xi_2} {g_2} \right ) &
         * \frac{\partial \xi_1}{\partial x_2} = -\frac{1}{J_{3D}} \left (
         *      \frac{\partial x_1}{\partial \xi_2} {g_3}
         *      - \frac{\partial x_3}{\partial \xi_2} {g_1} \right )\\
         * \frac{\partial \xi_1}{\partial x_3} = \frac{1}{J_{3D}} \left (
         *      \frac{\partial x_1}{\partial \xi_2} {g_2}
         *      - \frac{\partial x_2}{\partial \xi_2} {g_1} \right )&
         * \frac{\partial \xi_2}{\partial x_1} = -\frac{1}{J_{3D}} \left (
         *      \frac{\partial x_2}{\partial \xi_1} {g_3}
         *      - \frac{\partial x_3}{\partial \xi_1} {g_2} \right ) \\
         * \frac{\partial \xi_2}{\partial x_2} = \frac{1}{J_{3D}} \left (
         *      \frac{\partial x_1}{\partial \xi_1} {g_3}
         *      - \frac{\partial x_3}{\partial \xi_1} {g_1} \right ) &
         * \frac{\partial \xi_2}{\partial x_3} = -\frac{1}{J_{3D}} \left (
         *      \frac{\partial x_1}{\partial \xi_1} {g_2}
         *      - \frac{\partial x_2}{\partial \xi_1} {g_1} \right )
         * \end{array} \f$
         *
         * where
         *
         * \f$ J_{3D} =
         * {g_3} \left (
         * \frac{\partial x_1}{\partial\xi_1}\frac{\partial x_2}{\partial\xi_2}-
         * \frac{\partial x_1}{\partial\xi_2}\frac{\partial x_2}{\partial\xi_1}
         * \right ) + {g_2}\left (
         * \frac{\partial x_1}{\partial\xi_2}\frac{\partial x_3}{\partial\xi_1}-
         * \frac{\partial x_1}{\partial\xi_1}\frac{\partial x_3}{\partial\xi_2}
         * \right ) + {g_1} \left (
         * \frac{\partial x_2}{\partial\xi_1}\frac{\partial x_3}{\partial\xi_2}-
         * \frac{\partial x_2}{\partial\xi_2}\frac{\partial x_3}{\partial\xi_1}
         * \right ) \f$
         *
         * and the two-dimensional surface Jacobian  is given by
         * \f$ J = \left | \frac{\partial {\bf x}}{\partial \xi_1} \times
         * \frac{\partial {\bf x}}{\partial \xi_2} \right | = \sqrt{J_{3D}} \f$
         *
         * @param   gtype       Type of geometry.
         * @param   coordim     Dimension of coordinate system.
         * @param   Coords      Coordinate information.
         * @param   tbasis      Tangent basis.
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

            // Copy shared pointers.
            for (int i = 0; i < m_coordDim; ++i)
            {
                m_coords[i] = Coords[i];
            }

            // The quadrature points of the mapping
            // (as specified in Coords)
            LibUtilities::PointsKey pkey0_map(
                                        Coords[0]->GetBasis(0)->GetPointsKey());
            LibUtilities::PointsKey pkey1_map(
                                        Coords[0]->GetBasis(1)->GetPointsKey());
            int nquad0_map = pkey0_map.GetNumPoints();
            int nquad1_map = pkey1_map.GetNumPoints();
            int nqtot_map  = nquad0_map*nquad1_map;

            // The quadrature points at the points at which we
            // want to know the metrics (as specified in tbasis)
            LibUtilities::PointsKey pkey0_tbasis(tbasis[0]->GetPointsKey());
            LibUtilities::PointsKey pkey1_tbasis(tbasis[1]->GetPointsKey());
            int nquad0_tbasis = pkey0_tbasis.GetNumPoints();
            int nquad1_tbasis = pkey1_tbasis.GetNumPoints();
            int nqtot_tbasis  = nquad0_tbasis*nquad1_tbasis;

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            m_pointsKey[0] = pkey0_tbasis;
            m_pointsKey[1] = pkey1_tbasis;

            // setup temp storage
            Array<OneD, Array<OneD,NekDouble> > d1_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d2_map   (coordim);

            m_deriv = Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(m_expDim);
            m_deriv[0] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
            m_deriv[1] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);

            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                d1_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d2_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                m_deriv[0][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                m_deriv[1][i] = Array<OneD,NekDouble>(nqtot_tbasis);

                // Transform from coefficient space to physical space
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),
                                    Coords[i]->UpdatePhys());
                // Take the derivative (calculated at the points as specified
                // in 'Coords')
                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),
                                        d1_map[i],
                                        d2_map[i]);

                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics
                //   ('tbasis')
                if( (pkey0_map == pkey0_tbasis) &&
                    (pkey1_map == pkey1_tbasis) )
                {
                    m_deriv[0][i] = d1_map[i];
                    m_deriv[1][i] = d2_map[i];
                }
                else
                {
                    LibUtilities::Interp2D(pkey0_map, pkey1_map, d1_map[i],
                                    pkey0_tbasis, pkey1_tbasis, m_deriv[0][i]);
                    LibUtilities::Interp2D(pkey0_map, pkey1_map, d2_map[i],
                                    pkey0_tbasis, pkey1_tbasis, m_deriv[1][i]);
                }
            }

            // Setting up Surface Normal Vectors
            /*if(coordim == 3)
            {
                v_ComputeSurfaceNormals();
            }*/

            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation
            // metrics
            SetUpJacGmat2D();

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
            // Check the number of derivative dimensions matches the coordinate
            // space.
            ASSERTL1(m_deriv[0].num_elements()==m_coordDim,
                     "The dimension of array d1 does not match the coordinate "
                     "dimension");
            ASSERTL1(m_deriv[1].num_elements()==m_coordDim,
                     "The dimension of array d2 does not match the coordinate "
                     "dimension");

            unsigned int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();
            unsigned int pts = 1;
            unsigned int i, j, k;

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

            pts = (m_type == eRegular || m_type == eMovingRegular) ? 1 : nqtot;

            // Jacobian is constant across the element.
            m_jac     = Array<OneD, NekDouble>(pts,0.0);

            // Number of entries corresponds to twice the coordinate
            // dimension. Entries are constant across element so second
            // dimension is 1.
            m_gmat    = Array<TwoD, NekDouble>(m_expDim*m_expDim, pts, 0.0);

            m_derivFactors = Array<TwoD, NekDouble>(m_expDim*m_coordDim, pts, 0.0);

            Array<TwoD, NekDouble> tmp(m_expDim*m_expDim, pts, 0.0);

            // Compute g_{ij} as t_i \cdot t_j
            for (i = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j)
                {
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        Vmath::Vvtvp(pts, &m_deriv[i][k][0], 1, &m_deriv[j][k][0], 1, &tmp[m_expDim*i+j][0], 1, &tmp[m_expDim*i+j][0], 1);
                    }
                }
            }

            // Compute g (Jacobian squared)
            Vmath::Vvtvvtm(pts, &tmp[0][0], 1, &tmp[3][0], 1, &tmp[1][0], 1, &tmp[2][0], 1, &m_jac[0], 1);

            // Compute g^{ij}
            Vmath::Vdiv(pts, &tmp[3][0], 1, &m_jac[0], 1, &m_gmat[0][0], 1);
            Vmath::Vdiv(pts, &tmp[1][0], 1, &m_jac[0], 1, &m_gmat[1][0], 1);
            Vmath::Vdiv(pts, &tmp[2][0], 1, &m_jac[0], 1, &m_gmat[2][0], 1);
            Vmath::Vdiv(pts, &tmp[0][0], 1, &m_jac[0], 1, &m_gmat[3][0], 1);
            Vmath::Smul (pts, -1.0, &m_gmat[1][0], 1, &m_gmat[1][0], 1);
            Vmath::Smul (pts, -1.0, &m_gmat[2][0], 1, &m_gmat[2][0], 1);

            // Sqrt jacobian
            Vmath::Vsqrt(pts, &m_jac[0], 1, &m_jac[0], 1);

            for (i = 0; i < m_expDim; ++i)
            {
                for (j = 0; j < m_expDim; ++j)
                {
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        Vmath::Vvtvp(pts, &m_deriv[i][k][0], 1, &m_gmat[m_expDim*i+j][0], 1, &m_derivFactors[m_expDim*k+j][0], 1, &m_derivFactors[m_expDim*k+j][0], 1);
                    }
                }
            }
        }


        /**
         * Check if the element is valid. This check only applies to elements
         * in a 2D coordinate space.
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

            Array<OneD, NekDouble> jac(pts, 0.0);
            Vmath::Vvtvvtm(pts, &m_deriv[0][0][0], 1, &m_deriv[1][1][0], 1,
                                &m_deriv[1][0][0], 1, &m_deriv[0][1][0], 1,
                                &jac[0],           1);

            if(Vmath::Vmin(pts, &jac[0], 1) < 0)
            {
                m_valid = false;
            }
        }


        /**
         * The principle direction may be chosen to ensure continuity of the
         * tangent vector across elements. This will typically provide improved
         * performance, but at the cost of potentially introducing
         * singularities due to the geometry. The possible choices are
         * - unit vectors in x, y, or z directions.
         * - unit vectors encircling a hole in the surface.
         *
         * @param   output      Storage for principle direction vector.
         */
        void GeomFactors2D::ComputePrincipleDirection(
                        Array<OneD,Array<OneD,NekDouble> > &output)
        {
            int nq = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();

            output = Array<OneD,Array<OneD,NekDouble> >(m_coordDim);
            for (int i = 0; i < m_coordDim; ++i)
            {
                output[i] = Array<OneD, NekDouble> (nq, 0.0);
            }

            // Construction of Connection
            switch(m_tangentDir)
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
                case eTangentZ:
                {
                    Vmath::Fill(nq, 1.0, output[2], 1);
                    break;
                }
                case eTangentCircular:
                {
                    // Tangent direction depends on spatial location.
                    Array<OneD,NekDouble> x0(nq);
                    Array<OneD,NekDouble> x1(nq);
                    Array<OneD,NekDouble> x2(nq);

                    // m_coords are StdExpansions which store the mapping
                    // between the std element and the local element. Bwd
                    // transforming the std element minimum basis gives a
                    // minimum physical basis for geometry. Need to then
                    // interpolate this up to the quadrature basis.
                    LibUtilities::Interp2D(
                                    m_coords[0]->GetBasis(0)->GetPointsKey(),
                                    m_coords[0]->GetBasis(1)->GetPointsKey(),
                                    m_coords[0]->GetPhys(),
                                    m_pointsKey[0], m_pointsKey[1], x0);
                    LibUtilities::Interp2D(
                                    m_coords[0]->GetBasis(0)->GetPointsKey(),
                                    m_coords[0]->GetBasis(1)->GetPointsKey(),
                                    m_coords[1]->GetPhys(),
                                    m_pointsKey[0], m_pointsKey[1], x1);
                    LibUtilities::Interp2D(
                                    m_coords[0]->GetBasis(0)->GetPointsKey(),
                                    m_coords[0]->GetBasis(1)->GetPointsKey(),
                                    m_coords[2]->GetPhys(),
                                    m_pointsKey[0], m_pointsKey[1], x2);

                    // circular around the center of the domain
                    NekDouble radius, xc=0.0, yc=0.0, xdis, ydis;

                    if (m_tangentDirCentre.num_elements() == 2) {
                        xc = m_tangentDirCentre[0];
                        yc = m_tangentDirCentre[1];
                    }

                    for (int i = 0; i < nq; i++)
                    {
                        xdis = x0[i]-xc;
                        ydis = x1[i]-yc;
                        radius = sqrt(xdis*xdis+ydis*ydis);
                        output[0][i] = ydis/radius;
                        output[1][i] = -1.0*xdis/radius;
                    }
                }
                default:
                {
                    ASSERTL0(false, "Unsupported tangent direction.");
                    break;
                }
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
            // Number of dimensions.
            int coordim = m_deriv[0].num_elements();
            // Total number of points in the element.
            int nqtot = m_deriv[0][0].num_elements();

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
            Vmath::Vmul (nqtot, m_deriv[0][2],     1,
                                m_deriv[1][1],     1,
                                temp,             1);
            Vmath::Vvtvm(nqtot, m_deriv[0][1],     1,
                                m_deriv[1][2],     1,
                                temp,             1,
                                m_normal[0],          1);

            Vmath::Vmul (nqtot, m_deriv[0][0],     1,
                                m_deriv[1][2],     1,
                                temp,           1);
            Vmath::Vvtvm(nqtot, m_deriv[0][2],     1,
                                m_deriv[1][0],     1,
                                temp,           1,
                                m_normal[1],    1);

            Vmath::Vmul (nqtot, m_deriv[0][1],     1,
                                m_deriv[1][0],     1,
                                temp,           1);
            Vmath::Vvtvm(nqtot, m_deriv[0][0],     1,
                                m_deriv[1][1],     1,
                                temp,           1,
                                m_normal[2],  1);

            // Reset temp array
            temp = Array<OneD, NekDouble>(nqtot,0.0);

            // Normalization of Surface Normal
            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vvtvp(nqtot, m_normal[i],  1,
                                    m_normal[i],  1,
                                    temp,           1,
                                    temp,           1);
            }
            Vmath::Vsqrt(nqtot, temp, 1, temp, 1);

            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vdiv(nqtot,  m_normal[i],  1,
                                    temp,           1,
                                    m_normal[i],  1);
            }
        }


        /**
         * Computes the tangent vector at each quadrature point in the
         * expansion using the surface normal and a principle direction. The
         * surface normal, and thus the tangent plane is uniquely defined
         * (modulo sign of the normal). A chosen principle direction is first
         * made orthonormal to the surface normal via Gram-Schmidt to produce
         * the first tangent vector. The second tangent vector is computed as
         * the cross-product of the first with the surface normal.
         */
        void GeomFactors2D::v_ComputeTangents()
        {
            int i, j, k;
            int nq = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();

            if (m_normal.num_elements() == 0) {
                v_ComputeSurfaceNormals();
            }

            if(m_tangentDir==eLOCAL)
            {
                // Allocate temporary and tangent storage.
                Array<OneD, NekDouble> temp(nq, 0.0);
                m_tangents  =  Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);
                for (i = 0; i < 2; ++i)
                {
                    m_tangents[i] = Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        m_tangents[i][k] = Array<OneD, NekDouble>(nq, 0.0);
                    }
                }

                for (i = 0; i < m_coordDim; ++i)
                {
                    for(j =0; j < 2; ++j)
                    {
                        Vmath::Vcopy(nq, m_deriv[0][i], 1, m_tangents[0][i], 1);
                    }
                }

                // Normalise first tangent vectors.
                VectorNormalise(m_tangents[0]);

                // The second set of tangential vectors is obtained by cross-
                // producting the surface normal with the first tangential vectors.
                VectorCrossProd(m_tangents[0], m_normal, m_tangents[1]);

                // Normalise second tangent vectors.
                VectorNormalise(m_tangents[1]);
            }

            else
            {
                // Get a principle direction for the expansion.
                Array<OneD, Array<OneD,NekDouble> > PrincipleDir;
                ComputePrincipleDirection(PrincipleDir);

                // Allocate temporary and tangent storage.
                Array<OneD, NekDouble> temp(nq, 0.0);
                m_tangents  =  Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);
                for (i = 0; i < 2; ++i)
                {
                    m_tangents[i] = Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
                    for (k = 0; k < m_coordDim; ++k)
                    {
                        m_tangents[i][k] = Array<OneD, NekDouble>(nq, 0.0);
                    }
                }

                // Gram-schmidz process to make this principaldirection orthonormal
                // to surface normal vectors. Let u1 = v1 = SurfaceNormal,
                // v2 = Principaldirection, which should be orthogornalised to u2
                // inner12 = < u1, v2 >, norm2 = < u1, u1 > = 1 by default
                // ipro = - < u1, v2 > / < u1, u1 >
                for (i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vvtvp(nq, m_normal[i],        1,
                                 PrincipleDir[i],   1,
                                 temp,              1,
                                 temp,              1);
                }
                Vmath::Neg(nq, temp, 1);

                // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
                for (i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vvtvp(nq, temp,              1,
                                 m_normal[i],        1,
                                 PrincipleDir[i],   1,
                                 m_tangents[0][i],   1);
                }

                // Normalise first tangent vectors.
                VectorNormalise(m_tangents[0]);

                // The second set of tangential vectors is obtained by cross-
                // producting the surface normal with the first tangential vectors.
                VectorCrossProd(m_tangents[0], m_normal, m_tangents[1]);

                // Normalise second tangent vectors.
                VectorNormalise(m_tangents[1]);
            }
        }

    }
}
