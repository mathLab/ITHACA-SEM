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
         * @param   SetUpQuadratureMetrics  Compute quadrature metrics.
         * @param   SetUpLaplacianMetrics   Compute Laplacian metrics.
         */
        GeomFactors2D::GeomFactors2D(
                        const GeomType gtype,
                        const int coordim,
                        const Array<OneD, const StdRegions
                                            ::StdExpansion2DSharedPtr> &Coords,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis,
                        const bool QuadMetrics,
                        const bool LaplMetrics,
                        const bool CheckJacPositive ) :
            GeomFactors(gtype,2,coordim,QuadMetrics,LaplMetrics)
        {
            // Sanity checks.
            ASSERTL1((coordim == 2)||(coordim == 3),
                     "The coordinate dimension should be equal to two or three"
                     "for two-dimensional elements");
            ASSERTL1(tbasis.num_elements()==2,
                     "tbasis should be an array of size two");
            ASSERTL1(LaplMetrics?QuadMetrics:true,
                     "SetUpQuadratureMetrics should be true if "
                     "SetUpLaplacianMetrics is true");

            // Copy shared pointers.
            for (int i = 0; i < m_coordDim; ++i)
            {
                m_coords[i] = Coords[i];
            }

            LibUtilities::ShapeType shape = Coords[0]->DetShapeType();

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

            m_deriv = Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);
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
            SetUpJacGmat2D(CheckJacPositive);

            // 2. the jacobian muliplied with the quadrature weights
            if(QuadMetrics)
            {
                SetUpQuadratureMetrics(shape,tbasis);
            }
            // 3. A combination of the metrics above that allows
            //    for more efficient evaluation of the laplacian
            if(LaplMetrics)
            {
                SetUpLaplacianMetrics(shape,tbasis);
            }
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
        void GeomFactors2D::SetUpJacGmat2D(bool CheckJacPositive)
        {
            // Check the number of derivative dimensions matches the coordinate
            // space.
            ASSERTL1(m_deriv[0].num_elements()==m_coordDim,
                     "The dimension of array d1 does not match the coordinate "
                     "dimension");
            ASSERTL1(m_deriv[1].num_elements()==m_coordDim,
                     "The dimension of array d2 does not match the coordinate "
                     "dimension");

            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();

            // Check each derivative combination has the correct sized storage
            // for all quadrature points.
            for (int i = 0; i < 2; ++i)
            {
                for (int j = 0; j < m_coordDim; ++j)
                {
                    ASSERTL1(m_deriv[i][j].num_elements() == nqtot,
                             "Number of quadrature points do not match");
                }
            }

            // Proceed differently depending on whether the geometry is
            // regular or deformed.
            if((m_type == eRegular)||(m_type == eMovingRegular))
            {
                // Jacobian is constant across the element.
                m_jac     = Array<OneD, NekDouble>(1,0.0);

                // Number of entries corresponds to twice the coordinate
                // dimension. Entries are constant across element so second
                // dimension is 1.
                m_gmat    = Array<TwoD, NekDouble>(2*m_coordDim,1,0.0);

                if(m_coordDim == 2) // assume g = [0,0,1]
                {
                    // Compute Jacobian
                    m_jac[0] = m_deriv[0][0][0]*m_deriv[1][1][0]
                                            - m_deriv[1][0][0]*m_deriv[0][1][0];

                    if(CheckJacPositive)
                    {
                        ASSERTL1(m_jac[0] > 0, "2D Regular Jacobian is not positive");
                    }
                    // Spencer's book page 160
                    // Compute derivatives of standard coordinate with respect
                    // to local coordinates.
                    m_gmat[0][0] =  m_deriv[1][1][0]/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -m_deriv[0][1][0]/m_jac[0]; // d xi_2/d x_1
                    m_gmat[2][0] = -m_deriv[1][0][0]/m_jac[0]; // d xi_1/d x_2
                    m_gmat[3][0] =  m_deriv[0][0][0]/m_jac[0]; // d xi_2/d x_2
                }
                else
                {
                    NekDouble g[3];
                    g[0] = m_deriv[0][1][0]*m_deriv[1][2][0] - m_deriv[0][2][0]*m_deriv[1][1][0];
                    g[1] = m_deriv[0][2][0]*m_deriv[1][0][0] - m_deriv[0][0][0]*m_deriv[1][2][0];
                    g[2] = m_deriv[0][0][0]*m_deriv[1][1][0] - m_deriv[0][1][0]*m_deriv[1][0][0];

                    m_jac[0] = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
                    if(CheckJacPositive)
                    {
                        ASSERTL1(m_jac[0] > 0, "Regular Jacobian is not positive");
                    }

                    // d xi_1/d x_1
                    m_gmat[0][0] =  (m_deriv[1][1][0]*g[2] - m_deriv[1][2][0]*g[1])/m_jac[0];
                    // d xi_2/d x_1
                    m_gmat[1][0] = -(m_deriv[0][1][0]*g[2] - m_deriv[0][2][0]*g[1])/m_jac[0];
                    // d_xi_1/d x_2
                    m_gmat[2][0] = -(m_deriv[1][0][0]*g[2] - m_deriv[1][2][0]*g[0])/m_jac[0];
                    // d_xi_2/d x_2
                    m_gmat[3][0] =  (m_deriv[0][0][0]*g[2] - m_deriv[0][2][0]*g[0])/m_jac[0];
                    // d_xi_1/d x_3
                    m_gmat[4][0] =  (m_deriv[1][0][0]*g[1] - m_deriv[1][1][0]*g[0])/m_jac[0];
                    // d xi_2/d x_3
                    m_gmat[5][0] = -(m_deriv[0][0][0]*g[1] - m_deriv[0][1][0]*g[0])/m_jac[0];

                    m_jac[0] = sqrt(m_jac[0]);
                }
            }
            else
            {
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(2*m_coordDim,nqtot,0.0);

                if(m_coordDim == 2) // assume g = [0,0,1]
                {
                    // set up Jacobian
                    Vmath::Vmul (nqtot, &m_deriv[1][0][0], 1, &m_deriv[0][1][0], 1,
                                 &m_jac[0], 1);
                    Vmath::Vvtvm(nqtot, &m_deriv[0][0][0], 1, &m_deriv[1][1][0], 1,
                                 &m_jac[0], 1, &m_jac[0], 1);

                    if(CheckJacPositive)
                    {
                        static int cnt = 0;

                        ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0,
                                 "2D Deformed Jacobian is not positive (cnt = " + 
                                 boost::lexical_cast<std::string>(cnt) + ")");
                        cnt++;
                    }

                    // d xi_1/d x_1
                    Vmath::Vdiv(nqtot,&m_deriv[1][1][0],1,&m_jac[0],1,&m_gmat[0][0],1);
                    Vmath::Vdiv(nqtot,&m_deriv[0][1][0],1,&m_jac[0],1,&m_gmat[1][0],1);
                    // d xi_2/d x_1
                    Vmath::Neg(nqtot,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&m_deriv[1][0][0],1,&m_jac[0],1,&m_gmat[2][0],1);
                    // d xi_1/d x_2
                    Vmath::Neg(nqtot,&m_gmat[2][0],1);
                    // d xi_2/d x_2
                    Vmath::Vdiv(nqtot,&m_deriv[0][0][0],1,&m_jac[0],1,&m_gmat[3][0],1);
                }
                else
                {
                    Array<OneD,NekDouble> g[3] = {Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot)
                                                  };
                    // g[0]
                    Vmath::Vmul (nqtot, &m_deriv[0][2][0], 1, &m_deriv[1][1][0], 1, &g[0][0],1);
                    Vmath::Vvtvm(nqtot, &m_deriv[0][1][0], 1, &m_deriv[1][2][0], 1, &g[0][0],1,
                                 &g[0][0],1);
                    //g[1]
                    Vmath::Vmul (nqtot, &m_deriv[0][0][0], 1, &m_deriv[1][2][0], 1, &g[1][0],1);
                    Vmath::Vvtvm(nqtot, &m_deriv[0][2][0], 1, &m_deriv[1][0][0], 1, &g[1][0],1,
                                 &g[1][0],1);
                    //g[2]
                    Vmath::Vmul (nqtot, &m_deriv[0][1][0], 1, &m_deriv[1][0][0], 1, &g[2][0],1);
                    Vmath::Vvtvm(nqtot, &m_deriv[0][0][0], 1, &m_deriv[1][1][0], 1, &g[2][0],1,
                                 &g[2][0],1);

                    // J_3D
                    Vmath::Vmul (nqtot, &g[0][0], 1, &g[0][0], 1, &m_jac[0], 1);
                    Vmath::Vvtvp(nqtot, &g[1][0], 1, &g[1][0], 1, &m_jac[0], 1,
                                 &m_jac[0],1);
                    Vmath::Vvtvp(nqtot, &g[2][0], 1, &g[2][0],1, &m_jac[0], 1,
                                 &m_jac[0],1);

                    // d xi_1/d x_1
                    Vmath::Vmul (nqtot,&m_deriv[1][2][0],1,&g[1][0],1,&m_gmat[0][0],1);
                    Vmath::Vvtvm(nqtot,&m_deriv[1][1][0],1,&g[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);

                    // d xi_2/d x_1
                    Vmath::Vmul (nqtot,&m_deriv[0][1][0],1,&g[2][0],1,&m_gmat[1][0],1);
                    Vmath::Vvtvm(nqtot,&m_deriv[0][2][0],1,&g[1][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);

                    // d xi_1/d x_2
                    Vmath::Vmul (nqtot,&m_deriv[1][0][0],1,&g[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vvtvm(nqtot,&m_deriv[1][2][0],1,&g[0][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);

                    // d xi_2/d x_2
                    Vmath::Vmul (nqtot,&m_deriv[0][2][0],1,&g[0][0],1,&m_gmat[3][0],1);
                    Vmath::Vvtvm(nqtot,&m_deriv[0][0][0],1,&g[2][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);

                    // d xi_1/d x_3
                    Vmath::Vmul (nqtot,&m_deriv[1][1][0],1,&g[0][0],1,&m_gmat[4][0],1);
                    Vmath::Vvtvm(nqtot,&m_deriv[1][0][0],1,&g[1][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);

                    // d xi_2/d x_3
                    Vmath::Vmul (nqtot,&m_deriv[0][0][0],1,&g[1][0],1,&m_gmat[5][0],1);
                    Vmath::Vvtvm(nqtot,&m_deriv[0][1][0],1,&g[0][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);

                    // J = sqrt(J_3D)
                    Vmath::Vsqrt(nqtot,&m_jac[0],1,&m_jac[0],1);
                }
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
         * Generate Normal vectors at all quadature points specified to the
         * pointsKey "to_key" according to anticlockwise convention.
         * @param   shape       Expansion shape.
         * @param   edge        Edge to evaluate normals on.
         * @param   to_key      PointsKey describing quadrature points.
         * @returns Array of normals organised with consecutive concatenated
         *          \f$N_Q\f$ blocks for each coordinate.
         */
/*        void GeomFactors2D::v_ComputeEdgeNormals(
                            const int edge,
                            const LibUtilities::PointsKey &to_key,
                            Array<OneD, Array<OneD, NekDouble> > &returnval) const
        {
            int i;
            StdRegions::ExpansionType shape = m_coords[0]->DetExpansionType();
            int nqe = to_key.GetNumPoints();

            // Regular geometry case
            if((m_type == eRegular)||(m_type == eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch(m_coords[0]->DetExpansionType())
                {
                case StdRegions::eTriangle:
                    switch(edge)
                    {
                    case 0:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i+1][0],returnval[i],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i+1][0] + m_gmat[2*i][0],returnval[i],1);
                        }
                            break;
                    case 2:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i][0],returnval[i],1);
                        }
                        break;
                    default:
                        ASSERTL0(false,"Edge is out of range (edge < 3)");
                    }
                    break;
                case StdRegions::eQuadrilateral:
                    switch(edge)
                    {
                    case 0:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i+1][0],returnval[i],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i][0],returnval[i],1);
                        }
                        break;
                    case 2:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i+1][0],returnval[i],1);
                        }
                        break;
                    case 3:
                        for(i = 0; i < m_coordDim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i][0],returnval[i],1);
                        }
                        break;
                    default:
                        ASSERTL0(false,"edge is out of range (edge < 4)");
                    }
                    break;
                }

                // normalise
                fac = 0.0;
                for(i =0 ; i < m_coordDim; ++i)
                {
                    fac += returnval[i][0]*returnval[i][0];
                }
                fac = 1.0/sqrt(fac);
                for (i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Smul(nqe,fac,returnval[i],1,returnval[i],1);
                }
            }
            else   // Set up deformed normals
            {
                int j;

                int nquad0 = m_pointsKey[0].GetNumPoints();
                int nquad1 = m_pointsKey[1].GetNumPoints();

                LibUtilities::PointsKey from_key;

                Array<OneD,NekDouble> normals(m_coordDim*max(nquad0,nquad1),0.0);
                Array<OneD,NekDouble> jac    (m_coordDim*max(nquad0,nquad1),0.0);

                // Extract Jacobian along edges and recover local
                // derivates (dx/dr) for polynomial interpolation by
                // multiplying m_gmat by jacobian
                switch(shape)
                {
                case StdRegions::eTriangle:
                    {
                        switch(edge)
                        {
                        case 0:
                            for(j = 0; j < nquad0; ++j)
                            {
                                jac[j] = m_jac[j];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                    normals[i*nquad0+j] = -m_gmat[2*i+1][j]*jac[j];
                                }
                            }
                            from_key = m_pointsKey[0];
                            break;
                        case 1:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j+nquad0-1];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                    normals[i*nquad1+j] = (m_gmat[2*i][nquad0*j + nquad0-1] +  m_gmat[2*i+1][nquad0*j + nquad0-1])*jac[j];
                                }
                            }
                            from_key = m_pointsKey[1];
                            break;
                        case 2:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                    normals[i*nquad1+j] = -m_gmat[2*i][nquad0*j]*jac[j];
                                }
                            }
                            from_key = m_pointsKey[1];
                            break;
                        default:
                            ASSERTL0(false,"edge is out of range (edge < 3)");

                        }
                    }
                    break;
                case StdRegions::eQuadrilateral:
                    {
                        switch(edge)
                        {
                        case 0:
                            for(j = 0; j < nquad0; ++j)
                            {
                                jac[j] = m_jac[j];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                   normals[i*nquad0+j] = -m_gmat[2*i+1][j]*jac[j];
                                }
                           }
                            from_key = m_pointsKey[0];
                            break;
                        case 1:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j+nquad0-1];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                    normals[i*nquad1+j]  = m_gmat[2*i][nquad0*j + nquad0-1]*jac[j];
                                }
                            }
                            from_key = m_pointsKey[1];
                            break;
                        case 2:
                            for(j = 0; j < nquad0; ++j)
                            {
                                jac[j] = m_jac[nquad0*(nquad1-1)+j];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                    normals[i*nquad0+j] = (m_gmat[2*i+1][nquad0*(nquad1-1)+j])*jac[j];
                                }
                            }
                            from_key = m_pointsKey[0];
                            break;
                        case 3:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j];
                                for(i = 0; i < m_coordDim; ++i)
                                {
                                    normals[i*nquad1+j] = -m_gmat[2*i][nquad0*j]*jac[j];
                                }
                            }
                            from_key = m_pointsKey[1];
                            break;
                        default:
                            ASSERTL0(false,"edge is out of range (edge < 3)");
                        }
                    }
                    break;
                default:
                    break;
                }


                int nq  = from_key.GetNumPoints();
                Array<OneD,NekDouble> work(nqe,0.0);

                // interpolate Jacobian and invert
                LibUtilities::Interp1D(from_key,jac,to_key,work);
                Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                // interpolate
                for(i = 0; i < m_coordDim; ++i)
                {
                    LibUtilities::Interp1D(from_key,&normals[i*nq],to_key,&returnval[i][0]);
                    Vmath::Vmul(nqe,work,1,returnval[i],1,returnval[i],1);
                }

                //normalise normal vectors
                Vmath::Zero(nqe,work,1);
                for(i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vvtvp(nqe,returnval[i],1, returnval[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nqe,work,1,work,1);
                Vmath::Sdiv(nqe,1.0,work,1,work,1);

                for(i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vmul(nqe,returnval[i],1,work,1,returnval[i],1);
                }

                // Reverse direction so that points are in
                // anticlockwise direction if edge >=2
                if(edge >= 2)
                {
                    for(i = 0; i < m_coordDim; ++i)
                    {
                        Vmath::Reverse(nqe,returnval[i],1, returnval[i],1);
                    }
                }
            }            
        }
*/

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

        /**
         * @brief Set up the m_weightedjac array, which holds the Jacobian at
         * each quadrature point multipled by the quadrature weight.
         */
        void GeomFactors2D::v_SetUpQuadratureMetrics(
             LibUtilities::ShapeType                        shape,
            const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == m_expDim,
                     "Inappropriate dimension of tbasis");

            int i;
            int nquad0 = m_pointsKey[0].GetNumPoints();
            int nquad1 = m_pointsKey[1].GetNumPoints();
            int nqtot  = nquad0*nquad1;

            m_weightedjac        = Array<OneD, NekDouble>(nqtot);
            m_isUsingQuadMetrics = true;

            // Fill the array m_weighted jac with the values of the (already
            // computed) jacobian (=m_jac)
            if (m_type == eRegular || m_type == eMovingRegular)
            {
                Vmath::Fill(nqtot,m_jac[0],m_weightedjac.get(),1);
            }
            else
            {
                Vmath::Vcopy(nqtot,m_jac.get(),1,m_weightedjac.get(),1);
            }

            // Get hold of the quadrature weights
            const Array<OneD, const NekDouble>& w0 = tbasis[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = tbasis[1]->GetW();

            // Multiply the jacobian with the quadrature weights
            switch(shape)
            {
            case LibUtilities::eQuadrilateral:
                {
                    for (i = 0; i < nquad1; ++i)
                    {
                        Vmath::Vmul(nquad0,m_weightedjac.get()+i*nquad0,1,
                                    w0.get(),1,m_weightedjac.get()+i*nquad0,1);
                    }

                    for (i = 0; i < nquad0; ++i)
                    {
                        Vmath::Vmul(nquad1, m_weightedjac.get()+i, nquad0,
                                    w1.get(), 1, m_weightedjac.get()+i, nquad0);
                    }

                    break;
                }

            case LibUtilities::eTriangle:
                {
                    for(i = 0; i < nquad1; ++i)
                    {
                        Vmath::Vmul(nquad0,m_weightedjac.get()+i*nquad0,1,
                                    w0.get(),1,m_weightedjac.get()+i*nquad0,1);
                    }
                    
                    switch(tbasis[1]->GetPointsType())
                    {
                        case LibUtilities::ePolyEvenlySpaced:
                    case LibUtilities::eGaussLobattoLegendre:
                        {
                            const Array<OneD, const NekDouble>& z1 =
                                tbasis[1]->GetZ();
                            for (i = 0; i < nquad1; ++i)
                            {
                                Blas::Dscal(nquad0, 0.5*(1-z1[i])*w1[i],
                                            m_weightedjac.get()+i*nquad0, 1);
                            }
                            break;
                        }
                    case LibUtilities::eGaussRadauMAlpha1Beta0:
                        {
                            for (i = 0; i < nquad1; ++i)
                            {
                                Blas::Dscal(nquad0, 0.5*w1[i],
                                            m_weightedjac.get()+i*nquad0, 1);
                            }
                            break;
                        }
                    default:
                        {
                            m_isUsingQuadMetrics = false;
                            m_weightedjac = Array<OneD, NekDouble>();
                            return;
                        }
                    }
                    break;
                }
            default:
                {
                    ASSERTL0(false,"Invalid shape type");
                }
            }
        }

        /**
         *
         */
        void GeomFactors2D::v_SetUpLaplacianMetrics(LibUtilities::ShapeType shape,
                                                  const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == m_expDim,"Inappropriate dimension of tbasis");
            ASSERTL1((m_coordDim == 2)||(m_coordDim <= 3),
                     "The coordinate dimension should be equal to two or three"
                     "for two-dimensional elements");

            int i;
            int nquad0 = m_pointsKey[0].GetNumPoints();
            int nquad1 = m_pointsKey[1].GetNumPoints();
            int nqtot  = nquad0*nquad1;

            m_laplacianmetrics      = Array<TwoD, NekDouble>(3,nqtot);
            m_laplacianMetricIsZero = Array<OneD, bool>(3, false);
            m_isUsingLaplMetrics  = true;

            switch(shape)
            {
            case LibUtilities::eQuadrilateral:
                {
                    if(( m_type == eRegular)||
                       ( m_type == eMovingRegular))
                    {
                        NekDouble g0 = m_gmat[0][0]*m_gmat[0][0] + m_gmat[2][0]*m_gmat[2][0];
                        NekDouble g1 = m_gmat[0][0]*m_gmat[1][0] + m_gmat[2][0]*m_gmat[3][0];
                        NekDouble g2 = m_gmat[1][0]*m_gmat[1][0] + m_gmat[3][0]*m_gmat[3][0];

                        if(m_coordDim == 3)
                        {
                            g0 += m_gmat[4][0]*m_gmat[4][0];
                            g1 += m_gmat[4][0]*m_gmat[5][0];
                            g2 += m_gmat[5][0]*m_gmat[5][0];
                        }

                        if(fabs(g1) < NekConstants::kGeomFactorsTol)
                        {
                            m_laplacianMetricIsZero[1] = true;
                        }

                        Vmath::Fill(nqtot,g0,&m_laplacianmetrics[0][0],1);
                        Vmath::Fill(nqtot,g1,&m_laplacianmetrics[1][0],1);
                        Vmath::Fill(nqtot,g2,&m_laplacianmetrics[2][0],1);
                    }
                    else
                    {
                        Vmath::Vmul (nqtot,&m_gmat[0][0],1,&m_gmat[0][0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[2][0],1,&m_gmat[2][0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);

                        Vmath::Vmul (nqtot,&m_gmat[0][0],1,&m_gmat[1][0],1,&m_laplacianmetrics[1][0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[2][0],1,&m_gmat[3][0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);

                        Vmath::Vmul (nqtot,&m_gmat[1][0],1,&m_gmat[1][0],1,&m_laplacianmetrics[2][0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[3][0],1,&m_gmat[3][0],1,&m_laplacianmetrics[2][0],1,&m_laplacianmetrics[2][0],1);

                        if(m_coordDim == 3)
                        {
                            Vmath::Vvtvp(nqtot,&m_gmat[4][0],1,&m_gmat[4][0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                            Vmath::Vvtvp(nqtot,&m_gmat[4][0],1,&m_gmat[5][0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                            Vmath::Vvtvp(nqtot,&m_gmat[5][0],1,&m_gmat[5][0],1,&m_laplacianmetrics[2][0],1,&m_laplacianmetrics[2][0],1);
                        }
                    }

                    Vmath::Vmul(nqtot,&m_weightedjac[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                    Vmath::Vmul(nqtot,&m_weightedjac[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                    Vmath::Vmul(nqtot,&m_weightedjac[0],1,&m_laplacianmetrics[2][0],1,&m_laplacianmetrics[2][0],1);
                }
                break;
            case LibUtilities::eTriangle:
                {
                    Array<OneD, NekDouble> dEta_dXi[2] = {Array<OneD, NekDouble>(nqtot,1.0),
                                                          Array<OneD, NekDouble>(nqtot,1.0)};

                    const Array<OneD, const NekDouble>& z0 = tbasis[0]->GetZ();
                    const Array<OneD, const NekDouble>& z1 = tbasis[1]->GetZ();

                    for(i = 0; i < nquad1; i++)
                    {
                        Blas::Dscal(nquad0,2.0/(1-z1[i]),&dEta_dXi[0][0]+i*nquad0,1);
                        Blas::Dscal(nquad0,2.0/(1-z1[i]),&dEta_dXi[1][0]+i*nquad0,1);
                    }
                    for(i = 0; i < nquad0; i++)
                    {
                        Blas::Dscal(nquad1,0.5*(1+z0[i]),&dEta_dXi[1][0]+i,nquad0);
                    }

                    Array<OneD, NekDouble> tmp(nqtot);

                    if(( m_type == eRegular)||
                       ( m_type == eMovingRegular))
                    {
                        Vmath::Smul (nqtot,m_gmat[0][0],&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Svtvp(nqtot,m_gmat[1][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vmul (nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Smul (nqtot,m_gmat[1][0],&tmp[0],1,&m_laplacianmetrics[1][0],1);


                        Vmath::Smul (nqtot,m_gmat[2][0],&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Svtvp(nqtot,m_gmat[3][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vvtvp(nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Svtvp(nqtot,m_gmat[3][0],&tmp[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);

                        if(m_coordDim == 3)
                        {
                            Vmath::Smul (nqtot,m_gmat[4][0],&dEta_dXi[0][0],1,&tmp[0],1);
                            Vmath::Svtvp(nqtot,m_gmat[5][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                            Vmath::Vvtvp(nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                            Vmath::Svtvp(nqtot,m_gmat[5][0],&tmp[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                        }

                        NekDouble g2 = m_gmat[1][0]*m_gmat[1][0] + m_gmat[3][0]*m_gmat[3][0];
                        if(m_coordDim == 3)
                        {
                            g2 += m_gmat[5][0]*m_gmat[5][0];
                        }
                        Vmath::Fill(nqtot,g2,&m_laplacianmetrics[2][0],1);


                    }
                    else
                    {
                        Vmath::Vmul (nqtot,&m_gmat[0][0],1,&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[1][0],1,&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vmul (nqtot,&tmp[0],      1,&tmp[0],      1,&m_laplacianmetrics[0][0],1);
                        Vmath::Vmul (nqtot,&m_gmat[1][0],1,&tmp[0],      1,&m_laplacianmetrics[1][0],1);
                        Vmath::Vmul (nqtot,&m_gmat[1][0],1,&m_gmat[1][0],1,&m_laplacianmetrics[2][0],1);


                        Vmath::Vmul (nqtot,&m_gmat[2][0],1,&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[3][0],1,&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vvtvp(nqtot,&tmp[0],1,&tmp[0],            1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[3][0],1,&tmp[0],      1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                        Vmath::Vvtvp(nqtot,&m_gmat[3][0],1,&m_gmat[3][0],1,&m_laplacianmetrics[2][0],1,&m_laplacianmetrics[2][0],1);

                        if(m_coordDim == 3)
                        {
                            Vmath::Vmul (nqtot,&m_gmat[4][0],1,&dEta_dXi[0][0],1,&tmp[0],1);
                            Vmath::Vvtvp(nqtot,&m_gmat[5][0],1,&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                            Vmath::Vvtvp(nqtot,&tmp[0],1,&tmp[0],            1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                            Vmath::Vvtvp(nqtot,&m_gmat[5][0],1,&tmp[0],      1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                            Vmath::Vvtvp(nqtot,&m_gmat[5][0],1,&m_gmat[5][0],1,&m_laplacianmetrics[2][0],1,&m_laplacianmetrics[2][0],1);
                        }
                    }

                    Vmath::Vmul(nqtot,&m_weightedjac[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                    Vmath::Vmul(nqtot,&m_weightedjac[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                    Vmath::Vmul(nqtot,&m_weightedjac[0],1,&m_laplacianmetrics[2][0],1,&m_laplacianmetrics[2][0],1);
                }
                break;
            default:
                {
                    ASSERTL0(false,"Invalid shape type");
                }
            }
        }

    }
}
