#include <SpatialDomains/GeomFactors2D.h>

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
                        const bool LaplMetrics) :
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
            for (int i = 0; i < mCoordDim; ++i)
            {
                mCoords[i] = Coords[i];
            }

            StdRegions::ExpansionType shape = Coords[0]->DetExpansionType();

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
            mPointsKey[0] = pkey0_tbasis;
            mPointsKey[1] = pkey1_tbasis;

            // setup temp storage
            Array<OneD, Array<OneD,NekDouble> > d1_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d2_map   (coordim);

            mDeriv = Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);
            mDeriv[0] = Array<OneD, Array<OneD,NekDouble> >(mCoordDim);
            mDeriv[1] = Array<OneD, Array<OneD,NekDouble> >(mCoordDim);

            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                d1_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d2_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                mDeriv[0][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                mDeriv[1][i] = Array<OneD,NekDouble>(nqtot_tbasis);

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
                    mDeriv[0][i] = d1_map[i];
                    mDeriv[1][i] = d2_map[i];
                }
                else
                {
                    LibUtilities::Interp2D(pkey0_map, pkey1_map, d1_map[i],
                                    pkey0_tbasis, pkey1_tbasis, mDeriv[0][i]);
                    LibUtilities::Interp2D(pkey0_map, pkey1_map, d2_map[i],
                                    pkey0_tbasis, pkey1_tbasis, mDeriv[1][i]);
                }
            }

            // Setting up Surface Normal Vectors
            if(coordim == 3)
            {
                v_ComputeSurfaceNormals();
            }

            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation
            // metrics
            SetUpJacGmat2D();

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
        void GeomFactors2D::SetUpJacGmat2D()
        {
            ASSERTL1(mDeriv[0].num_elements()==mCoordDim,
                     "The dimension of array d1 does not match the coordinate "
                     "dimension");
            ASSERTL1(mDeriv[1].num_elements()==mCoordDim,
                     "The dimension of array d2 does not match the coordinate "
                     "dimension");

            int nqtot = mPointsKey[0].GetNumPoints() *
                        mPointsKey[1].GetNumPoints();

            ASSERTL1(mDeriv[0][0].num_elements() == nqtot,
                     "Number of quadrature points do not match");
            ASSERTL1(mDeriv[1][0].num_elements() == nqtot,
                     "Number of quadrature points do not match");

            if((mType == eRegular)||(mType == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(2*mCoordDim,1,0.0);

                if(mCoordDim == 2) // assume g = [0,0,1]
                {
                    m_jac[0] = mDeriv[0][0][0]*mDeriv[1][1][0] - mDeriv[1][0][0]*mDeriv[0][1][0];

                    ASSERTL1(m_jac[0] > 0,
                             "2D Regular Jacobian is not positive");

                    // Spencer's book page 160
                    m_gmat[0][0] =  mDeriv[1][1][0]/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -mDeriv[0][1][0]/m_jac[0]; // d xi_2/d x_1
                    m_gmat[2][0] = -mDeriv[1][0][0]/m_jac[0]; // d xi_1/d x_2
                    m_gmat[3][0] =  mDeriv[0][0][0]/m_jac[0]; // d xi_2/d x_2
                }
                else
                {
                    NekDouble g[3];
                    g[0] = mDeriv[0][1][0]*mDeriv[1][2][0] - mDeriv[0][2][0]*mDeriv[1][1][0];
                    g[1] = mDeriv[0][2][0]*mDeriv[1][0][0] - mDeriv[0][0][0]*mDeriv[1][2][0];
                    g[2] = mDeriv[0][0][0]*mDeriv[1][1][0] - mDeriv[0][1][0]*mDeriv[1][0][0];

                    m_jac[0] = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
                    ASSERTL1(m_jac[0] > 0, "Regular Jacobian is not positive");

                    // d xi_1/d x_1
                    m_gmat[0][0] =  (mDeriv[1][1][0]*g[2] - mDeriv[1][2][0]*g[1])/m_jac[0];
                    // d xi_2/d x_1
                    m_gmat[1][0] = -(mDeriv[0][1][0]*g[2] - mDeriv[0][2][0]*g[1])/m_jac[0];
                    // d_xi_1/d x_2
                    m_gmat[2][0] = -(mDeriv[1][0][0]*g[2] - mDeriv[1][2][0]*g[0])/m_jac[0];
                    // d_xi_2/d x_2
                    m_gmat[3][0] =  (mDeriv[0][0][0]*g[2] - mDeriv[0][2][0]*g[0])/m_jac[0];
                    // d_xi_1/d x_3
                    m_gmat[4][0] =  (mDeriv[1][0][0]*g[1] - mDeriv[1][1][0]*g[0])/m_jac[0];
                    // d xi_2/d x_3
                    m_gmat[5][0] = -(mDeriv[0][0][0]*g[1] - mDeriv[0][1][0]*g[0])/m_jac[0];

                    m_jac[0] = sqrt(m_jac[0]);
                }
            }
            else
            {
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(2*mCoordDim,nqtot,0.0);

                if(mCoordDim == 2) // assume g = [0,0,1]
                {
                    // set up Jacobian
                    Vmath::Vmul (nqtot, &mDeriv[1][0][0], 1, &mDeriv[0][1][0], 1,
                                 &m_jac[0], 1);
                    Vmath::Vvtvm(nqtot, &mDeriv[0][0][0], 1, &mDeriv[1][1][0], 1,
                                 &m_jac[0], 1, &m_jac[0], 1);

                    ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0,
                             "2D Deformed Jacobian is not positive");

                    // d xi_1/d x_1
                    Vmath::Vdiv(nqtot,&mDeriv[1][1][0],1,&m_jac[0],1,&m_gmat[0][0],1);
                    Vmath::Vdiv(nqtot,&mDeriv[0][1][0],1,&m_jac[0],1,&m_gmat[1][0],1);
                    // d xi_2/d x_1
                    Vmath::Neg(nqtot,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&mDeriv[1][0][0],1,&m_jac[0],1,&m_gmat[2][0],1);
                    // d xi_1/d x_2
                    Vmath::Neg(nqtot,&m_gmat[2][0],1);
                    // d xi_2/d x_2
                    Vmath::Vdiv(nqtot,&mDeriv[0][0][0],1,&m_jac[0],1,&m_gmat[3][0],1);
                }
                else
                {
                    Array<OneD,NekDouble> g[3] = {Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot)
                                                  };
                    // g[0]
                    Vmath::Vmul (nqtot, &mDeriv[0][2][0], 1, &mDeriv[1][1][0], 1, &g[0][0],1);
                    Vmath::Vvtvm(nqtot, &mDeriv[0][1][0], 1, &mDeriv[1][2][0], 1, &g[0][0],1,
                                 &g[0][0],1);
                    //g[1]
                    Vmath::Vmul (nqtot, &mDeriv[0][0][0], 1, &mDeriv[1][2][0], 1, &g[1][0],1);
                    Vmath::Vvtvm(nqtot, &mDeriv[0][2][0], 1, &mDeriv[1][0][0], 1, &g[1][0],1,
                                 &g[1][0],1);
                    //g[2]
                    Vmath::Vmul (nqtot, &mDeriv[0][1][0], 1, &mDeriv[1][0][0], 1, &g[2][0],1);
                    Vmath::Vvtvm(nqtot, &mDeriv[0][0][0], 1, &mDeriv[1][1][0], 1, &g[2][0],1,
                                 &g[2][0],1);

                    // J_3D
                    Vmath::Vmul (nqtot, &g[0][0], 1, &g[0][0], 1, &m_jac[0], 1);
                    Vmath::Vvtvp(nqtot, &g[1][0], 1, &g[1][0], 1, &m_jac[0], 1,
                                 &m_jac[0],1);
                    Vmath::Vvtvp(nqtot, &g[2][0], 1, &g[2][0],1, &m_jac[0], 1,
                                 &m_jac[0],1);

                    // d xi_1/d x_1
                    Vmath::Vmul (nqtot,&mDeriv[1][2][0],1,&g[1][0],1,&m_gmat[0][0],1);
                    Vmath::Vvtvm(nqtot,&mDeriv[1][1][0],1,&g[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);

                    // d xi_2/d x_1
                    Vmath::Vmul (nqtot,&mDeriv[0][1][0],1,&g[2][0],1,&m_gmat[1][0],1);
                    Vmath::Vvtvm(nqtot,&mDeriv[0][2][0],1,&g[1][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);

                    // d xi_1/d x_2
                    Vmath::Vmul (nqtot,&mDeriv[1][0][0],1,&g[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vvtvm(nqtot,&mDeriv[1][2][0],1,&g[0][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);

                    // d xi_2/d x_2
                    Vmath::Vmul (nqtot,&mDeriv[0][2][0],1,&g[0][0],1,&m_gmat[3][0],1);
                    Vmath::Vvtvm(nqtot,&mDeriv[0][0][0],1,&g[2][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);

                    // d xi_1/d x_3
                    Vmath::Vmul (nqtot,&mDeriv[1][1][0],1,&g[0][0],1,&m_gmat[4][0],1);
                    Vmath::Vvtvm(nqtot,&mDeriv[1][0][0],1,&g[1][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);

                    // d xi_2/d x_3
                    Vmath::Vmul (nqtot,&mDeriv[0][0][0],1,&g[1][0],1,&m_gmat[5][0],1);
                    Vmath::Vvtvm(nqtot,&mDeriv[0][1][0],1,&g[0][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
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
            int nq = mPointsKey[0].GetNumPoints() *
                        mPointsKey[1].GetNumPoints();

            output = Array<OneD,Array<OneD,NekDouble> >(mCoordDim);
            for (int i = 0; i < mCoordDim; ++i)
            {
                output[i] = Array<OneD, NekDouble> (nq, 0.0);
            }

            // Construction of Connection
            switch(mTangentDir)
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
                    Array<OneD,NekDouble> x0(nq);
                    Array<OneD,NekDouble> x1(nq);
                    Array<OneD,NekDouble> x2(nq);
                    mCoords[0]->GetCoords(x0,x1,x2);

                    // circular around the center of the domain
                    NekDouble radius,xmax, xmin, xmid, ymax, ymin, ymid, xdis, ydis;

                    xmax = Vmath::Vmax(nq,x0,1);
                    xmin = Vmath::Vmin(nq,x0,1);
                    ymax = Vmath::Vmax(nq,x1,1);
                    ymin = Vmath::Vmin(nq,x1,1);

                    xmid = 0.5*(xmax+xmin);
                    ymid = 0.5*(ymax+ymin);

                    for (int i = 0; i < nq; i++)
                    {
                        xdis = x0[i]-xmid;
                        ydis = x1[i]-ymid;

                        radius = sqrt(xdis*xdis+ydis*ydis);
                        output[0][i] = ydis/radius;
                        output[1][i] = -1.0*xdis/radius;
                    }
                    break;
                }
                case eTangentCircular2:
                {
                    cout << "Anisotropy layers" << endl;
                    Array<OneD,NekDouble> x0(nq);
                    Array<OneD,NekDouble> x1(nq);
                    Array<OneD,NekDouble> x2(nq);
                    mCoords[0]->GetCoords(x0,x1,x2);


                    // circular around the center of the domain
                    NekDouble radius, xc, yc, xdis, ydis;

                    xc = 25.0;
                    yc = 0.0;

                    for (int i = 0; i < nq; i++)
                    {

                        if(x0[i]<=xc)
                        {
                            output[0][i] = 1.0;
                            output[1][i] = 0.0;
                        }

                        else
                        {
                            xdis = x0[i]-xc;
                            ydis = x1[i]-yc;
                            radius = sqrt(xdis*xdis+ydis*ydis);
                            output[0][i] = ydis/radius;
                            output[1][i] = -1.0*xdis/radius;
                        }
                    }
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
            int coordim = mDeriv[0].num_elements();
            // Total number of points in the element.
            int nqtot = mDeriv[0][0].num_elements();

            // Allocate temporary storage.
            Array<OneD, NekDouble> temp(nqtot,0.0);

            // Initialization of storage for normals
            mNormal = Array<OneD, Array<OneD, NekDouble> >(mCoordDim);
            for (int i = 0; i < mCoordDim; ++i)
            {
                mNormal[i] = Array<OneD,NekDouble>(nqtot);
            }

            // Derive Geometric Normal vectors by cross product of two
            // tangential basis vectors.
            Vmath::Vmul (nqtot, mDeriv[0][2],     1,
                                mDeriv[1][1],     1,
                                temp,             1);
            Vmath::Vvtvm(nqtot, mDeriv[0][1],     1,
                                mDeriv[1][2],     1,
                                temp,             1,
                                mNormal[0],          1);

            Vmath::Vmul (nqtot, mDeriv[0][0],     1,
                                mDeriv[1][2],     1,
                                temp,           1);
            Vmath::Vvtvm(nqtot, mDeriv[0][2],     1,
                                mDeriv[1][0],     1,
                                temp,           1,
                                mNormal[1],    1);

            Vmath::Vmul (nqtot, mDeriv[0][1],     1,
                                mDeriv[1][0],     1,
                                temp,           1);
            Vmath::Vvtvm(nqtot, mDeriv[0][0],     1,
                                mDeriv[1][1],     1,
                                temp,           1,
                                mNormal[2],  1);

            // Reset temp array
            temp = Array<OneD, NekDouble>(nqtot,0.0);

            // Normalization of Surface Normal
            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vvtvp(nqtot, mNormal[i],  1,
                                    mNormal[i],  1,
                                    temp,           1,
                                    temp,           1);
            }
            Vmath::Vsqrt(nqtot, temp, 1, temp, 1);

            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vdiv(nqtot,  mNormal[i],  1,
                                    temp,           1,
                                    mNormal[i],  1);
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
        void GeomFactors2D::v_ComputeEdgeNormals(
                            const int edge,
                            const LibUtilities::PointsKey &to_key,
                            Array<OneD, Array<OneD, NekDouble> > &returnval) const
        {
            int i;
            StdRegions::ExpansionType shape = mCoords[0]->DetExpansionType();
            int nqe = to_key.GetNumPoints();

            // Regular geometry case
            if((mType == eRegular)||(mType == eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch(mCoords[0]->DetExpansionType())
                {
                case StdRegions::eTriangle:
                    switch(edge)
                    {
                    case 0:
                        for(i = 0; i < mCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i+1][0],returnval[i],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < mCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i+1][0] + m_gmat[2*i][0],returnval[i],1);
                        }
                            break;
                    case 2:
                        for(i = 0; i < mCoordDim; ++i)
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
                        for(i = 0; i < mCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i+1][0],returnval[i],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < mCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i][0],returnval[i],1);
                        }
                        break;
                    case 2:
                        for(i = 0; i < mCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i+1][0],returnval[i],1);
                        }
                        break;
                    case 3:
                        for(i = 0; i < mCoordDim; ++i)
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
                for(i =0 ; i < mCoordDim; ++i)
                {
                    fac += returnval[i][0]*returnval[i][0];
                }
                fac = 1.0/sqrt(fac);
                for (i = 0; i < mCoordDim; ++i)
                {
                    Vmath::Smul(nqe,fac,returnval[i],1,returnval[i],1);
                }
            }
            else   // Set up deformed normals
            {
                int j;

                int nquad0 = mPointsKey[0].GetNumPoints();
                int nquad1 = mPointsKey[1].GetNumPoints();

                LibUtilities::PointsKey from_key;

                Array<OneD,NekDouble> normals(mCoordDim*max(nquad0,nquad1),0.0);
                Array<OneD,NekDouble> jac    (mCoordDim*max(nquad0,nquad1),0.0);

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
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                    normals[i*nquad0+j] = -m_gmat[2*i+1][j]*jac[j];
                                }
                            }
                            from_key = mPointsKey[0];
                            break;
                        case 1:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j+nquad0-1];
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                    normals[i*nquad1+j] = (m_gmat[2*i][nquad0*j + nquad0-1] +  m_gmat[2*i+1][nquad0*j + nquad0-1])*jac[j];
                                }
                            }
                            from_key = mPointsKey[1];
                            break;
                        case 2:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j];
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                    normals[i*nquad1+j] = -m_gmat[2*i][nquad0*j]*jac[j];
                                }
                            }
                            from_key = mPointsKey[1];
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
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                   normals[i*nquad0+j] = -m_gmat[2*i+1][j]*jac[j];
                                }
                           }
                            from_key = mPointsKey[0];
                            break;
                        case 1:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j+nquad0-1];
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                    normals[i*nquad1+j]  = m_gmat[2*i][nquad0*j + nquad0-1]*jac[j];
                                }
                            }
                            from_key = mPointsKey[1];
                            break;
                        case 2:
                            for(j = 0; j < nquad0; ++j)
                            {
                                jac[j] = m_jac[nquad0*(nquad1-1)+j];
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                    normals[i*nquad0+j] = (m_gmat[2*i+1][nquad0*(nquad1-1)+j])*jac[j];
                                }
                            }
                            from_key = mPointsKey[0];
                            break;
                        case 3:
                            for(j = 0; j < nquad1; ++j)
                            {
                                jac[j] = m_jac[nquad0*j];
                                for(i = 0; i < mCoordDim; ++i)
                                {
                                    normals[i*nquad1+j] = -m_gmat[2*i][nquad0*j]*jac[j];
                                }
                            }
                            from_key = mPointsKey[1];
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
                for(i = 0; i < mCoordDim; ++i)
                {
                    LibUtilities::Interp1D(from_key,&normals[i*nq],to_key,&returnval[i][0]);
                    Vmath::Vmul(nqe,work,1,returnval[i],1,returnval[i],1);
                }

                //normalise normal vectors
                Vmath::Zero(nqe,work,1);
                for(i = 0; i < mCoordDim; ++i)
                {
                    Vmath::Vvtvp(nqe,returnval[i],1, returnval[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nqe,work,1,work,1);
                Vmath::Sdiv(nqe,1.0,work,1,work,1);

                for(i = 0; i < mCoordDim; ++i)
                {
                    Vmath::Vmul(nqe,returnval[i],1,work,1,returnval[i],1);
                }

                // Reverse direction so that points are in
                // anticlockwise direction if edge >=2
                if(edge >= 2)
                {
                    for(i = 0; i < mCoordDim; ++i)
                    {
                        Vmath::Reverse(nqe,returnval[i],1, returnval[i],1);
                    }
                }
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
            int nq = mPointsKey[0].GetNumPoints() *
                        mPointsKey[1].GetNumPoints();

            if (mNormal.num_elements() == 0) {
                v_ComputeSurfaceNormals();
            }

            // Get a principle direction for the expansion.
            Array<OneD, Array<OneD,NekDouble> > PrincipleDir;
            ComputePrincipleDirection(PrincipleDir);

            // Allocate temporary and tangent storage.
            Array<OneD, NekDouble> temp(nq, 0.0);
            mTangents  =  Array<OneD, Array<OneD, Array<OneD,NekDouble> > >(2);
            for (i = 0; i < 2; ++i)
            {
                mTangents[i] = Array<OneD, Array<OneD, NekDouble> >(mCoordDim);
                for (k = 0; k < mCoordDim; ++k)
                {
                    mTangents[i][k] = Array<OneD, NekDouble>(nq, 0.0);
                }
            }

            // Gram-schmidz process to make this principaldirection orthonormal
            // to surface normal vectors. Let u1 = v1 = SurfaceNormal,
            // v2 = Principaldirection, which should be orthogornalised to u2
            // inner12 = < u1, v2 >, norm2 = < u1, u1 > = 1 by default
            // ipro = - < u1, v2 > / < u1, u1 >
            for (i = 0; i < mCoordDim; ++i)
            {
                Vmath::Vvtvp(nq, mNormal[i],        1,
                                 PrincipleDir[i],   1,
                                 temp,              1,
                                 temp,              1);
            }
            Vmath::Neg(nq, temp, 1);

            // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
            for (i = 0; i < mCoordDim; ++i)
            {
                Vmath::Vvtvp(nq, temp,              1,
                                 mNormal[i],        1,
                                 PrincipleDir[i],   1,
                                 mTangents[0][i],   1);
            }

            // Normalise first tangent vectors.
            VectorNormalise(mTangents[0]);

            // The second set of tangential vectors is obtained by cross-
            // producting the surface normal with the first tangential vectors.
            VectorCrossProd(mTangents[0], mNormal, mTangents[1]);

            // Normalise second tangent vectors.
            VectorNormalise(mTangents[1]);
        }


        /**
         *
         */
        void GeomFactors2D::v_SetUpQuadratureMetrics(StdRegions::ExpansionType shape,
                                                   const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == mExpDim,"Inappropriate dimension of tbasis");

            int i;
            int nquad0 = mPointsKey[0].GetNumPoints();
            int nquad1 = mPointsKey[1].GetNumPoints();
            int nqtot  = nquad0*nquad1;

            m_weightedjac           = Array<OneD, NekDouble>(nqtot);
            mIsUsingQuadMetrics = true;

            // Fill the array m_weighted jac with the values
            // of the (already computed) jacobian (=m_jac)
            if((mType == eRegular)||(mType == eMovingRegular))
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
            case StdRegions::eQuadrilateral:
                {
                    for(i = 0; i < nquad1; ++i)
                    {
                        Vmath::Vmul(nquad0,m_weightedjac.get()+i*nquad0,1,
                                    w0.get(),1,m_weightedjac.get()+i*nquad0,1);
                    }

                    for(i = 0; i < nquad0; ++i)
                    {
                        Vmath::Vmul(nquad1,m_weightedjac.get()+i,nquad0,w1.get(),1,
                                    m_weightedjac.get()+i,nquad0);
                    }
                }
                break;
            case StdRegions::eTriangle:
                {
                    for(i = 0; i < nquad1; ++i)
                    {
                        Vmath::Vmul(nquad0,m_weightedjac.get()+i*nquad0,1,
                                    w0.get(),1,m_weightedjac.get()+i*nquad0,1);
                    }

                    switch(tbasis[1]->GetPointsType())
                    {
                    case LibUtilities::ePolyEvenlySpaced:
                    case LibUtilities::eGaussLobattoLegendre:  // Legendre inner product
                        for(i = 0; i < nquad1; ++i)
                        {
                            const Array<OneD, const NekDouble>& z1 = tbasis[1]->GetZ();
                            Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i],m_weightedjac.get()+i*nquad0,1);
                        }
                        break;
                    case LibUtilities::eGaussRadauMAlpha1Beta0: // (1,0) Jacobi Inner product
                        for(i = 0; i < nquad1; ++i)
                        {
                            Blas::Dscal(nquad0,0.5*w1[i],m_weightedjac.get()+i*nquad0,1);
                        }
                        break;
                    default:
                        {
                            ASSERTL0(false,"Currently no implementation for this PointsType");
                        }
                    }
                }
                break;
            default:
                {
                    ASSERTL0(false,"Invalid shape type");
                }
            }

        }

        /**
         *
         */
        void GeomFactors2D::v_SetUpLaplacianMetrics(StdRegions::ExpansionType shape,
                                                  const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == mExpDim,"Inappropriate dimension of tbasis");
            ASSERTL1((mCoordDim == 2)||(mCoordDim <= 3),
                     "The coordinate dimension should be equal to two or three"
                     "for two-dimensional elements");

            int i;
            int nquad0 = mPointsKey[0].GetNumPoints();
            int nquad1 = mPointsKey[1].GetNumPoints();
            int nqtot  = nquad0*nquad1;

            m_laplacianmetrics      = Array<TwoD, NekDouble>(3,nqtot);
            m_laplacianMetricIsZero = Array<OneD, bool>(3, false);
            mIsUsingLaplMetrics  = true;

            // Get hold of the quadrature weights
            const Array<OneD, const NekDouble>& w0 = tbasis[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = tbasis[1]->GetW();

            switch(shape)
            {
            case StdRegions::eQuadrilateral:
                {
                    if(( mType == eRegular)||
                       ( mType == eMovingRegular))
                    {
                        NekDouble g0 = m_gmat[0][0]*m_gmat[0][0] + m_gmat[2][0]*m_gmat[2][0];
                        NekDouble g1 = m_gmat[0][0]*m_gmat[1][0] + m_gmat[2][0]*m_gmat[3][0];
                        NekDouble g2 = m_gmat[1][0]*m_gmat[1][0] + m_gmat[3][0]*m_gmat[3][0];

                        if(mCoordDim == 3)
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

                        if(mCoordDim == 3)
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
            case StdRegions::eTriangle:
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

                    if(( mType == eRegular)||
                       ( mType == eMovingRegular))
                    {
                        Vmath::Smul (nqtot,m_gmat[0][0],&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Svtvp(nqtot,m_gmat[1][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vmul (nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Smul (nqtot,m_gmat[1][0],&tmp[0],1,&m_laplacianmetrics[1][0],1);


                        Vmath::Smul (nqtot,m_gmat[2][0],&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Svtvp(nqtot,m_gmat[3][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vvtvp(nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Svtvp(nqtot,m_gmat[3][0],&tmp[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);

                        if(mCoordDim == 3)
                        {
                            Vmath::Smul (nqtot,m_gmat[4][0],&dEta_dXi[0][0],1,&tmp[0],1);
                            Vmath::Svtvp(nqtot,m_gmat[5][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                            Vmath::Vvtvp(nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                            Vmath::Svtvp(nqtot,m_gmat[5][0],&tmp[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                        }

                        NekDouble g2 = m_gmat[1][0]*m_gmat[1][0] + m_gmat[3][0]*m_gmat[3][0];
                        if(mCoordDim == 3)
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

                        if(mCoordDim == 3)
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
