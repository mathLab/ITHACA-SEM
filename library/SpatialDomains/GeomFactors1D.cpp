#include <SpatialDomains/GeomFactors2D.h>
#include <SpatialDomains/Geometry2D.h>

#include <SpatialDomains/GeomFactors1D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * @class GeomFactors1D
         *
         * Computes and stores geometric factors and pointwise geometry
         * information for a 1D segment.
         */

        /**
         * The argument 'tbasis' contains the information about the quadrature
         * points on which the weighted metric terms should be specified.
         * This constructor
         * - Calculates the local factors \f$ d \xi/d x_i \f$ as a recipricol
         *   derivative \f$ d x_i/d\xi \f$;
         * - Calcuate the jacobian as \f$ J = \sqrt{\frac{d x_1}{d _\xi}.
         *   \frac{d x_1}{d _\xi} + \frac{d x_2}{d _\xi}.\frac{d x_2}{d _\xi}
         *   + \frac{d x_3}{d _\xi}\frac{d x_3}{d _\xi}} \f$
         *
         * @param   gtype       Type of geometry.
         * @param   coordim     Dimension of coordinate system.
         * @param   Coords      ?
         * @param   tbasis      Basis for derivatives
         * @param   SetUpQuadratureMetrics  ?
         * @param   SetUpLaplacianMetrics   ?
         */
        GeomFactors1D::GeomFactors1D(const GeomType gtype,
                        const int coordim,
                        const Array<OneD, const StdRegions
                                            ::StdExpansion1DSharedPtr> &Coords,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis,
                        const bool QuadMetrics,
                        const bool LaplMetrics):
            GeomFactors(gtype,1,coordim,QuadMetrics,LaplMetrics)
        {
            // Perform sanity checks
            ASSERTL1((coordim == 1)||(coordim == 2)||(coordim == 3),
                     "The coordinate dimension should be equal to one, two or "
                     "three for one-dimensional elements");
            ASSERTL1(tbasis.num_elements()==1,
                     "tbasis should be an array of size one");
            ASSERTL1(LaplMetrics?QuadMetrics:true,
                     "SetUpQuadratureMetrics should be true if "
                     "SetUpLaplacianMetrics is true");

            for (int i = 0; i < mCoordDim; ++i)
            {
                mCoords[i] = Coords[i];
            }

            // Get the shape of the expansion
            StdRegions::ExpansionType shape = Coords[0]->DetExpansionType();

            // The quadrature points of the mapping
            // (as specified in Coords)
            LibUtilities::PointsKey pkey_map(
                                        Coords[0]->GetBasis(0)->GetPointsKey());
            int nquad_map = pkey_map.GetNumPoints();

            // The quadrature points at the points at which we
            // want to know the metrics (as specified in tbasis)
            LibUtilities::PointsKey pkey_tbasis(tbasis[0]->GetPointsKey());
            int nquad_tbasis = pkey_tbasis.GetNumPoints();

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            mPointsKey[0] = pkey_tbasis;

            // setup temp storage
            Array<OneD, Array<OneD,NekDouble> > der_map   (coordim);

            mDeriv = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(1);
            mDeriv[0] = Array<OneD, Array<OneD, NekDouble> >(mCoordDim);

            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                der_map[i]    = Array<OneD,NekDouble>(nquad_map);
                mDeriv[0][i]  = Array<OneD,NekDouble>(nquad_tbasis);

                // Transform from coefficient space to physical space
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),
                                    Coords[i]->UpdatePhys());
                // Take the derivative (calculated at the points as specified
                // in 'Coords')
                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(), der_map[i]);

                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics
                //   ('tbasis')
                if( pkey_map == pkey_tbasis )
                {
                    mDeriv[0][i] = der_map[i];
                }
                else
                {
                    LibUtilities::Interp1D(pkey_map, der_map[i], pkey_tbasis,
                                           mDeriv[0][i]);
                }
            }

            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation
            // metrics
            SetUpJacGmat1D();

            // 2. the jacobian muliplied with the quadrature weights
            if(QuadMetrics)
            {
                SetUpQuadratureMetrics(shape,tbasis);
            }
            // 3. A combination of the metrics above that allows
            //    for more efficient evaluation of the laplacian
            if(LaplMetrics)
            {
                SetUpLaplacianMetrics (shape,tbasis);
            }
        }


        /**
         * Create a copy of an existing GeomFactors1D object.
         */
        GeomFactors1D::GeomFactors1D(const GeomFactors1D& S) :
            GeomFactors(S)
        {
        }


        /**
         *
         */
        GeomFactors1D::~GeomFactors1D()
        {
        }


        /**
         * Constructs the one-dimensional Jacobian matrix
         * \f[J_{1D} = J = \sqrt{\frac{d x_1}{d _\xi}.
         *   \frac{d x_1}{d _\xi} + \frac{d x_2}{d _\xi}.\frac{d x_2}{d _\xi}
         *   + \frac{d x_3}{d _\xi}\frac{d x_3}{d _\xi}} \f]
         * and local factors (\f$ g \f$ matrix)
         * \f[g_i = \frac{d_\xi}{x_i}\f]
         *
         * @param   der         Local derivative matrix of size Coordim x Nquad.
         *                      Derivatives at each quadrature point.
         */
        void GeomFactors1D::SetUpJacGmat1D()
        {
            ASSERTL1(mDeriv[0].num_elements()==mCoordDim,
                     "The dimension of array mDeriv does not match the "
                     "coordinate dimension");
            int i;
            int nquad = mPointsKey[0].GetNumPoints();
            ASSERTL1(mDeriv[0][0].num_elements() == nquad,
                     "Number of quadrature points do not match");

            // If regular or moving geometry.
            if(( mType == eRegular)||
               ( mType == eMovingRegular))
            {
                // Allocate storage.
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(mCoordDim,1,0.0);

                // Loop over each dimension in the coordinate system.
                for(i = 0; i < mCoordDim; ++i)
                {
                    // i-th component of normal = 1/derivative_i.
                    m_gmat[i][0] = (fabs(mDeriv[0][i][0])
                            > NekConstants::kNekZeroTol) ? 1.0/mDeriv[0][i][0]: 0.0;
                    // i-th component contribution to Jacobian.
                    m_jac[0]    += mDeriv[0][i][0]*mDeriv[0][i][0];
                }
                // take square root
                m_jac[0] = sqrt(m_jac[0]);
            }
            // If deformed geometry
            else
            {
                m_jac     = Array<OneD, NekDouble>(nquad,0.0);
                m_gmat    = Array<TwoD, NekDouble>(mCoordDim,nquad);

                // invert local derivative for gmat;
                for(i = 0; i < mCoordDim; ++i)
                {
                    for(int j = 0; j < nquad; ++j)
                    {
                        m_gmat[i][j] = (fabs(mDeriv[0][i][j])
                            > NekConstants::kNekZeroTol) ? 1.0/mDeriv[0][i][j] : 0.0;
                    }
                    // compute jacobian for this dimension.
                    Vmath::Vvtvp(nquad,mDeriv[0][i],1,mDeriv[0][i],1,m_jac,1,m_jac,1);
                }
                Vmath::Vsqrt(nquad,m_jac,1,m_jac,1);
            }
        }


        /**
         * It is assumed that this 1D segment forms the single edge \a edge of
         * a 2D geometry described by \a geom. It retrieves the edge normals
         * from the associated 2D expansion and orientates it correctly.
         * @param   geom        Geometry (2D) to which this segment forms an
         *                      edge.
         * @param   edge        Edge number of this segment on \a geom.
         * @param   to_key      PointsKey describing the quadrature points at
         *                      which to evaluate the normals (typically the
         *                      PointsKey of the associated segment).
         */
        void GeomFactors1D::v_ComputeNormals(
                            const GeometrySharedPtr &geom,
                            const int edge,
                            const LibUtilities::PointsKey &to_key)
        {
            int k;
            int nq      = to_key.GetNumPoints();

            // Ensure we have a 2D geometry.
            Geometry2DSharedPtr g;
            if (!(g = boost::dynamic_pointer_cast<Geometry2D>(geom)))
            {
                ASSERTL0(false, "FAIL");
            }
            GeomFactorsSharedPtr gf = geom->GetMetricInfo();

            // Retrieve the GeomFactors object describing the shape geometry
            // and generate the normals to the edge.
            mNormal = Array<OneD, Array<OneD, NekDouble> >(mCoordDim);
            for (k = 0; k < mCoordDim; ++k)
            {
                mNormal[k] = Array<OneD, NekDouble>(nq);
            }
            gf->ComputeEdgeNormals(edge, to_key, mNormal);

            if(g->GetEorient(edge) == StdRegions::eBackwards)
            {
                for(k = 0; k < mCoordDim; ++k)
                {
                    if(gf->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Reverse(nq, mNormal[k], 1, mNormal[k],1);
                    }
                    Vmath::Neg(nq,mNormal[k],1);
                }
            }
        }

    }
}


