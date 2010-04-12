#include <SpatialDomains/GeomFactors3D.h>

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
         * @param   SetUpQuadratureMetrics  ?
         * @param   SetUpLaplacianMetrics   ?
         */
        GeomFactors3D::GeomFactors3D(const GeomType gtype,
                          const int coordim,
                          const Array<OneD, const StdRegions
                                            ::StdExpansion3DSharedPtr> &Coords,
                          const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis,
                          const bool QuadMetrics,
                          const bool LaplMetrics) :
            GeomFactors(gtype,3,coordim,QuadMetrics,LaplMetrics)
        {
            ASSERTL1((coordim == 3),
                     "The coordinate dimension should be to three"
                     "for three-dimensional elements");
            ASSERTL1(tbasis.num_elements()==3,"tbasis should be an array of size three");
            ASSERTL1(LaplMetrics?QuadMetrics:true,
                     "SetUpQuadratureMetrics should be true if SetUpLaplacianMetrics is true");

            StdRegions::ExpansionType shape = Coords[0]->DetExpansionType();

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
            mPointsKey[0] = pkey0_tbasis;
            mPointsKey[1] = pkey1_tbasis;
            mPointsKey[2] = pkey2_tbasis;

            // setup temp storage
            Array<OneD, Array<OneD,NekDouble> > d1_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d2_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d3_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > d1_tbasis(coordim);
            Array<OneD, Array<OneD,NekDouble> > d2_tbasis(coordim);
            Array<OneD, Array<OneD,NekDouble> > d3_tbasis(coordim);

            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                d1_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d2_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d3_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d1_tbasis[i] = Array<OneD,NekDouble>(nqtot_tbasis);
                d2_tbasis[i] = Array<OneD,NekDouble>(nqtot_tbasis);
                d3_tbasis[i] = Array<OneD,NekDouble>(nqtot_tbasis);

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
                    d1_tbasis[i] = d1_map[i];
                    d2_tbasis[i] = d2_map[i];
                    d3_tbasis[i] = d3_map[i];
                }
                else
                {
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,             pkey2_map,    d1_map[i],
                                           pkey0_tbasis, pkey1_tbasis,          pkey2_tbasis, d1_tbasis[i]);
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,             pkey2_map,    d2_map[i],
                                           pkey0_tbasis, pkey1_tbasis,          pkey2_tbasis, d2_tbasis[i]);
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,             pkey2_map,    d3_map[i],
                                           pkey0_tbasis, pkey1_tbasis,          pkey2_tbasis, d3_tbasis[i]);
                }
            }

            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation      metrics
            SetUpJacGmat3D(d1_tbasis,d2_tbasis,d3_tbasis);

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
         *
         */
        void GeomFactors3D::SetUpJacGmat3D(
                        const Array<OneD, Array<OneD, NekDouble> > d1,
                        const Array<OneD, Array<OneD, NekDouble> > d2,
                        const Array<OneD, Array<OneD, NekDouble> > d3)
        {
            ASSERTL1(d1.num_elements()==mCoordDim,"The dimension of array d1 does not"
                     "match the coordinate dimension");
            ASSERTL1(d2.num_elements()==mCoordDim,"The dimension of array d2 does not"
                     "match the coordinate dimension");
            ASSERTL1(d3.num_elements()==mCoordDim,"The dimension of array d3 does not"
                     "match the coordinate dimension");

            int nqtot = mPointsKey[0].GetNumPoints() *
                        mPointsKey[1].GetNumPoints() *
                        mPointsKey[2].GetNumPoints();

            ASSERTL1(d1[0].num_elements() == nqtot,"Number of quadrature points do not match");
            ASSERTL1(d2[0].num_elements() == nqtot,"Number of quadrature points do not match");
            ASSERTL1(d3[0].num_elements() == nqtot,"Number of quadrature points do not match");

            if((mType == eRegular)||(mType == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(3*mCoordDim,1,0.0);

                m_jac[0] =  d1[0][0]*(d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])
                           -d2[0][0]*(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0])
                           +d3[0][0]*(d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0]);

                ASSERTL1(m_jac[0] > 0, "3D Regular Jacobian is not positive");

                // Spen's book page 160
                m_gmat[0][0] =  (d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])/m_jac[0];  // d xi_1/d x_1
                m_gmat[1][0] = -(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0])/m_jac[0];  // d xi_2/d x_1
                m_gmat[2][0] =  (d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0])/m_jac[0];  // d xi_3/d x_1
                m_gmat[3][0] = -(d2[0][0]*d3[2][0] - d3[0][0]*d2[2][0])/m_jac[0];  // d xi_1/d x_2
                m_gmat[4][0] =  (d1[0][0]*d3[2][0] - d3[0][0]*d1[2][0])/m_jac[0];  // d xi_2/d x_2
                m_gmat[5][0] = -(d1[0][0]*d2[2][0] - d2[0][0]*d1[2][0])/m_jac[0];  // d xi_3/d x_2
                m_gmat[6][0] =  (d2[0][0]*d3[1][0] - d3[0][0]*d2[1][0])/m_jac[0];  // d xi_1/d x_3
                m_gmat[7][0] = -(d1[0][0]*d3[1][0] - d3[0][0]*d1[1][0])/m_jac[0];  // d xi_2/d x_3
                m_gmat[8][0] =  (d1[0][0]*d2[1][0] - d2[0][0]*d1[1][0])/m_jac[0];  // d xi_3/d x_3
            }
            else // Deformed case
            {
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(3*mCoordDim,nqtot,0.0);

                // Derivatives of the form: dj[i] = dx_(i+1)/dxi_j

                // Spencers book page 160
                // g[0] = d xi_1/d x_1
                Vmath::Vmul (nqtot,&d3[1][0],1,&d2[2][0],1,&m_gmat[0][0],1);
                Vmath::Vvtvm(nqtot,&d2[1][0],1,&d3[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
                // g[1] = d xi_2/d x_1
                Vmath::Vmul (nqtot,&d1[1][0],1,&d3[2][0],1,&m_gmat[1][0],1);
                Vmath::Vvtvm(nqtot,&d3[1][0],1,&d1[2][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
                // g[2] = d xi_3/d x_1
                Vmath::Vmul (nqtot,&d2[1][0],1,&d1[2][0],1,&m_gmat[2][0],1);
                Vmath::Vvtvm(nqtot,&d1[1][0],1,&d2[2][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
                // g[3] = d xi_2/d x_1
                Vmath::Vmul (nqtot,&d2[0][0],1,&d3[2][0],1,&m_gmat[3][0],1);
                Vmath::Vvtvm(nqtot,&d3[0][0],1,&d2[2][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
                // g[4] = d xi_2/d x_2
                Vmath::Vmul (nqtot,&d3[0][0],1,&d1[2][0],1,&m_gmat[4][0],1);
                Vmath::Vvtvm(nqtot,&d1[0][0],1,&d3[2][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
                // g[5] = d xi_2/d x_3
                Vmath::Vmul (nqtot,&d1[0][0],1,&d2[2][0],1,&m_gmat[5][0],1);
                Vmath::Vvtvm(nqtot,&d2[0][0],1,&d1[2][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
                // g[6] = d xi_3/d x_1
                Vmath::Vmul (nqtot,&d3[0][0],1,&d2[1][0],1,&m_gmat[6][0],1);
                Vmath::Vvtvm(nqtot,&d2[0][0],1,&d3[1][0],1,&m_gmat[6][0],1,&m_gmat[6][0],1);
                // g[7] = d xi_3/d x_2
                Vmath::Vmul (nqtot,&d1[0][0],1,&d3[1][0],1,&m_gmat[7][0],1);
                Vmath::Vvtvm(nqtot,&d3[0][0],1,&d1[1][0],1,&m_gmat[7][0],1,&m_gmat[7][0],1);
                // g[8] = d xi_3/d x_3
                Vmath::Vmul (nqtot,&d2[0][0],1,&d1[1][0],1,&m_gmat[8][0],1);
                Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&m_gmat[8][0],1,&m_gmat[8][0],1);

                // J3D - Spencers book page 158
                Vmath::Vmul (nqtot,&d1[0][0],1,&m_gmat[0][0],1,&m_jac[0],1);
                Vmath::Vvtvp(nqtot,&d2[0][0],1,&m_gmat[1][0],1,&m_jac[0],1,&m_jac[0],1);
                Vmath::Vvtvp(nqtot,&d3[0][0],1,&m_gmat[2][0],1,&m_jac[0],1,&m_jac[0],1);

                ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0, "3D Deformed Jacobian is not positive");

                // Scale g[i] by 1/J3D
                Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[6][0],1,&m_jac[0],1,&m_gmat[6][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[7][0],1,&m_jac[0],1,&m_gmat[7][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[8][0],1,&m_jac[0],1,&m_gmat[8][0],1);
            }
        }


        /**
         *
         */
        void GeomFactors3D::v_SetUpQuadratureMetrics(StdRegions::ExpansionType shape,
                                                   const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == mExpDim,"Inappropriate dimension of tbasis");

            int i,j,k;
            int nquad0 = mPointsKey[0].GetNumPoints();
            int nquad1 = mPointsKey[1].GetNumPoints();
            int nquad2 = mPointsKey[2].GetNumPoints();
            int nqtot  = nquad0*nquad1*nquad2;

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
            const Array<OneD, const NekDouble>& w2 = tbasis[2]->GetW();

            // Multiply the jacobian with the quadrature weights
            switch(shape)
            {
            case StdRegions::eHexahedron:
                {
                    for(k = 0; k < nquad2; ++k)
                    {
                        for(j = 0; j < nquad1; ++j)
                        {
                            Vmath::Vmul(nquad0,m_weightedjac.get()+(k*nquad1+j)*nquad0,1,
                                    w0.get(),1,m_weightedjac.get()+(k*nquad1+j)*nquad0,1);
                        }
                    }

                    for(k = 0; k < nquad2; ++k)
                    {
                        for(i = 0; i < nquad0; ++i)
                        {
                            Vmath::Vmul(nquad1,m_weightedjac.get()+k*nquad0*nquad1+i,nquad0,
                                    w1.get(),1,m_weightedjac.get()+k*nquad0*nquad1+i,nquad0);
                        }
                    }

                    for(j = 0; j < nquad1; ++j)
                    {
                        for(i = 0; i < nquad0; ++i)
                        {
                            Vmath::Vmul(nquad2,m_weightedjac.get()+j*nquad0+i,nquad0*nquad1,
                                    w1.get(),1,m_weightedjac.get()+j*nquad0+i,nquad0*nquad1);
                        }
                    }
                }
                break;
            case StdRegions::eTetrahedron:
            case StdRegions::ePrism:
            case StdRegions::ePyramid:
                {
                    ASSERTL0(false, "SetUpQuadratureWeights: Need to implement quadrature weights for this shape.");

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

    }
}
