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

            // The jacobian seems to be calculated wrongly:
            // Rather than the formula:
            //    m_jac[0] =  d1[0][0]*(d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0])
            //               -d1[1][0]*(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0])
            //               +d1[2][0]*(d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0]);
            // I think this should be (According to Spencer's book page 158):
            //    m_jac[0] =  d1[0][0]*(d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])
            //               -d2[0][0]*(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0])
            //               +d3[0][0]*(d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0]);
            // Please verify and update this, also for the deformed case...
            //
            // In addition, m_gmat[2][0] seems to be calculated differently
            // as in Spencer's book page 160)
            // Currently, it is:
            // m_gmat[2][0] =  (d1[1][0]*d2[2][0] - d1[2][0]*d2[1][0])/m_jac[0];
            // but Spencer's book would suggest:
            // m_gmat[2][0] =  (d1[1][0]*d2[2][0] - d1[2][0]*d3[1][0])/m_jac[0];
            // I am not sure which version is right. please verify!
            // Also check the deformed case.
            //
            // Update:
            // Checked both expressions on Spencer's book:
            // - J3D from pg 158 is fine, so the implementation.
            // - There is a typo on d xi_3/dx_1. The third term is *not*
            //   dx_2/dxi_3, but dx_2/dxi_2.
            // - I guess terms or commentaries are swaped below; where you read
            //   d xi_M/d x_N should be d xi_N/d x_M. In other words,
            //   transposed.
            // - Deformed case not checked.
            //
            // Update 2 (pvos):
            // I did change the formulation of the jacobian from the first to the
            // second version. I think this should be correct
            // (certainly if you know that dj[i] = dx_i/dxi_j)
            // I only updated the regular geometry. Deformed still to be done

            if((mType == eRegular)||(mType == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(3*mCoordDim,1,0.0);

                // J3D: Determinant of three-dimensional Jacobian
//                 m_jac[0] = d1[0][0]*( d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0] )
//                           -d1[1][0]*( d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0] )
//                           +d1[2][0]*( d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0] );
                m_jac[0] =  d1[0][0]*(d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])
                           -d2[0][0]*(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0])
                           +d3[0][0]*(d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0]);

                ASSERTL1(m_jac[0] > 0, "3D Regular Jacobian is not positive");
                // Spen's book page 160
                m_gmat[0][0] =  (d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0])/m_jac[0];  // d xi_1/d x_1
                m_gmat[1][0] = -(d1[1][0]*d3[2][0] - d1[2][0]*d3[1][0])/m_jac[0];  // d xi_2/d x_1
                m_gmat[2][0] =  (d1[1][0]*d2[2][0] - d1[2][0]*d2[1][0])/m_jac[0];  // d xi_3/d x_1
                m_gmat[3][0] = -(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0])/m_jac[0];  // d xi_1/d x_2
                m_gmat[4][0] =  (d1[0][0]*d3[2][0] - d1[2][0]*d3[0][0])/m_jac[0];  // d xi_2/d x_2
                m_gmat[5][0] = -(d1[0][0]*d2[2][0] - d1[2][0]*d2[0][0])/m_jac[0];  // d xi_3/d x_2
                m_gmat[6][0] =  (d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0])/m_jac[0];  // d xi_1/d x_3
                m_gmat[7][0] = -(d1[0][0]*d3[1][0] - d1[1][0]*d3[0][0])/m_jac[0];  // d xi_2/d x_3
                m_gmat[8][0] =  (d1[0][0]*d2[1][0] - d1[1][0]*d2[0][0])/m_jac[0];  // d xi_3/d x_3
            }
            else // Deformed case
            {
                ASSERTL0(false,"This routine needs corrections. Please see notes in the code...");
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(3*mCoordDim,nqtot,0.0);

                // set up Jacobian
                Array<OneD,NekDouble> tmp[3] = {Array<OneD, NekDouble>(nqtot),
                                                Array<OneD, NekDouble>(nqtot),
                                                Array<OneD, NekDouble>(nqtot)};
                // g[0]
                Vmath::Vmul (nqtot,&d2[2][0],1,&d3[1][0],1,&tmp[0][0],1);
                Vmath::Vvtvm(nqtot,&d2[1][0],1,&d3[2][0],1,&tmp[0][0],1,&tmp[0][0],1);
                //g[1]
                Vmath::Vmul (nqtot,&d2[0][0],1,&d3[2][0],1,&tmp[1][0],1);
                Vmath::Vvtvm(nqtot,&d2[2][0],1,&d3[0][0],1,&tmp[1][0],1,&tmp[1][0],1);
                //g[2]
                Vmath::Vmul (nqtot,&d2[1][0],1,&d3[0][0],1,&tmp[2][0],1);
                Vmath::Vvtvm(nqtot,&d2[0][0],1,&d3[1][0],1,&tmp[2][0],1,&tmp[2][0],1);

                // J3D
                Vmath::Vmul (nqtot,&d1[0][0],1,&tmp[0][0],1,&m_jac[0],1);
                Vmath::Vvtvp(nqtot,&d1[1][0],1,&tmp[1][0],1,&m_jac[0],1,&m_jac[0],1);
                Vmath::Vvtvp(nqtot,&d1[2][0],1,&tmp[2][0],1,&m_jac[0],1,&m_jac[0],1);

                ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0, "3D Deformed Jacobian is not positive");

                // d xi_1/d x_1
                Vmath::Vmul (nqtot,&d2[2][0],1,&d3[1][0],1,&m_gmat[0][0],1);
                Vmath::Vvtvm(nqtot,&d2[1][0],1,&d3[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);

                // d xi_1/d x_2
                Vmath::Vmul (nqtot,&d1[1][0],1,&d3[2][0],1,&m_gmat[1][0],1);
                Vmath::Vvtvm(nqtot,&d1[2][0],1,&d3[1][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);

                // d xi_1/d x_3
                Vmath::Vmul (nqtot,&d1[2][0],1,&d2[1][0],1,&m_gmat[2][0],1);
                Vmath::Vvtvm(nqtot,&d1[1][0],1,&d2[2][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);

                // d xi_2/d x_1
                Vmath::Vmul (nqtot,&d2[0][0],1,&d3[2][0],1,&m_gmat[3][0],1);
                Vmath::Vvtvm(nqtot,&d2[2][0],1,&d3[0][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);

                // d xi_2/d x_2
                Vmath::Vmul (nqtot,&d1[2][0],1,&d3[0][0],1,&m_gmat[4][0],1);
                Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[2][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);

                // d xi_2/d x_3
                Vmath::Vmul (nqtot,&d1[0][0],1,&d2[2][0],1,&m_gmat[5][0],1);
                Vmath::Vvtvm(nqtot,&d1[2][0],1,&d2[0][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);

                // d xi_3/d x_1
                Vmath::Vmul (nqtot,&d2[1][0],1,&d3[0][0],1,&m_gmat[6][0],1);
                Vmath::Vvtvm(nqtot,&d2[0][0],1,&d3[1][0],1,&m_gmat[6][0],1,&m_gmat[6][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[6][0],1,&m_jac[0],1,&m_gmat[6][0],1);

                // d xi_3/d x_2
                Vmath::Vmul (nqtot,&d1[0][0],1,&d3[1][0],1,&m_gmat[7][0],1);
                Vmath::Vvtvm(nqtot,&d1[1][0],1,&d3[0][0],1,&m_gmat[7][0],1,&m_gmat[7][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[7][0],1,&m_jac[0],1,&m_gmat[7][0],1);

                // d xi_3/d x_3
                Vmath::Vmul (nqtot,&d1[1][0],1,&d2[0][0],1,&m_gmat[8][0],1);
                Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&m_gmat[8][0],1,&m_gmat[8][0],1);
                Vmath::Vdiv(nqtot,&m_gmat[8][0],1,&m_jac[0],1,&m_gmat[8][0],1);
            }
        }

    }
}
