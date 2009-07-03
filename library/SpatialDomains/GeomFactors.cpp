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
//  License for the specific language governing rights and limitations under
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
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry2D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        GeomFactors::GeomFactors(void):
            m_gtype(eRegular),
            m_expdim(0), 
            m_coordim(0),
            m_quadratureMetricsFlag(false),
            m_laplacianMetricsFlag(false)
        {
        }

        GeomFactors::GeomFactors(const GeomType gtype,
                                 const int expdim, 
                                 const int coordim):
            m_gtype(gtype),
            m_expdim(expdim), 
            m_coordim(coordim),
            m_quadratureMetricsFlag(false),
            m_laplacianMetricsFlag(false),
            m_pointsKey(expdim)
        {
        }

        /** \brief One dimensional geometric factors based on one, two or three
        dimensional coordinate description 

        - Calculate the local factors \f$ d \xi/d x_i \f$ 
        as a recipricol derivative \f$ d x_i/d\xi \f$

        - Calcuate the jacobian as \f$ J = \sqrt{\frac{d x_1}{d _\xi}.
        \frac{d x_1}{d _\xi} + \frac{d x_2}{d _\xi}.\frac{d x_2}{d _\xi}  + 
        \frac{d x_3}{d _\xi}\frac{d x_3}{d _\xi}} \f$

        **/

        GeomFactors::GeomFactors(const GeomType gtype, 
                                 const int coordim, 
                                 const Array<OneD, const StdRegions::StdExpansion1DSharedPtr>  &Coords,
                                 const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis,
                                 const bool SetUpQuadratureMetrics,
                                 const bool SetUpLaplacianMetrics):
            m_gtype(gtype), 
            m_coordim(coordim), 
            m_expdim(1),
            m_quadratureMetricsFlag(false),
            m_laplacianMetricsFlag(false),
            m_pointsKey(1)
        {

            ASSERTL1((coordim == 1)||(coordim == 2)||(coordim == 3),
                     "The coordinate dimension should be equal to one, two or three"
                     "for one-dimensional elements");
            ASSERTL1(tbasis.num_elements()==1,"tbasis should be an array of size one");
            ASSERTL1(SetUpLaplacianMetrics?SetUpQuadratureMetrics:true,
                     "SetUpQuadratureMetrics should be true if SetUpLaplacianMetrics is true");

            StdRegions::ExpansionType shape = Coords[0]->DetExpansionType();

            // The quadrature points of the mapping 
            // (as specified in Coords)
            LibUtilities::PointsKey pkey_map(Coords[0]->GetBasis(0)->GetPointsKey());
            int nquad_map = pkey_map.GetNumPoints();

            // The quadrature points at the points at which we
            // want to know the metrics (as specified in tbasis)
            LibUtilities::PointsKey pkey_tbasis(tbasis[0]->GetPointsKey());
            int nquad_tbasis = pkey_tbasis.GetNumPoints();

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            m_pointsKey[0] = pkey_tbasis;

            // setup temp storage
            Array<OneD, Array<OneD,NekDouble> > der_map   (coordim);
            Array<OneD, Array<OneD,NekDouble> > der_tbasis(coordim);
            
            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                der_map[i]    = Array<OneD,NekDouble>(nquad_map);
                der_tbasis[i] = Array<OneD,NekDouble>(nquad_tbasis);

                // Transform from coefficient space to physical space
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),Coords[i]->UpdatePhys());
                // Take the derivative (calculated at the points as specified in 'Coords')
                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),der_map[i]);

                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics ('tbasis')
                if( pkey_map == pkey_tbasis )
                {
                    der_tbasis[i] = der_map[i];
                }
                else
                {
                    LibUtilities::Interp1D(pkey_map, der_map[i], pkey_tbasis, der_tbasis[i]);
                }
            }

            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation metrics
            SetUpJacGmat1D(der_tbasis);

            /* The routines below need implementation

            // 2. the jacobian muliplied with the quadrature weights
            if(SetUpQuadratureMetrics)
            {
            SetUpQuadratureMetrics1D(shape,tbasis);
            }
            // 3. A combination of the metrics above that allows
            //    for more efficient evaluation of the laplacian
            if(SetUpLaplacianMetrics)
            {
            SetUpLaplacianMetrics1D (shape,tbasis);
            }

            */
        }
        
        /**
        \brief Two dimensional geometric factors based on two or three
        dimensional coordinate description

        The geometric factors are evaluated by considering the
        mapping to a coordinate system based on the local tangent
        vectors (which are the local derivatives of the global
        coordinates) and the normal \f$ \bf g \f$ to these two
        tangent vectors. We therefore use the 3 x 3 relationships but
        assume that \f$ \partial x_1/\partial \xi_3 = { g_1},\,
        \partial x_2/\partial \xi_3 = { g_2},\, \partial x_3/\partial
        \xi_3 = { g_3} \f$ i.e.

        \f$ {\bf g }= \left [ \begin{array}{c} g_1 \\ g_2 \\ g_3 \end{array} 
        \right ] = \frac{\partial {\bf x}}{\partial \xi_1} \times 
        \frac{\partial {\bf x}}{\partial \xi_2} =   \left [ \begin{array}{c}
        \frac{\partial x_2}{\partial \xi_1}   \frac{\partial x_3}{\partial \xi_2} - 
        \frac{\partial x_3}{\partial \xi_1}   \frac{\partial x_2}{\partial \xi_2} \\
        \frac{\partial x_3}{\partial \xi_1}   \frac{\partial x_1}{\partial \xi_2} - 
        \frac{\partial x_1}{\partial \xi_1}   \frac{\partial x_3}{\partial \xi_2} \\
        \frac{\partial x_1}{\partial \xi_1}   \frac{\partial x_2}{\partial \xi_2} - 
        \frac{\partial x_2}{\partial \xi_1}   \frac{\partial x_1}{\partial \xi_2}
        \end{array} \right ] \f$
    
        The geometric factors are then given by:

        \f$ \begin{array}{cc}
        \frac{\partial \xi_1}{\partial x_1} = \frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_2} {g_3} - \frac{\partial x_3}{\partial \xi_2} {g_2} \right ) & 
        \frac{\partial \xi_1}{\partial x_2} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} {g_3} - \frac{\partial x_3}{\partial \xi_2} {g_1} \right )\\
        \frac{\partial \xi_1}{\partial x_3} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} {g_2} - \frac{\partial x_2}{\partial \xi_2} {g_1} \right )&
        \frac{\partial \xi_2}{\partial x_1} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_1} {g_3} - \frac{\partial x_3}{\partial \xi_1} {g_2} \right ) \\
        \frac{\partial \xi_2}{\partial x_2} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} {g_3} - \frac{\partial x_3}{\partial \xi_1} {g_1} \right ) &
        \frac{\partial \xi_2}{\partial x_3} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} {g_2} - \frac{\partial x_2}{\partial \xi_1} {g_1} \right )
        \end{array} \f$

        where

        \f$ J_{3D} = 
        {g_3} \left ( 
        \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_2} - 
        \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_2}{\partial \xi_1}\right )
        + {g_2}\left (  
        \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_1}  -
        \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2}\right )
        + {g_1} \left ( 
        \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2} -
        \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_1} 
        \right ) \f$

        and the two-dimensional surface Jacobian  is given by
        \f$ J = \left | \frac{\partial {\bf x}}{\partial \xi_1} \times 
        \frac{\partial {\bf x}}{\partial \xi_2} \right | = \sqrt{J_{3D}} \f$

        **/
        GeomFactors::GeomFactors(const GeomType gtype, 
                                 const int coordim,
                                 const Array<OneD, const StdRegions::StdExpansion2DSharedPtr> &Coords,
                                 const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis,
                                 const bool SetUpQuadratureMetrics,
                                 const bool SetUpLaplacianMetrics):
            m_gtype(gtype), 
            m_coordim(coordim), 
            m_expdim(2),
            m_quadratureMetricsFlag(false),
            m_laplacianMetricsFlag(false),
            m_pointsKey(2)
        {
            ASSERTL1((coordim == 2)||(coordim == 3),
                     "The coordinate dimension should be equal to two or three"
                     "for two-dimensional elements");
            ASSERTL1(tbasis.num_elements()==2,"tbasis should be an array of size two");
            ASSERTL1(SetUpLaplacianMetrics?SetUpQuadratureMetrics:true,
                     "SetUpQuadratureMetrics should be true if SetUpLaplacianMetrics is true");

            StdRegions::ExpansionType shape = Coords[0]->DetExpansionType();

            // The quadrature points of the mapping 
            // (as specified in Coords)
            LibUtilities::PointsKey pkey0_map(Coords[0]->GetBasis(0)->GetPointsKey());
            LibUtilities::PointsKey pkey1_map(Coords[0]->GetBasis(1)->GetPointsKey());
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
            Array<OneD, Array<OneD,NekDouble> > d1_tbasis(coordim);
            Array<OneD, Array<OneD,NekDouble> > d2_tbasis(coordim);
            
            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
                d1_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d2_map[i]    = Array<OneD,NekDouble>(nqtot_map);
                d1_tbasis[i] = Array<OneD,NekDouble>(nqtot_tbasis);
                d2_tbasis[i] = Array<OneD,NekDouble>(nqtot_tbasis);

                // Transform from coefficient space to physical space
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),Coords[i]->UpdatePhys());
                // Take the derivative (calculated at the points as specified in 'Coords')
                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1_map[i],d2_map[i]);

                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics ('tbasis')
                if( (pkey0_map == pkey0_tbasis) &&
                    (pkey1_map == pkey1_tbasis) )
                {
                    d1_tbasis[i] = d1_map[i];
                    d2_tbasis[i] = d2_map[i];
                }
                else
                {
                    LibUtilities::Interp2D(pkey0_map, pkey1_map, d1_map[i], pkey0_tbasis, pkey1_tbasis, d1_tbasis[i]);
                    LibUtilities::Interp2D(pkey0_map, pkey1_map, d2_map[i], pkey0_tbasis, pkey1_tbasis, d2_tbasis[i]);
                }
            }

            // Setting up two tangential basis vectors
            SetUpTangentialbasis(d1_tbasis, d2_tbasis);
              
            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation metrics
            SetUpJacGmat2D(d1_tbasis,d2_tbasis);
            // 2. the jacobian muliplied with the quadrature weights
            if(SetUpQuadratureMetrics)
            {
                SetUpQuadratureMetrics2D(shape,tbasis);
            }
            // 3. A combination of the metrics above that allows
            //    for more efficient evaluation of the laplacian
            if(SetUpLaplacianMetrics)
            {
                SetUpLaplacianMetrics2D (shape,tbasis);
            }
        }
        

        GeomFactors::GeomFactors(const GeomType gtype, 
                                 const int coordim,
                                 const Array<OneD, const StdRegions::StdExpansion3DSharedPtr> &Coords,
                                 const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis,
                                 const bool SetUpQuadratureMetrics,
                                 const bool SetUpLaplacianMetrics):
            m_gtype(gtype), 
            m_coordim(coordim), 
            m_expdim(3),
            m_quadratureMetricsFlag(false),
            m_laplacianMetricsFlag(false),
            m_pointsKey(3)
        {
            ASSERTL1((coordim == 3),
                     "The coordinate dimension should be to three"
                     "for three-dimensional elements");
            ASSERTL1(tbasis.num_elements()==3,"tbasis should be an array of size three");
            ASSERTL1(SetUpLaplacianMetrics?SetUpQuadratureMetrics:true,
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
            m_pointsKey[0] = pkey0_tbasis;
            m_pointsKey[1] = pkey1_tbasis;
            m_pointsKey[2] = pkey2_tbasis;

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
                // Take the derivative (calculated at the points as specified in 'Coords')
                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1_map[i],d2_map[i],d3_map[i]);

                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics ('tbasis')
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
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,    pkey2_map,    d1_map[i], 
                                           pkey0_tbasis, pkey1_tbasis, pkey2_tbasis, d1_tbasis[i]);
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,    pkey2_map,    d2_map[i], 
                                           pkey0_tbasis, pkey1_tbasis, pkey2_tbasis, d2_tbasis[i]);
                    LibUtilities::Interp3D(pkey0_map,    pkey1_map,    pkey2_map,    d3_map[i], 
                                           pkey0_tbasis, pkey1_tbasis, pkey2_tbasis, d3_tbasis[i]);
                }
            }

            // Based upon these derivatives, calculate:
            // 1. The (determinant of the) jacobian and the differentation metrics
            SetUpJacGmat3D(d1_tbasis,d2_tbasis,d3_tbasis);

            /* The routines below need implementation

            // 2. the jacobian muliplied with the quadrature weights
            if(SetUpQuadratureMetrics)
            {
            SetUpQuadratureMetrics3D(shape,tbasis);
            }
            // 3. A combination of the metrics above that allows
            //    for more efficient evaluation of the laplacian
            if(SetUpLaplacianMetrics)
            {
            SetUpLaplacianMetrics3D (shape,tbasis);
            }

            */
        }  
   
        void GeomFactors::SetUpJacGmat1D(const Array<OneD, Array<OneD, NekDouble> > der)
        {
            ASSERTL1(der.num_elements()==m_coordim,"The dimension of array der does not"
                     "match the coordinate dimension");
            int i;
            int nquad = m_pointsKey[0].GetNumPoints();
            ASSERTL1(der[0].num_elements() == nquad,"Number of quadrature points do not match");

            NekDouble fac0,fac1;

            if(( m_gtype == eRegular)||
               ( m_gtype == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(m_coordim,1,0.0);
                
                fac0 = fac1 = 0.0;
                for(i = 0; i < m_coordim; ++i)
                {
                    m_gmat[i][0] = (fabs(der[i][0]) > NekConstants::kNekZeroTol)? 1.0/der[i][0]: 0.0;
                    m_jac[0]    += der[i][0]*der[i][0];

                    fac0 += m_gmat[i][0]*m_gmat[i][0];
                }
                m_jac[0] = sqrt(m_jac[0]);
                
                fac0 = fac1 = sqrt(fac0);
            }
            else
            {
                m_jac     = Array<OneD, NekDouble>(nquad,0.0);
                m_gmat    = Array<TwoD, NekDouble>(m_coordim,nquad);

                // invert local derivative for gmat;
                fac0 = fac1 = 0.0;
                for(i = 0; i < m_coordim; ++i)
                {
                    for(int j = 0; j < nquad; ++j)
                    {
                        m_gmat[i][j] = (fabs(der[i][j]) > NekConstants::kNekZeroTol)? 1.0/der[i][j]: 0.0;
                    }
                    Vmath::Vvtvp(nquad,der[i],1,der[i],1,m_jac,1,m_jac,1);

                    fac0 += m_gmat[i][0]*m_gmat[i][0];
                    fac1 += m_gmat[i][nquad-1]*m_gmat[i][nquad-1];
                }
                Vmath::Vsqrt(nquad,m_jac,1,m_jac,1);

                fac0 = sqrt(fac0);
                fac1 = sqrt(fac1);
            } 
        }
     
        void GeomFactors::SetUpJacGmat2D(const Array<OneD, Array<OneD, NekDouble> > d1,
                                         const Array<OneD, Array<OneD, NekDouble> > d2)
        {
            ASSERTL1(d1.num_elements()==m_coordim,"The dimension of array d1 does not"
                     "match the coordinate dimension");
            ASSERTL1(d2.num_elements()==m_coordim,"The dimension of array d2 does not"
                     "match the coordinate dimension");

            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints();

            ASSERTL1(d1[0].num_elements() == nqtot,"Number of quadrature points do not match");
            ASSERTL1(d2[0].num_elements() == nqtot,"Number of quadrature points do not match");

            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(2*m_coordim,1,0.0);
                
                if(m_coordim == 2) // assume g = [0,0,1]
                {
                    m_jac[0] = d1[0][0]*d2[1][0] - d2[0][0]*d1[1][0];

                    ASSERTL1(m_jac[0] > 0, "2D Regular Jacobian is not positive");
                    // Spencer's book page 160
                    m_gmat[0][0] =  d2[1][0]/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -d1[1][0]/m_jac[0]; // d xi_2/d x_1
                    m_gmat[2][0] = -d2[0][0]/m_jac[0]; // d xi_1/d x_2
                    m_gmat[3][0] =  d1[0][0]/m_jac[0]; // d xi_2/d x_2
                }
                else
                {
                    NekDouble g[3];
                    g[0] = d1[1][0]*d2[2][0] - d1[2][0]*d2[1][0];
                    g[1] = d1[2][0]*d2[0][0] - d1[0][0]*d2[2][0];
                    g[2] = d1[0][0]*d2[1][0] - d1[1][0]*d2[0][0];

                    m_jac[0] = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
                    ASSERTL1(m_jac[0] > 0, "Regular Jacobian is not positive");

                    m_gmat[0][0] =  (d2[1][0]*g[2] - d2[2][0]*g[1])/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -(d1[1][0]*g[2] - d1[2][0]*g[1])/m_jac[0]; // d xi_2/d x_1
                    m_gmat[2][0] = -(d2[0][0]*g[2] - d2[2][0]*g[0])/m_jac[0]; // d xi_1/d x_2
                    m_gmat[3][0] =  (d1[0][0]*g[2] - d1[2][0]*g[0])/m_jac[0]; // d xi_2/d x_2
                    m_gmat[4][0] =  (d2[0][0]*g[1] - d2[1][0]*g[0])/m_jac[0]; // d xi_1/d x_3
                    m_gmat[5][0] = -(d1[0][0]*g[1] - d1[1][0]*g[0])/m_jac[0]; // d xi_2/d x_3

                    m_jac[0] = sqrt(m_jac[0]);
                }
            }
            else
            {
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(2*m_coordim,nqtot,0.0);

                if(m_coordim == 2) // assume g = [0,0,1]
                {
                    // set up Jacobian 
                    Vmath::Vmul (nqtot,&d2[0][0],1,&d1[1][0],1,&m_jac[0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&m_jac[0],1,&m_jac[0],1);

                    ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0, "2D Deformed Jacobian is not positive");
                    
                    Vmath::Vdiv(nqtot,&d2[1][0],1,&m_jac[0],1,&m_gmat[0][0],1); // d xi_1/d x_1
                    Vmath::Vdiv(nqtot,&d1[1][0],1,&m_jac[0],1,&m_gmat[1][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[1][0],1);                          // d xi_2/d x_1
                    Vmath::Vdiv(nqtot,&d2[0][0],1,&m_jac[0],1,&m_gmat[2][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[2][0],1);                          // d xi_1/d x_2
                    Vmath::Vdiv(nqtot,&d1[0][0],1,&m_jac[0],1,&m_gmat[3][0],1); // d xi_2/d x_2
                }
                else
                {
                    Array<OneD,NekDouble> g[3] = {Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot)};
                    // g[0]
                    Vmath::Vmul (nqtot,&d1[2][0],1,&d2[1][0],1,&g[0][0],1);
                    Vmath::Vvtvm(nqtot,&d1[1][0],1,&d2[2][0],1,&g[0][0],1,&g[0][0],1);
                    //g[1]
                    Vmath::Vmul (nqtot,&d1[0][0],1,&d2[2][0],1,&g[1][0],1);
                    Vmath::Vvtvm(nqtot,&d1[2][0],1,&d2[0][0],1,&g[1][0],1,&g[1][0],1);
                    //g[2]
                    Vmath::Vmul (nqtot,&d1[1][0],1,&d2[0][0],1,&g[2][0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&g[2][0],1,&g[2][0],1);

                    // J_3D
                    Vmath::Vmul (nqtot,&g[0][0],1,&g[0][0],1,&m_jac[0],1);
                    Vmath::Vvtvp(nqtot,&g[1][0],1,&g[1][0],1,&m_jac[0],1,&m_jac[0],1);
                    Vmath::Vvtvp(nqtot,&g[2][0],1,&g[2][0],1,&m_jac[0],1,&m_jac[0],1);

                    // d xi_1/d x_1
                    Vmath::Vmul (nqtot,&d2[2][0],1,&g[1][0],1,&m_gmat[0][0],1);
                    Vmath::Vvtvm(nqtot,&d2[1][0],1,&g[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);

                    // d xi_2/d x_1
                    Vmath::Vmul (nqtot,&d1[1][0],1,&g[2][0],1,&m_gmat[1][0],1);
                    Vmath::Vvtvm(nqtot,&d1[2][0],1,&g[1][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);

                    // d xi_1/d x_2
                    Vmath::Vmul (nqtot,&d2[0][0],1,&g[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vvtvm(nqtot,&d2[2][0],1,&g[0][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);

                    // d xi_2/d x_2
                    Vmath::Vmul (nqtot,&d1[2][0],1,&g[0][0],1,&m_gmat[3][0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&g[2][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);

                    // d xi_1/d x_3
                    Vmath::Vmul (nqtot,&d2[1][0],1,&g[0][0],1,&m_gmat[4][0],1);
                    Vmath::Vvtvm(nqtot,&d2[0][0],1,&g[1][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);

                    // d xi_2/d x_3
                    Vmath::Vmul (nqtot,&d1[0][0],1,&g[1][0],1,&m_gmat[5][0],1);
                    Vmath::Vvtvm(nqtot,&d1[1][0],1,&g[0][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);

                    // J = sqrt(J_3D)
                    Vmath::Vsqrt(nqtot,&m_jac[0],1,&m_jac[0],1);
                }
            }
        }
 
        void GeomFactors::SetUpJacGmat3D(const Array<OneD, Array<OneD, NekDouble> > d1,
                                         const Array<OneD, Array<OneD, NekDouble> > d2,
                                         const Array<OneD, Array<OneD, NekDouble> > d3)
        {
            ASSERTL1(d1.num_elements()==m_coordim,"The dimension of array d1 does not"
                     "match the coordinate dimension");
            ASSERTL1(d2.num_elements()==m_coordim,"The dimension of array d2 does not"
                     "match the coordinate dimension");
            ASSERTL1(d3.num_elements()==m_coordim,"The dimension of array d3 does not"
                     "match the coordinate dimension");
            
            int nqtot = m_pointsKey[0].GetNumPoints() *
                        m_pointsKey[1].GetNumPoints() *
                        m_pointsKey[2].GetNumPoints();

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

            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(3*m_coordim,1,0.0);
                
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
                m_gmat = Array<TwoD, NekDouble>(3*m_coordim,nqtot,0.0);

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


        void GeomFactors::SetUpLaplacianMetrics2D(StdRegions::ExpansionType shape,
                                                  const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == m_expdim,"Inappropriate dimension of tbasis");
            ASSERTL1((m_coordim == 2)||(m_coordim <= 3),
                     "The coordinate dimension should be equal to two or three"
                     "for two-dimensional elements");

            int i;
            int nquad0 = m_pointsKey[0].GetNumPoints();
            int nquad1 = m_pointsKey[1].GetNumPoints();
            int nqtot  = nquad0*nquad1;

            m_laplacianmetrics      = Array<TwoD, NekDouble>(3,nqtot); 
            m_laplacianMetricIsZero = Array<OneD, bool>(3, false);
            m_laplacianMetricsFlag  = true;
 
            // Get hold of the quadrature weights
            const Array<OneD, const NekDouble>& w0 = tbasis[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = tbasis[1]->GetW();

            switch(shape)
            {
            case StdRegions::eQuadrilateral:
                {      
                    if(( m_gtype == eRegular)||
                       ( m_gtype == eMovingRegular))
                    {
                        NekDouble g0 = m_gmat[0][0]*m_gmat[0][0] + m_gmat[2][0]*m_gmat[2][0];
                        NekDouble g1 = m_gmat[0][0]*m_gmat[1][0] + m_gmat[2][0]*m_gmat[3][0];
                        NekDouble g2 = m_gmat[1][0]*m_gmat[1][0] + m_gmat[3][0]*m_gmat[3][0];

                        if(m_coordim == 3)
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

                        if(m_coordim == 3)
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
  
                    if(( m_gtype == eRegular)||
                       ( m_gtype == eMovingRegular))
                    {
                        Vmath::Smul (nqtot,m_gmat[0][0],&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Svtvp(nqtot,m_gmat[1][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);
                                             
                        Vmath::Vmul (nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Smul (nqtot,m_gmat[1][0],&tmp[0],1,&m_laplacianmetrics[1][0],1);


                        Vmath::Smul (nqtot,m_gmat[2][0],&dEta_dXi[0][0],1,&tmp[0],1);
                        Vmath::Svtvp(nqtot,m_gmat[3][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                        Vmath::Vvtvp(nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                        Vmath::Svtvp(nqtot,m_gmat[3][0],&tmp[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);

                        if(m_coordim == 3)
                        {
                            Vmath::Smul (nqtot,m_gmat[4][0],&dEta_dXi[0][0],1,&tmp[0],1);
                            Vmath::Svtvp(nqtot,m_gmat[5][0],&dEta_dXi[1][0],1,&tmp[0],1,&tmp[0],1);

                            Vmath::Vvtvp(nqtot,&tmp[0],1,   &tmp[0],1,&m_laplacianmetrics[0][0],1,&m_laplacianmetrics[0][0],1);
                            Vmath::Svtvp(nqtot,m_gmat[5][0],&tmp[0],1,&m_laplacianmetrics[1][0],1,&m_laplacianmetrics[1][0],1);
                        }

                        NekDouble g2 = m_gmat[1][0]*m_gmat[1][0] + m_gmat[3][0]*m_gmat[3][0];
                        if(m_coordim == 3)
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

                        if(m_coordim == 3)
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
            
        void GeomFactors::SetUpQuadratureMetrics2D(StdRegions::ExpansionType shape,
                                                   const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            ASSERTL1(tbasis.num_elements() == m_expdim,"Inappropriate dimension of tbasis");

            int i;
            int nquad0 = m_pointsKey[0].GetNumPoints();
            int nquad1 = m_pointsKey[1].GetNumPoints();
            int nqtot  = nquad0*nquad1;

            m_weightedjac           = Array<OneD, NekDouble>(nqtot); 
            m_quadratureMetricsFlag = true;

            // Fill the array m_weighted jac with the values
            // of the (already computed) jacobian (=m_jac)
            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
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
                    case LibUtilities::eGaussLobattoLegendre:
		    case LibUtilities::eGaussLobattoKronrodLegendre: // Legendre inner product 
                        for(i = 0; i < nquad1; ++i)
                        {
                            const Array<OneD, const NekDouble>& z1 = tbasis[1]->GetZ();
                            Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i],m_weightedjac.get()+i*nquad0,1);
                        }
                        break;
                    case LibUtilities::eGaussRadauMAlpha1Beta0:
		    case LibUtilities::eGaussRadauKronrodMAlpha1Beta0: // (1,0) Jacobi Inner product 
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

        // Post-process two tangential basis of 2D Manifold (differentiatino x_map in eta and xi directions )
        // to orthonormal vectors.

        void GeomFactors::SetUpTangentialbasis(const Array<OneD, Array<OneD,NekDouble> > d1_tbasis,
                                               const Array<OneD, Array<OneD,NekDouble> > d2_tbasis)
          {
              int coordim = d1_tbasis.num_elements();
              int nqtot = d1_tbasis[0].num_elements();
              //   m_coordim = coordim;

              // Initialization of tangential basis
	    m_tbasis1 = Array<OneD, Array<OneD, NekDouble> > (coordim);
	    m_tbasis2 = Array<OneD, Array<OneD, NekDouble> > (coordim);
            
            // Calculate local derivatives
            for(int i = 0; i < coordim; ++i)
            {
		m_tbasis1[i] = Array<OneD,NekDouble>(nqtot);
		m_tbasis2[i] = Array<OneD,NekDouble>(nqtot);
            }

	    // Assign them as global variables
	    for(int i = 0; i < coordim; ++i)
	      {
		Vmath::Vcopy(nqtot, &d1_tbasis[i][0], 1, &m_tbasis1[i][0], 1);
		Vmath::Vcopy(nqtot, &d2_tbasis[i][0], 1, &m_tbasis2[i][0], 1);
	      }

             // Normalization of each tangent vector
            Array<OneD, NekDouble> norm1(nqtot,0.0);
            Array<OneD, NekDouble> norm2(nqtot,0.0);
            Array<OneD, NekDouble> inner12(nqtot, 0.0);

            for (int i = 0; i < coordim; ++i)
            {
                 Vmath::Vvtvp(nqtot, m_tbasis1[i], 1, m_tbasis1[i], 1, norm1, 1, norm1, 1);
                 Vmath::Vvtvp(nqtot, m_tbasis2[i], 1, m_tbasis2[i], 1, norm2, 1, norm2, 1);
            }
            Vmath::Vsqrt(nqtot, norm1, 1, norm1, 1);
            Vmath::Vsqrt(nqtot, norm2, 1, norm2, 1);

            for (int i = 0; i < coordim; ++i)
            {
                Vmath::Vdiv(nqtot, m_tbasis1[i], 1, norm1, 1, m_tbasis1[i], 1);
                Vmath::Vdiv(nqtot, m_tbasis2[i], 1, norm2, 1, m_tbasis2[i], 1);
            }

            // Gram-Schmitz orthogonalization of two tangent vectors
            norm1 = Array<OneD, NekDouble>(nqtot, 0.0);
            norm2 = Array<OneD, NekDouble>(nqtot, 0.0);
            for (int i = 0; i < coordim; ++i)
             {
                     Vmath::Vvtvp(nqtot, m_tbasis1[i], 1, m_tbasis2[i], 1, inner12, 1, inner12, 1);
                     Vmath::Vvtvp(nqtot, m_tbasis2[i], 1, m_tbasis2[i], 1, norm2, 1, norm2, 1);
             }
             Vmath::Vdiv(nqtot, inner12, 1, norm2, 1, inner12, 1);
             Vmath::Neg(nqtot, inner12, 1);

             for (int i = 0; i < coordim; ++i)
              {
                  Vmath::Vvtvp(nqtot, inner12, 1, m_tbasis1[i], 1, m_tbasis2[i], 1, m_tbasis2[i], 1);
              }

              norm2 = Array<OneD, NekDouble>(nqtot, 0.0);
             for (int i = 0; i < coordim; ++i)
              {
                     Vmath::Vvtvp(nqtot, m_tbasis2[i], 1, m_tbasis2[i], 1, norm2, 1, norm2, 1);
              }
              Vmath::Vsqrt(nqtot, norm2, 1, norm2, 1);

             for (int i = 0; i < coordim; ++i)
              {
                     Vmath::Vdiv(nqtot, m_tbasis2[i], 1, norm2, 1, m_tbasis2[i], 1);
              }
          }


        // Generate Normal vectors at all quadature points specified
        // to the pointsKey "to_key" according to anticlockwise
        // convention 
        Array<OneD, NekDouble> GeomFactors::GenNormals2D(enum StdRegions::ExpansionType shape,
                                                         const int edge,  
                                                         const LibUtilities::PointsKey &to_key)
        {
            int i; 
            int nqe = to_key.GetNumPoints();
            Array<OneD, NekDouble> returnval(m_coordim*nqe,0.0);
            
            // Regular geometry case
            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch(shape)
                {
                case StdRegions::eTriangle:
                    switch(edge)
                    {
                    case 0:
                        for(i = 0; i < m_coordim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i+1][0],&returnval[i*nqe],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < m_coordim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i+1][0] + m_gmat[2*i][0],&returnval[i*nqe],1);
                        }
                            break;
                    case 2:
                        for(i = 0; i < m_coordim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i][0],&returnval[i*nqe],1);
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
                        for(i = 0; i < m_coordim; ++i)
                        {
                            Vmath::Fill(nqe,-m_gmat[2*i+1][0],&returnval[i*nqe],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < m_coordim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i][0],&returnval[i*nqe],1);
                        }
                        break;
                    case 2:
                        for(i = 0; i < m_coordim; ++i)
                        {
                            Vmath::Fill(nqe,m_gmat[2*i+1][0],&returnval[i*nqe],1);
                        }
                        break;
                    case 3:
                        for(i = 0; i < m_coordim; ++i)
                        {                            
                            Vmath::Fill(nqe,-m_gmat[2*i][0],&returnval[i*nqe],1);
                        }
                        break;
                    default:
                        ASSERTL0(false,"edge is out of range (edge < 4)");
                    }
                    break;
                }
                
                // normalise 
                fac = 0.0;
                for(i =0 ; i < m_coordim; ++i)
                {
                    fac += returnval[i*nqe]*returnval[i*nqe];
                }
                fac = 1.0/sqrt(fac);
                Vmath::Smul(m_coordim*nqe,fac,returnval,1,returnval,1);
            }
            else   // Set up deformed normals
            {
                int j;
                
                int nquad0 = m_pointsKey[0].GetNumPoints();
                int nquad1 = m_pointsKey[1].GetNumPoints();
                
               LibUtilities::PointsKey from_key;

                Array<OneD,NekDouble> normals(m_coordim*max(nquad0,nquad1),0.0);
                Array<OneD,NekDouble> jac    (m_coordim*max(nquad0,nquad1),0.0);

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
                                for(i = 0; i < m_coordim; ++i)
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
                                for(i = 0; i < m_coordim; ++i)
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
                                for(i = 0; i < m_coordim; ++i)
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
                                for(i = 0; i < m_coordim; ++i)
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
                                for(i = 0; i < m_coordim; ++i)
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
                                for(i = 0; i < m_coordim; ++i)
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
                                for(i = 0; i < m_coordim; ++i)
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
                for(i = 0; i < m_coordim; ++i)
                {
                    LibUtilities::Interp1D(from_key,&normals[i*nq],to_key,&returnval[i*nqe]);
                    Vmath::Vmul(nqe,&work[0],1,&returnval[i*nqe],1,&returnval[i*nqe],1);
                }
                
                //normalise normal vectors
                Vmath::Zero(nqe,work,1);
                for(i = 0; i < m_coordim; ++i)
                {
                    Vmath::Vvtvp(nqe,&returnval[i*nqe],1, &returnval[i*nqe],1,&work[0],1,&work[0],1);
                }

                Vmath::Vsqrt(nqe,&work[0],1,&work[0],1);
                Vmath::Sdiv(nqe,1.0,&work[0],1,&work[0],1);
                
                for(i = 0; i < m_coordim; ++i)
                {
                    Vmath::Vmul(nqe,&returnval[i*nqe],1,&work[0],1,&returnval[i*nqe],1);
                }                
            
                // Reverse direction so that points are in
                // anticlockwise direction if edge >=2
                if(edge >= 2) 
                {
                    for(i = 0; i < m_coordim; ++i)
                    {
                        Vmath::Reverse(nqe,&returnval[i*nqe],1, &returnval[i*nqe],1);
                    }
                }
            }
            
            return returnval;
        }
        
        GeomFactors::~GeomFactors(){
        }

        bool operator==(const GeomFactors &lhs, const GeomFactors &rhs)
        {
            if(!(lhs.m_gtype == rhs.m_gtype))
            {
                return false;
            }

            if(!(lhs.m_expdim == rhs.m_expdim))
            {
                return false;
            }

            if(!(lhs.m_coordim == rhs.m_coordim))
            {
                return false;
            }

            if(!(lhs.m_quadratureMetricsFlag == rhs.m_quadratureMetricsFlag))
            {
                return false;
            }

            if(!(lhs.m_laplacianMetricsFlag == rhs.m_laplacianMetricsFlag))
            {
                return false;
            }

            for(int i = 0; i < lhs.m_expdim; i++)
            {
                if(!(lhs.m_pointsKey[i] == rhs.m_pointsKey[i]))
                {
                    return false;
                }
            }
            
            if(!IsEqual(lhs.m_jac,rhs.m_jac,NekConstants::kGeomFactorsTol))
            {
                return false;
            }

            if(!IsEqual(lhs.m_gmat,rhs.m_gmat,NekConstants::kGeomFactorsTol))
            {
                return false;
            }

            return true;
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: GeomFactors.cpp,v $
// Revision 1.43  2009/06/15 01:59:21  claes
// Gauss-Kronrod updates
//
// Revision 1.42  2009/05/15 14:38:41  pvos
// Changed check for regular quads so that it also includes parallellograms
//
// Revision 1.41  2009/05/01 13:23:21  pvos
// Fixed various bugs
//
// Revision 1.40  2009/04/20 16:13:23  sherwin
// Modified Import and Write functions and redefined how Expansion is used
//
// Revision 1.39  2009/04/04 00:29:37  rcantao
// Made a few checks on J3D and companion. Remarks included.
//
// Revision 1.38  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.37  2009/01/01 02:33:29  ehan
// cleaned up the code
//
// Revision 1.36  2008/12/18 14:08:58  pvos
// NekConstants update
//
// Revision 1.35  2008/12/17 16:57:20  pvos
// Performance updates
//
// Revision 1.34  2008/12/16 14:09:07  pvos
// Performance updates
//
// Revision 1.33  2008/12/03 23:42:14  ehan
// Fixed some error for 3D Geomfactors.
//
// Revision 1.32  2008/11/25 23:02:15  ehan
// Fixed level 1 assertion violation: GeomFactors(..) only works for 1D and 2D, but not it works for 3D.
//
// Revision 1.31  2008/11/24 18:33:34  ehan
// Added 3D routines for GeomFactors()
//
// Revision 1.30  2008/09/23 22:07:32  ehan
// Modified 3D GeomFactor constructor
//
// Revision 1.29  2008/09/23 18:19:56  pvos
// Updates for working ProjectContField3D demo
//
// Revision 1.28  2008/09/20 11:32:24  ehan
// Rearranged 3D geomfactors
//
// Revision 1.27  2008/09/15 10:21:32  ehan
// Fixed some errors.
//
// Revision 1.26  2008/09/09 14:18:02  sherwin
// Removed m_normals from GeomFactor. Added GenNormals2D and additional copy type constructor
//
// Revision 1.25  2008/07/28 22:25:43  ehan
// Fixed error for deformed quads.
//
// Revision 1.24  2008/07/19 21:20:36  sherwin
// Changed normal orientation to anticlockwise
//
// Revision 1.23  2008/07/17 19:25:55  ehan
// Added 3D GeomFactors(..) .
//
// Revision 1.22  2008/07/09 11:40:25  sherwin
// Fixed the initialisation of m_expdim
//
// Revision 1.21  2008/06/13 18:07:30  ehan
// Commented out the function GeomFactors(..)
//
// Revision 1.20  2008/06/09 21:34:28  jfrazier
// Added some code for 3d.
//
// Revision 1.19  2008/05/30 00:33:48  delisi
// Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
//
// Revision 1.18  2008/04/06 22:32:53  bnelson
// Fixed gcc compiler warnings.
//
// Revision 1.17  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.16  2007/12/17 20:27:19  sherwin
// Added normals to GeomFactors
//
// Revision 1.15  2007/12/06 22:47:15  pvos
// 2D Helmholtz solver updates
//
// Revision 1.14  2007/12/03 21:30:35  sherwin
// Added normal details
//
// Revision 1.13  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.12  2007/07/13 09:02:24  sherwin
// Mods for Helmholtz solver
//
// Revision 1.11  2007/07/10 22:20:55  jfrazier
// Revision of geo fac manager to test for equality.
//
// Revision 1.10  2007/07/10 17:06:30  jfrazier
// Added method and underlying structure to manage geomfactors.
//
// Revision 1.9  2007/06/07 15:54:19  pvos
// Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
// Also made corrections to various ASSERTL2 calls
//
// Revision 1.8  2007/05/31 19:13:12  pvos
// Updated NodalTriExp + LocalRegions/Project2D + some other modifications
//
// Revision 1.7  2007/05/28 21:48:41  sherwin
// Update for 2D functionality
//
// Revision 1.6  2007/05/28 08:35:25  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.5  2007/05/25 17:52:01  jfrazier
// Updated to use new Array classes.
//
// Revision 1.4  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.3  2007/03/29 19:23:59  bnelson
// Replaced boost::shared_array with SharedArray
//
// Revision 1.2  2007/03/25 15:48:22  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.1  2007/03/20 09:17:39  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.5  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.4  2007/03/02 12:01:58  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.3  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.2  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
// Revision 1.18  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.17  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.16  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.15  2006/03/13 18:04:07  sherwin
//
// Corrected silly error in calling new
//
// Revision 1.14  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.13  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.12  2006/02/26 21:19:42  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.11  2006/02/19 01:37:32  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
