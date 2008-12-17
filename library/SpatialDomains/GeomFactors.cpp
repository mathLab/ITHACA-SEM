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
        GeomFactors::GeomFactors(void):m_gtype(eRegular),
                                       m_expdim(0), m_coordim(0)
        {
        }

        GeomFactors::GeomFactors(const GeomType gtype,
            const int expdim, const int coordim):m_gtype(gtype),
            m_expdim(expdim), m_coordim(coordim)
        {
            m_pointsKey = Array<OneD, LibUtilities::PointsKey>(m_expdim);
        }

        /** \brief One dimensional geometric factors based on one, two or three
        dimensional coordinate description 

        - Calculate the local factors \f$ d \xi/d x_i \f$ 
        as a recipricol derivative \f$ d x_i/d\xi \f$

        - Calcuate the jacobian as \f$ J = \sqrt{\frac{d x_1}{d _\xi}.
        \frac{d x_1}{d _\xi} + \frac{d x_2}{d _\xi}.\frac{d x_2}{d _\xi}  + 
        \frac{d x_3}{d _\xi}\frac{d x_3}{d _\xi}} \f$

        **/

        GeomFactors::GeomFactors(const GeomType gtype, const int coordim, 
            const Array<OneD, const StdRegions::StdExpansion1DSharedPtr>  &Coords):
            m_gtype(gtype), m_coordim(coordim), m_expdim(1)
        {
            int        i,nquad;
            LibUtilities::PointsType  ptype;

            ASSERTL1(coordim <= 3, "Only understand up to 3 Coordinate for this routine");

            nquad = Coords[0]->GetNumPoints(0);
            ptype = Coords[0]->GetPointsType(0);

            m_pointsKey = Array<OneD, LibUtilities::PointsKey> (1);

            m_pointsKey[0] = Coords[0]->GetBasis(0)->GetPointsKey();
            
            Array<OneD,NekDouble> der[3] = {Array<OneD, NekDouble>(nquad),
                                            Array<OneD, NekDouble>(nquad), 
                                            Array<OneD, NekDouble>(nquad)};

            // Calculate local derivatives using physical space storage
            for(i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetNumPoints(0)  == nquad,
                    "Points order are different for each coordinate");
                ASSERTL2(Coords[i]->GetPointsType(0) == ptype,
                    "Points type are different for each coordinate");
        
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(), Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(), der[i]);
            }

            SetUpJacGmat(nquad,der);
        }
        
        void GeomFactors::SetUpJacGmat(const int nquad,
                                       const Array<OneD,NekDouble> der[3])
        {
            int i;
            NekDouble fac0,fac1;

            if(( m_gtype == eRegular)||
               ( m_gtype == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(m_coordim,1,0.0);
                
                fac0 = fac1 = 0.0;
                for(i = 0; i < m_coordim; ++i)
                {
                    m_gmat[i][0] = (fabs(der[i][0]) > NekZeroTol)? 1.0/der[i][0]: 0.0;
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
                        m_gmat[i][j] = (fabs(der[i][j]) > NekZeroTol)? 1.0/der[i][j]: 0.0;
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
        GeomFactors::GeomFactors(const GeomType gtype, const int coordim,
                                 const Array<OneD, const StdRegions::StdExpansion2DSharedPtr> &Coords):
        m_gtype(gtype), m_coordim(coordim), m_expdim(2)
        {
            int nquad0,nquad1,nqtot;
            LibUtilities::PointsType  ptype0, ptype1;

            ASSERTL1((coordim >= 2)&&(coordim <= 3),
                "Only understand up to three Coordinate and must have "
                "at least two coordinates for this routine");

            m_pointsKey = Array<OneD, LibUtilities::PointsKey> (m_expdim);

            m_pointsKey[0] = Coords[0]->GetBasis(0)->GetPointsKey();
            nquad0 = m_pointsKey[0].GetNumPoints();
            ptype0 = m_pointsKey[0].GetPointsType();

            m_pointsKey[1] = Coords[0]->GetBasis(1)->GetPointsKey();
            nquad1 = m_pointsKey[1].GetNumPoints();
            ptype1 = m_pointsKey[1].GetPointsType();

            nqtot = nquad0*nquad1;

            // setup temp storage
            Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot)};

            // Calculate local derivatives using physical space storage
            for(int i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetNumPoints(0)  == nquad0,
                    "Points order are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetNumPoints(1)  == nquad1,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(0) == ptype0,
                    "Points type are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsType(1) == ptype1,
                    "Points type are different for coordinate 1 ");

                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),
                    Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1[i],d2[i]);
            }


            SetUpJacGmat(Coords[0]->DetExpansionType(),nquad0,nquad1,d1,d2);
        }

        GeomFactors::GeomFactors(const GeomType gtype, const int coordim,
                                 const Array<OneD, const StdRegions::StdExpansion3DSharedPtr> &Coords):
        m_gtype(gtype), m_coordim(coordim), m_expdim(3)
        {
            int nquad0,nquad1,nquad2,nqtot;
            LibUtilities::PointsType  ptype0, ptype1, ptype2;

            ASSERTL1((coordim >= 2)&&(coordim <= 3),
                "Only understand up to three Coordinate and must have "
                "at least two coordinates for this routine");

            m_pointsKey = Array<OneD, LibUtilities::PointsKey> (m_expdim);

            m_pointsKey[0] = Coords[0]->GetBasis(0)->GetPointsKey();
            nquad0 = m_pointsKey[0].GetNumPoints();
            ptype0 = m_pointsKey[0].GetPointsType();

            m_pointsKey[1] = Coords[0]->GetBasis(1)->GetPointsKey();
            nquad1 = m_pointsKey[1].GetNumPoints();
            ptype1 = m_pointsKey[1].GetPointsType();

            m_pointsKey[2] = Coords[0]->GetBasis(2)->GetPointsKey();
            nquad2 = m_pointsKey[2].GetNumPoints();
            ptype2 = m_pointsKey[2].GetPointsType();

            nqtot = nquad0*nquad1*nquad2;

            // setup temp storage
            Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d3[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot)};                                           

            // Calculate local derivatives using physical space storage
            for(int i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetNumPoints(0)  == nquad0,
                    "Points order are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetNumPoints(1)  == nquad1,
                    "Points order are different for coordinate 1 ");
               ASSERTL2(Coords[i]->GetNumPoints(2)  == nquad2,
                    "Points order are different for coordinate 2 ");
                ASSERTL2(Coords[i]->GetPointsType(0) == ptype0,
                    "Points type are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsType(1) == ptype1,
                    "Points type are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(2) == ptype2,
                    "Points type are different for coordinate 2 ");

                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(),Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1[i],d2[i],d3[i]);

//                 for( int j = 0; j < nqtot; ++j ) {
//                     cout << "d1[" << i << "]["<< j  << "]= " << d1[i][j] << endl;
//                     cout << "d2[" << i << "]["<< j  << "]= " << d2[i][j] << endl;
//                     cout << "d3[" << i << "]["<< j  << "]= " << d3[i][j] << endl;
//                 }
            }

            cout << "D_1[1:" << coordim << "][1:" << nqtot << "] = " << endl;
            for( int i = 0; i < coordim; ++i ) {
                for( int j = 0; j < nqtot; ++j ) {
                    cout << d1[i][j] << "\t\t";
                }
                cout << endl;
            }

            cout << "\nD_2[1:" << coordim << "][1:" << nqtot << "] = " << endl;
            for( int i = 0; i < coordim; ++i ) {
                for( int j = 0; j < nqtot; ++j ) {
                    cout << d2[i][j] << "\t\t";
                }
                cout << endl;
            }

            cout << "\nD_3[1:" << coordim << "][1:" << nqtot << "] = " << endl;
            for( int i = 0; i < coordim; ++i ) {
                for( int j = 0; j < nqtot; ++j ) {
                    cout << d3[i][j] << "\t\t";
                }
                cout << endl;
            }

            SetUpJacGmat(Coords[0]->DetExpansionType(),nquad0,nquad1,nquad2,d1,d2,d3);
        }  

        GeomFactors::GeomFactors(enum StdRegions::ExpansionType shape,
                                 const GeomFactors &Xgfac,
                                 const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {

            int nq0,nq1,nq2,nqtot;
            int Xnq0, Xnq1,Xnq2, Xnqtot;
            LibUtilities::PointsKey fpoints0;
            LibUtilities::PointsKey fpoints1;
            LibUtilities::PointsKey fpoints2;

                m_gtype   = Xgfac.m_gtype; 
                m_expdim  = Xgfac.m_expdim;
                m_coordim = Xgfac.m_coordim;

                m_pointsKey = Array<OneD,LibUtilities::PointsKey>(m_expdim);

                ASSERTL0(m_gtype == eDeformed,
                        "This routine assumes element geometry is deformed");

            if(m_coordim == 2)
            {
             
                fpoints0 = Xgfac.m_pointsKey[0];
                fpoints1 = Xgfac.m_pointsKey[1];

                Xnq0 = fpoints0.GetNumPoints();
                Xnq1 = fpoints1.GetNumPoints();

                Xnqtot = Xnq0*Xnq1;

                ASSERTL1((m_coordim >= 2)&&(m_coordim <= 3),
                        "Only understand up to three Coordinate and must have "
                        "at least two coordinates for this routine");

                m_pointsKey[0] = tbasis[0]->GetPointsKey();
                m_pointsKey[1] = tbasis[1]->GetPointsKey();


                nq0 = m_pointsKey[0].GetNumPoints();
                nq1 = m_pointsKey[1].GetNumPoints();

                
                nqtot = nq0*nq1;
                
                // setup temp storage
                Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                            Array<OneD, NekDouble>(nqtot), 
                                            Array<OneD, NekDouble>(nqtot)};
                Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                            Array<OneD, NekDouble>(nqtot), 
                                            Array<OneD, NekDouble>(nqtot)};

                                                                         
                // interpolate Geometric derivatives which are polynomials
                Array<OneD,NekDouble> dxdxi(Xnqtot);
            
//             if(m_coordim == 2)
//             {
                // Interpolate d2[1];
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[0])[0],1,&dxdxi[0],1);
                LibUtilities::Interp2D(fpoints0, fpoints1, dxdxi, m_pointsKey[0], m_pointsKey[1], d2[1]);
                
                // Interpolate d1[1];
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[1])[0],1,&dxdxi[0],1);
                LibUtilities::Interp2D(fpoints0, fpoints1, dxdxi,m_pointsKey[0], m_pointsKey[1], d1[1]);
                Vmath::Neg(nqtot,d1[1],1);
                
                
                // Interpolate d2[0];
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[2])[0],1,&dxdxi[0],1);
                LibUtilities::Interp2D(fpoints0, fpoints1, dxdxi,m_pointsKey[0], m_pointsKey[1], d2[0]);
                Vmath::Neg(nqtot,d2[0],1);
                
                        // Interpolate d1[0];
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[3])[0],1,&dxdxi[0],1);
                LibUtilities::Interp2D(fpoints0, fpoints1, dxdxi, m_pointsKey[0], m_pointsKey[1], d1[0]);
                
                SetUpJacGmat(shape, nq0,nq1,d1,d2);
                
            }
            else if(m_coordim == 3)
            {

                fpoints0 = Xgfac.m_pointsKey[0];
                fpoints1 = Xgfac.m_pointsKey[1];
                fpoints2 = Xgfac.m_pointsKey[2];

                

                Xnq0 = fpoints0.GetNumPoints();
                Xnq1 = fpoints1.GetNumPoints();
                Xnq2 = fpoints2.GetNumPoints();
                Xnqtot = Xnq0*Xnq1*Xnq2;

                ASSERTL1((m_coordim >= 2)&&(m_coordim <= 3),
                        "Only understand up to three Coordinate and must have "
                        "at least two coordinates for this routine");

                m_pointsKey[0] = tbasis[0]->GetPointsKey();
                m_pointsKey[1] = tbasis[1]->GetPointsKey();
                m_pointsKey[2] = tbasis[2]->GetPointsKey();

                nq0 = m_pointsKey[0].GetNumPoints();
                nq1 = m_pointsKey[1].GetNumPoints();
                nq2 = m_pointsKey[2].GetNumPoints();

                cout << "nq0 = " << nq0 << endl;
                cout << "nq1 = " << nq1 << endl;
                cout << "nq2 = " << nq2 << endl;
                
                nqtot = nq0*nq1*nq2;
                
                // setup temp storage
                Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                            Array<OneD, NekDouble>(nqtot), 
                                            Array<OneD, NekDouble>(nqtot)};
                Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                            Array<OneD, NekDouble>(nqtot), 
                                            Array<OneD, NekDouble>(nqtot)};
                Array<OneD,NekDouble> d3[3] = {Array<OneD, NekDouble>(nqtot),
                                            Array<OneD, NekDouble>(nqtot), 
                                            Array<OneD, NekDouble>(nqtot)};
                                            
                               
                // interpolate Geometric derivatives which are polynomials
                Array<OneD,NekDouble> dxdxi(Xnqtot);                   

            
               // Implementation based on Spend's book page 160

                // Compute the RHS of the d xi_1/d x_1 = m_gmat[0] : J3D * m_gmat[0] = dxdi[0], then interpolate d3[2]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[0])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints1, dxdxi, m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d3[2]);

                 // Compute the RHS of the d xi_1/d x_2 = m_gmat[1] : J3D * m_gmat[1] = dxdi[0], then interpolate d2[2]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[1])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d2[2]);
                //Vmath::Neg(nqtot,d2[2],1);

                // Compute the RHS of the d xi_1/d x_3 = m_gmat[2] : J3D * m_gmat[2] = dxdi[0], then interpolate d1[2]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[2])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d1[2]);
  
                // Compute the RHS of the d xi_2/d x_1 = m_gmat[3] : J3D * m_gmat[3] = dxdi[0], then interpolate d3[1]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[3])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d3[1]);
                //Vmath::Neg(nqtot,d3[1],1);

                // Compute the RHS of the d xi_2/d x_2 = m_gmat[4] : J3D * m_gmat[4] = dxdi[0], then interpolate d2[1]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[4])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d2[1]);

               // Compute the RHS of the d xi_2/d x_3 = m_gmat[5] : J3D * m_gmat[5] = dxdi[0], then interpolate d1[1]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[5])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d1[1]);
               // Vmath::Neg(nqtot,d1[1],1);

                // Compute the RHS of the d xi_3/d x_1 = m_gmat[6] : J3D * m_gmat[6] = dxdi[0], then interpolate d3[0]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[6])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d3[0]);

                // Compute the RHS of the d xi_3/d x_2 = m_gmat[7] : J3D * m_gmat[7] = dxdi[0], then interpolate d2[0]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[7])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d2[0]);
               // Vmath::Neg(nqtot,d2[0],1);

                 // Compute the RHS of the d xi_3/d x_3 = m_gmat[8] : J3D * m_gmat[8] = dxdi[0], then interpolate d1[0]
                Vmath::Vmul(Xnqtot,&(Xgfac.m_jac)[0],1,&(Xgfac.m_gmat[8])[0],1,&dxdxi[0],1);
                LibUtilities::Interp3D(fpoints0, fpoints1,fpoints2, dxdxi,m_pointsKey[0], m_pointsKey[1],m_pointsKey[2], d1[0]);
              
                SetUpJacGmat(shape, nq0,nq1,nq2,d1,d2,d3);
                
            } else
            {
                ASSERTL0(false,"Invalid dimension");
            }
            
            //SetUpJacGmat(shape, nq0,nq1,d1,d2);

        }
        
        void GeomFactors::SetUpJacGmat(enum StdRegions::ExpansionType shape,
                                       const int nquad0, 
                                       const int nquad1, 
                                       const Array<OneD, NekDouble> d1[3],
                                       const Array<OneD, NekDouble> d2[3])
        {
            int i,j,nqtot;

            nqtot = nquad0*nquad1;

            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(2*m_coordim,1,0.0);
                
                if(m_coordim == 2) // assume g = [0,0,1]
                {
                    m_jac[0] = d1[0][0]*d2[1][0] - d2[0][0]*d1[1][0];

                    ASSERTL1(m_jac[0] > 0, "2D Regular Jacobian is not positive");
                    // Spen's book page 160
                    m_gmat[0][0] =  d2[1][0]/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -d1[1][0]/m_jac[0]; // d xi_1/d x_2
                    m_gmat[2][0] = -d2[0][0]/m_jac[0]; // d xi_2/d x_1
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
                    m_gmat[1][0] = -(d1[1][0]*g[2] - d1[2][0]*g[1])/m_jac[0]; // d xi_1/d x_2
                    m_gmat[2][0] = -(d2[0][0]*g[2] - d2[2][0]*g[0])/m_jac[0]; // d xi_2/d x_1
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
                    Vmath::Neg(nqtot,&m_gmat[1][0],1);                          // d xi_1/d x_2
                    Vmath::Vdiv(nqtot,&d2[0][0],1,&m_jac[0],1,&m_gmat[2][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[2][0],1);                          // d xi_2/d x_1
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

        void GeomFactors::SetUpJacGmat(enum StdRegions::ExpansionType shape,
                                       const int nquad0, 
                                       const int nquad1,
                                       const int nquad2,
                                       const Array<OneD, NekDouble> d1[3],
                                       const Array<OneD, NekDouble> d2[3],
                                       const Array<OneD, NekDouble> d3[3])
        {
            
            int nqtot = nquad0*nquad1*nquad2;

            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
            {
                m_jac     = Array<OneD, NekDouble>(1,0.0);
                m_gmat    = Array<TwoD, NekDouble>(3*m_coordim,1,0.0);
                
                if(m_coordim == 3) // assume g = [0,0,1]
                {
                    // J3D: Determinent of three-dimensional Jacobian
                    m_jac[0] =  d1[0][0]*(d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0])
                               -d1[1][0]*(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0])
                               +d1[2][0]*(d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0]);
                               
                                cout << "J3D1 = " << m_jac[0] << endl;

                                cout << "d1[0] = " << d1[0][0] << endl;
                                cout << "d1[1] = " << d1[1][0] << endl;
                                cout << "d1[2] = " << d1[2][0] << endl;
                                cout << "d2[0] = " << d2[0][0] << endl;
                                cout << "d2[1] = " << d2[1][0] << endl;
                                cout << "d2[2] = " << d2[2][0] << endl;
                                cout << "d3[0] = " << d3[0][0] << endl;
                                cout << "d3[1] = " << d3[1][0] << endl;
                                cout << "d3[2] = " << d3[2][0] << endl;

                                cout << " d1[0][0]*(d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0]) = " <<  d1[0][0]*(d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0]) << endl;

                                cout << " d1[1][0]*(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0]) = " << d1[1][0]*(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0]) << endl;

                                cout << " d1[2][0]*(d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0]) = " << d1[2][0]*(d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0]) << endl;

                    ASSERTL1(m_jac[0] > 0, "3D Regular Jacobian is not positive");
                    // Spen's book page 160
                    m_gmat[0][0] =  (d2[1][0]*d3[2][0] - d2[2][0]*d3[1][0])/m_jac[0];  // d xi_1/d x_1
                    m_gmat[1][0] = -(d1[1][0]*d3[2][0] - d1[2][0]*d3[1][0])/m_jac[0];  // d xi_1/d x_2
                    m_gmat[2][0] =  (d1[1][0]*d2[2][0] - d1[2][0]*d2[1][0])/m_jac[0];  // d xi_1/d x_3
                    m_gmat[3][0] = -(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0])/m_jac[0];  // d xi_2/d x_1
                    m_gmat[4][0] =  (d1[0][0]*d3[2][0] - d1[2][0]*d3[0][0])/m_jac[0];  // d xi_2/d x_2
                    m_gmat[5][0] = -(d1[0][0]*d2[2][0] - d1[2][0]*d2[0][0])/m_jac[0];  // d xi_2/d x_3
                    m_gmat[6][0] =  (d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0])/m_jac[0];  // d xi_3/d x_1
                    m_gmat[7][0] = -(d1[0][0]*d3[1][0] - d1[1][0]*d3[0][0])/m_jac[0];  // d xi_3/d x_2
                    m_gmat[8][0] =  (d1[0][0]*d2[1][0] - d1[1][0]*d2[0][0])/m_jac[0];  // d xi_3/d x_3

                }
                else
                {
                    ASSERTL0(false,"Routine not yet implemented");
                }

            }
            else // Deformed case
            {
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(3*m_coordim,nqtot,0.0);

                if(m_coordim == 3) // assume g = [0,0,1]
                {
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
                    cout << "J3D1 = " << m_jac[0] << endl;
                    Vmath::Vvtvp(nqtot,&d1[1][0],1,&tmp[1][0],1,&m_jac[0],1,&m_jac[0],1);
                    cout << "J3D2 = " << m_jac[0] << endl;
                    Vmath::Vvtvp(nqtot,&d1[2][0],1,&tmp[2][0],1,&m_jac[0],1,&m_jac[0],1);
                    cout << "J3D3 = " << m_jac[0] << endl;


                   
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
                else
                {
                    ASSERTL0(false,"Routine not yet implemented");

                }
            }
        }

        // Generate Normal vectors at all quadature points specified
        // to the pointsKey "to_key" according to anticlockwise
        // convention 
        Array<OneD, NekDouble> GeomFactors::GenNormals2D(enum StdRegions::ExpansionType shape, const int edge,  const LibUtilities::PointsKey &to_key)
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



        /**
        \brief Three dimensional geometric factors and Jacobian

        The geometric factors are related to the local derivatives of the
        global coordinates by:

        \f$ \begin{array}{lll}
        \frac{\partial \xi_1}{\partial x_1} = \frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_3} - \frac{\partial x_2}{\partial \xi_3} \frac{\partial x_3}{\partial \xi_2} \right ) & 
        \frac{\partial \xi_1}{\partial x_2} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_3}  - \frac{\partial x_3}{\partial \xi_2} \frac{\partial x_1}{\partial \xi_3}  \right )&
        \frac{\partial \xi_1}{\partial x_3} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_2}{\partial \xi_3}  - \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_1}{\partial \xi_3}  \right )\\

        \frac{\partial \xi_2}{\partial x_1} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_3} - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_3}  \right ) &
        \frac{\partial \xi_2}{\partial x_2} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_3}  - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_3} \right ) &
        \frac{\partial \xi_2}{\partial x_3} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_3}  - \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_3}  \right )\\

        \frac{\partial \xi_3}{\partial x_1} = \frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2} - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_3}  \right ) &
        \frac{\partial \xi_3}{\partial x_2} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2}  - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_2} \right ) &
        \frac{\partial \xi_3}{\partial x_3} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_2}  - \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_2}  \right )
        \end{array} \f$

        where

        \f$ J_{3D} = 
        \frac{\partial x_1}{\partial \xi_1} \left ( 
        \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_3} - 
        \frac{\partial x_2}{\partial \xi_3} \frac{\partial x_3}{\partial \xi_2}\right )
        - \frac{\partial x_1}{\partial \xi_2} \left (  
        \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_3}  -
        \frac{\partial x_2}{\partial \xi_3} \frac{\partial x_3}{\partial \xi_1}\right )
        + \frac{\partial x_1}{\partial \xi_3} \left ( 
        \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2} -
        \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_1} 
        \right ) \f$

        **/
        
//         GeomFactors::GeomFactorsOld(const GeomType gtype, const int coordim,
//                                  const Array<OneD, const StdRegions::StdExpansion3DSharedPtr> &Coords):
//         m_gtype(gtype), m_coordim(coordim), m_expdim(3)
//         {
//             ASSERTL1((coordim<=3), "Only understand up to three coordinate");
// 
//             m_pointsKey = Array<OneD, LibUtilities::PointsKey> (m_expdim);
// 
//             int nquad0, nquad1, nquad2, nqtot; 
//             LibUtilities::PointsType  ptype0, ptype1, ptype2;
// 
//             m_pointsKey[0] = Coords[0]->GetBasis(0)->GetPointsKey();
//             
//             nquad0 = m_pointsKey[0].GetNumPoints();
//             ptype0 = m_pointsKey[0].GetPointsType();
//               
//             m_pointsKey[1] = Coords[0]->GetBasis(1)->GetPointsKey();
//             nquad1 = m_pointsKey[1].GetNumPoints();
//             ptype1 = m_pointsKey[1].GetPointsType();
//                
//             m_pointsKey[2] = Coords[0]->GetBasis(2)->GetPointsKey();;
//             nquad2 = m_pointsKey[2].GetNumPoints();
//             ptype2 = m_pointsKey[2].GetPointsType();
//              
//  
//              nqtot = nquad0*nquad1*nquad2;
// 
//              cout << "nqtot = " << nqtot << endl;
// 
//             // setup temp storage
//             Array<OneD, NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
//                                             Array<OneD, NekDouble>(nqtot),
//                                             Array<OneD, NekDouble>(nqtot)};
//             Array<OneD, NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
//                                             Array<OneD, NekDouble>(nqtot),
//                                             Array<OneD, NekDouble>(nqtot)};
//             Array<OneD, NekDouble> d3[3] = {Array<OneD, NekDouble>(nqtot),
//                                             Array<OneD, NekDouble>(nqtot),
//                                             Array<OneD, NekDouble>(nqtot)};
// 
//         // Calculate local derivatives using physical space storage
//         for(int i=0; i<coordim; ++i)
//         {
//             ASSERTL2(Coords[i]->GetNumPoints(0) == nquad0,
//                     "Points order are different for coordinate 0");
//             ASSERTL2(Coords[i]->GetNumPoints(1) == nquad1,
//                     "Points order are different for coordinate 1");
//             ASSERTL2(Coords[i]->GetNumPoints(2) == nquad2,
//                     "Points order are different for coordinate 2");
//            
//             ASSERTL2(Coords[i]->GetPointsType(0) == ptype0,
//                     "Points type are different for coordinate 0");
//             ASSERTL2(Coords[i]->GetPointsType(1) == ptype1,
//                     "Points type are different for coordinate 1");
//             ASSERTL2(Coords[i]->GetPointsType(2) == ptype2,
//                     "Points type are different for coordinate 2");
// 
//             Coords[i]->BwdTrans(Coords[i]->GetCoeffs(), Coords[i]->UpdatePhys());
//             Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(), d1[i], d2[i], d3[i]);
// 
//         }
// 
//         if((m_gtype == eRegular) || (m_gtype == eMovingRegular) )
//         {
//             m_jac = Array<OneD, NekDouble>(1,0.0);
//             m_gmat = Array<TwoD, NekDouble>(3*coordim, 1, 0.0);
//             
//                 // The three-dimensional Jacobian form (J_3d from the page 158, Spen's book)
//                 m_jac[0] = d1[0][0]*(d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])
//                          - d2[0][0]*(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0])
//                          + d3[0][0]*(d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0]);
//                 
//                 ASSERTL1(m_jac[0] > 0, "3D Regular Jacobian is not positive");
// 
//                 // Partial derivatives with respect to x_1, x_2, and x_3
//                 // in terms of derivative with respect to xi_1, xi_2, and x_3.
//                 // (Spen's book page 160)
//                 m_gmat[0][0] =  (d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])/m_jac[0];// d xi_1/d x_1
//                 m_gmat[1][0] = -(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0])/m_jac[0];// d xi_2/d x_1
//                 m_gmat[2][0] =  (d1[1][0]*d2[2][0] - d3[1][0]*d1[2][0])/m_jac[0];// d xi_3/d x_1
// 
//                 m_gmat[3][0] =  -(d2[0][0]*d3[2][0] - d3[0][0]*d2[2][0])/m_jac[0];// d xi_1/d x_2
//                 m_gmat[4][0] =   (d1[0][0]*d3[2][0] - d3[0][0]*d1[2][0])/m_jac[0];// d xi_2/d x_2
//                 m_gmat[5][0] =  -(d1[0][0]*d2[2][0] - d2[0][0]*d1[2][0])/m_jac[0];// d xi_3/d x_2
// 
//                 m_gmat[6][0] =   (d2[0][0]*d3[1][0] - d3[0][0]*d2[1][0])/m_jac[0];// d xi_1/d x_3
//                 m_gmat[7][0] =  -(d1[0][0]*d3[1][0] - d3[0][0]*d1[1][0])/m_jac[0];// d xi_2/d x_3
//                 m_gmat[8][0] =   (d1[0][0]*d2[1][0] - d2[0][0]*d1[1][0])/m_jac[0];// d xi_3/d x_3
// 
//                 m_jac[0] = sqrt(m_jac[0]);
//                                               
//             }
//             else
//             {
//                 m_jac = Array<OneD, NekDouble>(nqtot, 0.0);
//                 m_gmat = Array<TwoD, NekDouble>(3*coordim, nqtot, 0.0);
//                 
//                 Array<OneD, NekDouble> g[3] = {Array<OneD, NekDouble>(nqtot),
//                                                Array<OneD, NekDouble>(nqtot),
//                                                Array<OneD, NekDouble>(nqtot)};
// 
//                 cout << "*********Geom type = " <<m_gtype << endl;                                             
//                                                 
//                 // set up Jacobian
//                 // g[0] = (d x_2/d xi_2)*(d x_3/d xi_3) - (d x_2/d xi_3)*(d x_3/d xi_2)
//                 //g[0] = (d x_2/d xi_3)*(d x_3/d xi_2)
//                 Vmath::Vmul(nqtot,&d3[1][0],1,&d2[2][0],1,&g[0][0],1);
//                 //g[0] = (d x_2/d xi_2)*(d x_3/d xi_3) - g[0]
//                 Vmath::Vvtvm(nqtot,&d2[1][0],1,&d3[2][0],1,&g[0][0],1,&g[0][0],1);
// 
//                 // g[1] = (d x_2/d xi_1)*(d x_3/d xi_3) - (d x_2/d xi_3)*(d x_3/d xi_1)
//                 //g[1] = (d x_2/d xi_3)*(d x_3/d xi_1)
//                 Vmath::Vmul(nqtot,&d3[1][0],1,&d1[2][0],1,&g[1][0],1);
//                 //g[1] = (d x_2/d xi_1)*(d x_3/d xi_3) - g[1]
//                 Vmath::Vvtvm(nqtot,&d1[1][0],1,&d3[2][0],1,&g[1][0],1,&g[1][0],1);
// 
//                 // g[2] = (d x_2/d xi_1)*(d x_3/d xi_2) - (d x_2/d xi_2)*(d x_3/d xi_1)
//                 //g[2] = (d x_2/d xi_2)*(d x_3/d xi_1)
//                 Vmath::Vmul(nqtot,&d2[1][0],1,&d1[2][0],1,&g[2][0],1);
//                 //g[2] = (d x_2/d xi_1)*(d x_3/d xi_2) - g[2]
//                 Vmath::Vvtvm(nqtot,&d1[1][0],1,&d2[2][0],1,&g[2][0],1,&g[2][0],1);
// 
//                 // J_3D 
//                 Vmath::Vmul(nqtot,&d2[0][0],1,&g[1][0],1,&m_jac[0],1);
//                 Vmath::Vvtvm(nqtot,&d1[0][0],1,&g[0][0],1,&m_jac[0],1,&m_jac[0],1);
//                 Vmath::Vvtvp(nqtot,&d3[0][0],1,&g[2][0],1,&m_jac[0],1,&m_jac[0],1);
// 
//                 // d xi_1/d x_1
//                 Vmath::Vmul(nqtot,&d3[1][0],1,&d2[2][0],1,&m_gmat[0][0],1);
//                 Vmath::Vvtvm(nqtot,&d2[1][0],1,&d3[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);
// 
//                 // d xi_2/d x_1
//                 Vmath::Vmul(nqtot,&d1[1][0],1,&d3[2][0],1,&m_gmat[1][0],1);
//                 Vmath::Vvtvm(nqtot,&d3[1][0],1,&d1[2][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);
// 
//                 // d xi_3/d x_1
//                 Vmath::Vmul(nqtot,&d3[1][0],1,&d1[2][0],1,&m_gmat[2][0],1);
//                 Vmath::Vvtvm(nqtot,&d1[1][0],1,&d2[2][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);
// 
//                 // d xi_1/d x_2
//                 Vmath::Vmul(nqtot,&d2[0][0],1,&d3[2][0],1,&m_gmat[3][0],1);
//                 Vmath::Vvtvm(nqtot,&d3[0][0],1,&d2[2][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);
//                 
//                 // d xi_2/d x_2
//                 Vmath::Vmul(nqtot,&d3[0][0],1,&d1[2][0],1,&m_gmat[4][0],1);
//                 Vmath::Vvtvm(nqtot,&d1[0][0],1,&d3[2][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);
//                 
//                 // d xi_3/d x_2
//                 Vmath::Vmul(nqtot,&d1[0][0],1,&d2[2][0],1,&m_gmat[5][0],1);
//                 Vmath::Vvtvm(nqtot,&d2[0][0],1,&d1[2][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);
// 
//                 // d xi_1/d x_3
//                 Vmath::Vmul(nqtot,&d3[0][0],1,&d2[1][0],1,&m_gmat[6][0],1);
//                 Vmath::Vvtvm(nqtot,&d2[0][0],1,&d3[1][0],1,&m_gmat[6][0],1,&m_gmat[6][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[6][0],1,&m_jac[0],1,&m_gmat[6][0],1);
//                 
//                 // d xi_2/d x_3
//                 Vmath::Vmul(nqtot,&d1[0][0],1,&d3[1][0],1,&m_gmat[7][0],1);
//                 Vmath::Vvtvm(nqtot,&d3[0][0],1,&d1[1][0],1,&m_gmat[7][0],1,&m_gmat[7][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[7][0],1,&m_jac[0],1,&m_gmat[7][0],1);
//                 
//                 // d xi_3/d x_3
//                 Vmath::Vmul(nqtot,&d2[0][0],1,&d1[1][0],1,&m_gmat[8][0],1);
//                 Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&m_gmat[8][0],1,&m_gmat[8][0],1);
//                 Vmath::Vdiv(nqtot,&m_gmat[8][0],1,&m_jac[0],1,&m_gmat[8][0],1);
// 
//                 // J = sqrt(J_3D)
//                 Vmath::Vsqrt(nqtot,&m_jac[0],1,&m_jac[0],1);
// 
//             }
// 
//      }
        
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

            if(!(lhs.m_jac == rhs.m_jac))
            {
                return false;
            }

            if(!(lhs.m_gmat == rhs.m_gmat))
            {
                return false;
            }

            return true;
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: GeomFactors.cpp,v $
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
