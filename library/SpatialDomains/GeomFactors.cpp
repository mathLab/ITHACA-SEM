////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/GeomFactors.cpp,v $
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

namespace Nektar
{
    namespace SpatialDomains
    {
        GeomFactors::GeomFactors(void):m_gtype(eNoGeomType)
        {
        }

        GeomFactors::GeomFactors(const GeomType gtype,
            const int expdim, const int coordim):m_gtype(gtype)
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

        GeomFactors::GeomFactors(const GeomType gtype, const int coordim, 
				 const ConstArray<OneD,StdRegions::StdExpansion1DSharedPtr>  &Coords)
        {
            int        i,nquad;
            LibUtilities::PointsType  ptype;

            ASSERTL1(coordim <= 3, "Only understand up to 3 Coordinate for this routine");

            m_gtype  = gtype;

            nquad = Coords[0]->GetNumPoints(0);
            ptype = Coords[0]->GetPointsType(0);

            Array<OneD,NekDouble> der[3] = {Array<OneD, NekDouble>(nquad),
                                            Array<OneD, NekDouble>(nquad), 
                                            Array<OneD, NekDouble>(nquad)};

            // Calculate local derivatives using physical space storage
            for(i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetPointsOrder(0) == nquad,
                    "Points order are different for each coordinate");
                ASSERTL2(Coords[i]->GetPointsType(0)  == ptype,
                    "Points type are different for each coordinate");
		
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(), 
				    Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(), der[i]);
            }

            if((m_gtype == eRegular)||
               (m_gtype == eMovingRegular))
            {
                m_jac  = Array<OneD, NekDouble>(1,0.0);
                m_gmat = Array<TwoD, NekDouble>(coordim,1,0.0);

                for(i = 0; i < coordim; ++i)
                {
                    m_gmat[i][0] = 1.0/der[i][0];
                    m_jac[0] += der[i][0]*der[i][0];
                }

                m_jac[0] = sqrt(m_jac[0]);
            }
            else
            {
                m_jac  = Array<OneD, NekDouble>(nquad);
                m_gmat = Array<TwoD, NekDouble>(coordim,nquad);

                // invert local derivative for gmat;
                for(i = 0; i < coordim; ++i)
                {
                    Vmath::Sdiv(nquad,1.0,&der[i][0],1,&m_gmat[i][0],1);
                    Vmath::Vmul(nquad,&der[i][0],1,&der[i][0],1,&m_jac[0],1);
                }
                Vmath::Vsqrt(nquad,&m_jac[0],1,&m_jac[0],1);
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
                                 const ConstArray<OneD,StdRegions::StdExpansion2DSharedPtr> &Coords)
        {
            int           i,nquad0,nquad1,nqtot;
            LibUtilities::PointsType  ptype0, ptype1;

            ASSERTL1((coordim >= 2)&&(coordim <= 3),
                "Only understand up to three Coordinate and must have "
                "at least two coordinates for this routine");

            m_gtype = gtype;

            nquad0 = Coords[0]->GetNumPoints(0);
            nquad1 = Coords[0]->GetNumPoints(1);
            ptype0 = Coords[0]->GetPointsType (0);
            ptype1 = Coords[0]->GetPointsType (1);

            nqtot = nquad0*nquad1;

            // setup temp storage
            Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};

            // Calculate local derivatives using physical space storage
            for(i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetPointsOrder(0) == nquad0,
                    "Points order are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsOrder(1) == nquad1,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(0)  == ptype0,
                    "Points type are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointstTpe(1)  == ptype1,
                    "Points type are different for coordinate 1 ");

                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(), 
				    Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1[i],d2[i]);
            }

            if((m_gtype == eRegular)||
               (m_gtype == eMovingRegular))
            {
                m_jac  = Array<OneD, NekDouble>(1,0.0);
                m_gmat = Array<TwoD, NekDouble>(2*coordim,1,0.0);

                if(coordim == 2) // assume g = [0,0,1]
                {
                    m_jac[0] = d1[0][0]*d2[1][0] - d2[0][0]*d1[1][0];

                    ASSERTL1(m_jac[0] > 0, "2D Regular Jacobian is not positive");

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
                m_gmat = Array<TwoD, NekDouble>(2*coordim,nqtot,0.0);

                if(coordim == 2) // assume g = [0,0,1]
                {
                    // set up Jacobian 
                    Vmath::Vmul (nqtot,&d2[0][0],1,&d1[1][0],1,&m_jac[0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&m_jac[0],1,
                                 &m_jac[0],1);

                    ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0, "2D Deformed Jacobian is not positive");
                    
                    Vmath::Vdiv(nqtot,&d2[1][0],1,&m_jac[0],1,&m_gmat[0][0],1); // d xi_1/d x_1
                    Vmath::Vdiv(nqtot,&d1[1][0],1,&m_jac[0],1,&m_gmat[1][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[1][0],1);                  // d xi_2/d x_1
                    Vmath::Vdiv(nqtot,&d2[0][0],1,&m_jac[0],1,&m_gmat[2][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[2][0],1);                  // d xi_1/d x_2
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
                    Vmath::Vvtvm(nqtot,&d1[2][0],1,&g[1][0],1,&m_gmat[1][0],1
                                 ,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);

                    // d xi_1/d x_2
                    Vmath::Vmul (nqtot,&d2[0][0],1,&g[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vvtvm(nqtot,&d2[2][0],1,&g[0][0],1,&m_gmat[2][0],1,
                                 &m_gmat[2][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);

                    // d xi_2/d x_2
                    Vmath::Vmul (nqtot,&d1[2][0],1,&g[0][0],1,&m_gmat[3][0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&g[2][0],1,&m_gmat[3][0],1,
                                 &m_gmat[3][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);

                    // d xi_1/d x_3
                    Vmath::Vmul (nqtot,&d2[1][0],1,&g[0][0],1,&m_gmat[4][0],1);
                    Vmath::Vvtvm(nqtot,&d2[0][0],1,&g[1][0],1,&m_gmat[4][0],1,
                                 &m_gmat[4][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);

                    // d xi_2/d x_3
                    Vmath::Vmul (nqtot,&d1[0][0],1,&g[1][0],1,&m_gmat[5][0],1);
                    Vmath::Vvtvm(nqtot,&d1[1][0],1,&g[0][0],1,&m_gmat[5][0],1,
                                 &m_gmat[5][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);

                    // J = sqrt(J_3D)
                    Vmath::Vsqrt(nqtot,&m_jac[0],1,&m_jac[0],1);

                }
            }
        }
#ifdef HIGH_D_FUNCTIONS


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

        GeomFactors::GeomFactors(const GeomType gtype, 
            const StdRegions::StdExpansion3D **Coords)  
        {

            int        i,nquad0,nquad1,nquad2,nqtot;
            StdRegions::PointsType  ptype0, ptype1, ptype2;
            double      *d1[3], *d2[3], *d3[3],*tmp;

            m_gtype = gtype;
            m_gmat  = new double*[9];

            nquad0 = Coords[0]->GetPointsOrder(0);
            nquad1 = Coords[0]->GetPointsOrder(1);
            nquad2 = Coords[0]->GetPointsOrder(2);
            ptype0 = Coords[0]->GetPointsType (0);
            ptype1 = Coords[0]->GetPointsType (1);
            ptype2 = Coords[0]->GetPointsType (2);

            nqtot = nquad0*nquad1*nquad2;

            d1[0] = new double [3*nqtot];
            d1[1] = d1[0] + nqtot;
            d1[2] = d1[1] + nqtot;

            d2[0] = new double [3*nqtot];
            d2[1] = d2[0] + nqtot;
            d2[2] = d2[1] + nqtot;

            d3[0] = new double [3*nqtot];
            d3[1] = d3[0] + nqtot;
            d3[2] = d3[1] + nqtot;

            // Calculate local derivatives using physical space storage
            for(i = 0; i < 3; ++i)
            {
                ASSERTL2(Coords[i]->GetPointsOrder(0) == nquad0,
                    "Points order are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsOrder(1) == nquad1,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsOrder(2) == nquad2,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(0)  == ptype0,
                    "Points type are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsType(1)  == ptype1,
                    "Points type are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(2)  == ptype2,
                    "Points type are different for coordinate 1 ");

                ((StdRegions::StdExpansion3D **) Coords)[i]->StdDeriv(d1[i],d2[i],d3[i]);
            }

            if((m_gtype == eRegular)||
               (m_gtype == eMovingRegular))
            {
                m_jac      = new double [1];
                m_gmat[0]  = new double [9];

                for(i = 1; i < 9; ++i)
                {
                    m_gmat[i] = m_gmat[i-1]+1;
                }

                m_jac[0] = d1[0][0]*(d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])
                    - d2[0][0]*(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0]) 
                    + d3[0][0]*(d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0]);

                // d xi_1/d x_1
                m_gmat[0][0] =  (d2[1][0]*d3[2][0] - d2[2][0]*d2[2][0])/m_jac[0]; 
                // d xi_2/d x_1
                m_gmat[1][0] = -(d1[1][0]*d3[2][0] - d1[2][0]*d3[2][0])/m_jac[0]; 
                // d xi_3/d x_1
                m_gmat[2][0] =  (d1[1][0]*d3[1][0] - d1[2][0]*d3[1][0])/m_jac[0]; 

                // d xi_1/d x_2
                m_gmat[3][0] = -(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0])/m_jac[0]; 
                // d xi_2/d x_2
                m_gmat[4][0] =  (d1[0][0]*d3[2][0] - d1[2][0]*d3[0][0])/m_jac[0]; 
                // d xi_3/d x_2
                m_gmat[5][0] = -(d1[0][0]*d2[2][0] - d1[2][0]*d2[0][0])/m_jac[0]; 

                // d xi_1/d x_3
                m_gmat[6][0] =  (d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0])/m_jac[0]; 
                // d xi_2/d x_3
                m_gmat[5][0] = -(d1[0][0]*d3[1][0] - d1[1][0]*d3[0][0])/m_jac[0]; 
                // d x7_3/d x_3
                m_gmat[8][0] =  (d1[0][0]*d2[1][0] - d1[1][0]*d2[0][0])/m_jac[0]; 

            }
            else
            {
                m_jac     = new double [nqtot];
                m_gmat[0] = new double [9*nqtot];

                for(i = 1; i < 9; ++i)
                {
                    m_gmat[i] = m_gmat[i-1] + nqtot;
                }

                tmp = new double [nqtot]; 

                // Jacobian
                Vmath::Vmul (nqtot,d3[1],1,d2[2],1,tmp,1);
                Vmath::Vvtvm(nqtot,d2[1],1,d3[2],1,tmp,1,tmp,1);
                Vmath::Vmul (nqtot,d1[0],1,tmp,1,m_jac,1);

                Vmath::Vmul (nqtot,d1[1],1,d3[2],1,tmp,1);
                Vmath::Vvtvm(nqtot,d3[1],1,d1[2],1,tmp,1,tmp,1);
                Vmath::Vvtvp(nqtot,d2[0],1,tmp,1,m_jac,1,m_jac,1);

                Vmath::Vmul (nqtot,d2[1],1,d1[2],1,tmp,1);
                Vmath::Vvtvm(nqtot,d1[1],1,d2[2],1,tmp,1,tmp,1);
                Vmath::Vvtvp(nqtot,d3[0],1,tmp,1,m_jac,1,m_jac,1);

                // d xi_1/d x_1
                Vmath::Vmul (nqtot,d2[2],1,d3[1],1,m_gmat[0],1);
                Vmath::Vvtvm(nqtot,d2[1],1,d3[2],1,m_gmat[0],1,m_gmat[0],1);
                Vmath::Vdiv(nqtot,m_gmat[0],1,m_jac,1,m_gmat[0],1);

                // d xi_2/d x_1
                Vmath::Vmul (nqtot,d1[1],1,d3[2],1,m_gmat[1],1);
                Vmath::Vvtvm(nqtot,d1[2],1,d3[1],1,m_gmat[1],1,m_gmat[1],1);
                Vmath::Vdiv(nqtot,m_gmat[1],1,m_jac,1,m_gmat[1],1);

                // d xi_3/d x_1
                Vmath::Vmul (nqtot,d1[1],1,d2[2],1,m_gmat[2],1);
                Vmath::Vvtvm(nqtot,d1[2],1,d3[1],1,m_gmat[2],1,m_gmat[2],1);
                Vmath::Vdiv(nqtot,m_gmat[2],1,m_jac,1,m_gmat[2],1);

                // d xi_1/d x_2
                Vmath::Vmul (nqtot,d2[0],1,d3[2],1,m_gmat[3],1);
                Vmath::Vvtvm(nqtot,d2[2],1,d3[0],1,m_gmat[3],1,m_gmat[3],1);
                Vmath::Vdiv(nqtot,m_gmat[3],1,m_jac,1,m_gmat[3],1);

                // d xi_2/d x_2
                Vmath::Vmul (nqtot,d1[2],1,d3[0],1,m_gmat[4],1);
                Vmath::Vvtvm(nqtot,d1[0],1,d3[2],1,m_gmat[4],1,m_gmat[4],1);
                Vmath::Vdiv(nqtot,m_gmat[4],1,m_jac,1,m_gmat[4],1);

                // d xi_3/d x_2
                Vmath::Vmul (nqtot,d1[0],1,d2[2],1,m_gmat[5],1);
                Vmath::Vvtvm(nqtot,d1[2],1,d3[0],1,m_gmat[5],1,m_gmat[5],1);
                Vmath::Vdiv(nqtot,m_gmat[5],1,m_jac,1,m_gmat[5],1);

                // d xi_1/d x_3
                Vmath::Vmul (nqtot,d2[1],1,d3[0],1,m_gmat[6],1);
                Vmath::Vvtvm(nqtot,d2[0],1,d3[1],1,m_gmat[6],1,m_gmat[6],1);
                Vmath::Vdiv(nqtot,m_gmat[6],1,m_jac,1,m_gmat[6],1);

                // d xi_2/d x_3
                Vmath::Vmul (nqtot,d1[0],1,d3[2],1,m_gmat[7],1);
                Vmath::Vvtvm(nqtot,d1[1],1,d3[0],1,m_gmat[7],1,m_gmat[7],1);
                Vmath::Vdiv(nqtot,m_gmat[7],1,m_jac,1,m_gmat[7],1);

                // d xim_3/d x_3
                Vmath::Vmul (nqtot,d1[1],1,d2[0],1,m_gmat[8],1);
                Vmath::Vvtvm(nqtot,d1[0],1,d2[1],1,m_gmat[8],1,m_gmat[8],1);
                Vmath::Vdiv(nqtot,m_gmat[8],1,m_jac,1,m_gmat[8],1);

                delete[] tmp;
            }

            delete [] d1[0];
            delete [] d2[0];
            delete [] d3[0];
        }

#endif

        GeomFactors::~GeomFactors(){
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: GeomFactors.cpp,v $
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
